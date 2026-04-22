function save_observables!(params::SimulationParams, ctx::SimulationContext)
    empty!(ctx.measurement) # re-use ctx.measurement array to prevent inefficient re-allocations
    for o = params.observables
        if o == :mpi_rank # useful for debugging parallel tempering
            push!(ctx.measurement, MPI_RANK())
        elseif o == :energy # measure energy (extensive)
            push!(ctx.measurement, E_total(params, ctx))
        elseif o == :mag2 # measure magnetizaton^2 of entire system (intensive)
            push!(ctx.measurement, mag2_total(ctx))
        elseif o == :mag_layer # average absolute single layer magnetizaton (intensive)
            push!(ctx.measurement, mag_layer(ctx))
        elseif o == :mag_phase_diff # measure scalar of magnetizations of adjacent layers
            push!(ctx.measurement, mag_phase_diff(ctx))
        elseif o == :Jdx # measure diamagnetic term of superfluid stiffness (along x-axis)
            push!(ctx.measurement, J_d(params, ctx; axis=1))
        elseif o == :Jdz # ... along z-axis
            push!(ctx.measurement, J_d(params, ctx; axis=3))
        elseif o == :Jpx # measure paramagnetic term of superfluid stiffness (along x-axis)
            push!(ctx.measurement, J_p(params, ctx; axis=1))
        elseif o == :Jpz # ... along z-axis
            push!(ctx.measurement, J_p(params, ctx; axis=3))
        else
            error("Encountered unknown observable: $o, stopping simulation...")
        end
    end
    # write to disk
    serialize(
        ctx.file[],
        ctx.measurement
    )
end

# swapping probability for the parallel tempering
function swap_probability(E_i::Float64, E_j::Float64, β_i::Float64, β_j::Float64)
    return min(1.0, exp( (E_i - E_j)*(β_i - β_j) ))
end

function parallel_tempering!(params::SimulationParams, ctx::SimulationContext, iteration::Int; reopen_file::Bool=false)
    local partner_energy::Float64
    local partner_β::Float64
    local p::Float64
    local swap::Bool
    partner_rank::Int = get_partner_rank(iteration)

    if 0 <= partner_rank < MPI_SIZE() # open BC’s
        # smaller MPI rank makes swap decision
        if MPI_RANK() < partner_rank
            partner_energy, partner_β = MPI.recv(MPI_COMM; source=partner_rank)
            p = swap_probability(E_total(params, ctx), partner_energy, ctx.β[], partner_β)
            swap = ( rand() <= p )
            MPI.send((swap, ctx.β[]), MPI_COMM; dest=partner_rank)

            # DEBUG message for parallel tempering
            println("iteration $iteration, β: $(round(ctx.β[]; digits=3))<->$(round(partner_β; digits=3)), E: $(round(E_total(params, ctx); digits=3))<->$(round(partner_energy; digits=3)) with p = $(round(p; digits=3))");
        else # larger MPI rank follows
            MPI.send((E_total(params, ctx), ctx.β[]), MPI_COMM; dest=partner_rank)
            swap, partner_β = MPI.recv(MPI_COMM; source=partner_rank)
        end

        if swap
            # Parallel tempering happens only during thermalization!
            ctx.β[] = partner_β
            ctx.swaps[] += 1
            if reopen_file
                # close the file and reopen it to append new data
                close(ctx.file[])
                sync_with_partner(partner_rank)
                ctx.file[] = open(joinpath(params.save_dir, "$(ctx.β[]).dat"), "a")
            end
        end
    end
end

# main entry point to run MC simulation
function run_simulation(params::SimulationParams)
    # Print MPI debug info
    (MPI_RANK() == MPI_ROOT) && MPI.versioninfo()

    printroot("Running simulation on $(MPI_SIZE()) ranks!")
    printroot("Observables: $(params.observables)")

    # Contains all dynamical things, different for each rank
    β_init::Float64 = 0.0
    if params.nonlinear
        # Linear spacing of temperatures (T)
        β_init = params.βmin * params.βmax / ( params.βmin + (params.βmax - params.βmin) * (MPI_RANK() / (MPI_SIZE() - 1)) )
    else
        # Linear spacing of inverse temperatures (β)
        β_init = params.βmin + (params.βmax - params.βmin) * (MPI_RANK() / (MPI_SIZE() - 1))
    end
    
    ctx = SimulationContext(
        β_init,
        params.L
    )
    println("Initial value: β = $(ctx.β[])")

    if isnothing(params.init_config_path)
        println("Initializing spins randomly")
        init_lattice!(ctx.spins)
    else
        println("Loading spins from file: $(params.init_config_path)")
        init_lattice!(ctx.spins, deserialize(params.init_config_path))
    end

    MPI.Barrier(MPI_COMM)
    printroot("Running $(params.relax_steps) thermalization steps with parallel tempering")
    for i = 1:params.relax_steps
        # Propose parallel tempering swap at every swap_interval steps
        if i % params.swap_interval == 0
            parallel_tempering!(params, ctx, i ÷ params.swap_interval; reopen_file=false)
        end
        update_heatbath!(params, ctx)
    end
    MPI.Barrier(MPI_COMM)
    println("$(ctx.swaps[]) swaps")
    #println("Finished thermalization in $(time()-ctx.last_log[]) seconds, $(ctx.swaps[]) swaps")
    #println("Final value: β = $(ctx.β[])")

    # if params.save_dir does not exist, create it
    !isdir(params.save_dir) && mkpath(params.save_dir)
    # if params.output_config_path does not exist, create it
    !isnothing(params.output_config_path) && !isdir(params.output_config_path) && mkpath(params.output_config_path)
    ctx.file[] = open(joinpath(params.save_dir, "$(ctx.β[]).dat"), "w")

    serialize(ctx.file[], (params.L, ctx.β[]))
    serialize(ctx.file[], params.observables)

    ctx.last_log[] = time()
    for i = 1:params.sim_steps
        if i % params.log_interval == 0
            println("Progress = $(i / params.sim_steps), avg. Wolff cluster size = $(ctx.cluster_size_sum[]/params.log_interval), speed = $(params.log_interval/(time()-ctx.last_log[])) steps/second")
            ctx.last_log[] = time()
            ctx.cluster_size_sum[] = 0
        end

        if !isnothing(params.output_config_path) && (i % params.config_save_interval == 0)
            spin_path = joinpath(params.output_config_path, "config_L$(params.L)_BETA$(ctx.β[])_STEP$i.dat")
            println("Saving spins to file: $spin_path")
            update!(ctx.spins, A())
            serialize(spin_path, ctx.spins.angles)
        end

        if params.swap_measurements && (i % params.swap_interval == 0)
            parallel_tempering!(params, ctx, i ÷ params.swap_interval; reopen_file=true)
        end

        # local updates
        update_heatbath!(params, ctx)
        update_wolff!(params, ctx)
        ctx.cluster_size_sum[] += length(ctx.cluster)
        
        # these are not so efficient:
        #update_metropolis!(params, ctx)
        #update_over_relaxation!(params, ctx)

        if i % params.save_interval == 0
            save_observables!(params, ctx)
        end
    end

    # close the file
    close(ctx.file[])
    println("Finished simulation")

    nothing
end

;
