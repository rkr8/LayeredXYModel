# metropolis sweep for black or white checkerboard sublattice
function update_metropolis!(params::SimulationParams,
                            ctx::SimulationContext,
                            update_black::Bool)
    if !ctx.spins.vectors_up_to_date[]
        error("Metropolis algorithm requires vector representation of spins")
    end

    L_x::Int, L_y::Int, L_z::Int = size(ctx.spins)

    for x=1:L_x, y=1:L_y, z=1:L_z
        if update_black ⊻ iseven(x+y+z)
            pos::Coord3D = Coord3D(x,y,z)
            φ::Float64 = 2π * rand(Float64)
            v::Spin2D = Spin2D(cos(φ), sin(φ))
            ΔE::Float64 = E_local(params, ctx, pos, v) - E_local(params, ctx, pos, ctx.spins[pos,V()])
            if ΔE < 0.0
                ctx.spins[pos,V()] = v
            elseif rand() <= exp(-ctx.β[] * ΔE)
                ctx.spins[pos,V()] = v
            end
        end
    end
    nothing
end

# full metropolis sweep
function update_metropolis!(params::SimulationParams, ctx::SimulationContext)
    update_metropolis!(params, ctx, true) # update black
    update_metropolis!(params, ctx, false) # update white
    nothing
end

# heatbath sweep for black or white checkerboard sublattice
function update_heatbath!(params::SimulationParams,
                          ctx::SimulationContext,
                          update_black::Bool)
    if !ctx.spins.vectors_up_to_date[]
        error("Heatbath algorithm requires vector representation of spins")
    end

    L_x::Int, L_y::Int, L_z::Int = size(ctx.spins)
    for x=1:L_x, y=1:L_y, z=1:L_z
        if update_black ⊻ iseven(x+y+z)
            pos::Coord3D = Coord3D(x,y,z)
            h::Spin2D = local_field(params, ctx, pos)
            h_norm::Float64 = norm(h)
            h_angle::Float64 = atan(h[2], h[1])
            φ::Float64 = rand_von_mises(h_angle, ctx.β[] * h_norm)
            ctx.spins[pos,V()] = Spin2D(cos(φ), sin(φ))
        end
    end

    nothing
end

# full heatbath sweep
function update_heatbath!(params::SimulationParams, ctx::SimulationContext)
    update_heatbath!(params, ctx, true) # update black
    update_heatbath!(params, ctx, false) # update white
end

# bond prob. for cluster formation
function bond_prob(spin_a::Spin2D, spin_b::Spin2D, spin_ref::Spin2D, β::Float64, J::Float64)
    1.0 - min(1.0, exp(-2.0 * β * J * (spin_ref ⋅ spin_a) * (spin_ref ⋅ spin_b)))
end

# grow cluster from random point and flip all spins accordingly
function update_wolff!(params::SimulationParams, ctx::SimulationContext)
    if !ctx.spins.vectors_up_to_date[]
        error("Wolff algorithm requires vector representation of spins")
    end

    L_x::Int, L_y::Int, L_z::Int = size(ctx.spins)
    dim = Coord3D(L_x, L_y, L_z)

    # Clean up
    empty!(ctx.cluster)
    empty!(ctx.to_check)

    # Initialize
    initial_pos = Coord3D(rand(1:L_x), rand(1:L_y), rand(1:L_z))
    push!(ctx.cluster, initial_pos)
    push!(ctx.to_check, initial_pos)
    φ_r = 2 * π * rand(Float64)
    ref = Spin2D(cos(φ_r), sin(φ_r))

    while !isempty(ctx.to_check) # grow cluster
        # very important: pos = pop!(ctx.to_check) is super slow!
        # for some reason this is 2 orders of magnitude faster than pop!(ctx.to_check):
        pos = rand(ctx.to_check)
        delete!(ctx.to_check, pos)
        for i = 1:3 # check nn’s
            for dir = -1:2:1
                next = shift_coordinates(pos + dir*basis[i], dim) # periodic boundary conditions
                if next ∉ ctx.cluster
                    if (rand(Float64) <= bond_prob(ctx.spins[pos,V()], ctx.spins[next,V()], ref, ctx.β[], params.J[i]))
                        push!(ctx.cluster, next)
                        push!(ctx.to_check, next)
                    end
                end
            end
        end
    end

    # reflect spins in cluster
    for pos = ctx.cluster
        ctx.spins[pos,V()] = ctx.spins[pos,V()] .- 2.0 * (ctx.spins[pos,V()] ⋅ ref) .* ref
    end

    nothing
end

# full over-relaxation sweep
function update_over_relaxation!(params::SimulationParams, ctx::SimulationContext)
    if !ctx.spins.vectors_up_to_date[]
        error("Over-relaxation algorithm requires vector representation of spins")
    end

    L_x::Int, L_y::Int, L_z::Int = size(ctx.spins)
    for x=1:L_x, y=1:L_y, z=1:L_z # reflect spins
        pos = Coord3D(x,y,z)
        v::Spin2D = ctx.spins[pos,V()]
        h::Spin2D = local_field(params, ctx, pos)
        # reflect:
        ctx.spins[pos,V()] = ( 2.0 * (v ⋅ h) / (h ⋅ h) ) .* h .- v
    end
end
