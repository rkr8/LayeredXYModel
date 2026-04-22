function E_local(params::SimulationParams, ctx::SimulationContext, pos::Coord3D, v_prop::Spin2D)
    return sum(
        -params.J[i] * (
            v_prop ⋅ ctx.spins[pos+basis[i],V()] +
            v_prop ⋅ ctx.spins[pos-basis[i],V()]
        )
        for i = 1:3
    )
end

function E_total(params::SimulationParams, ctx::SimulationContext)
    if !ctx.spins.vectors_up_to_date[]
        error("Energy calculation requires vector representation of spins")
    end

    L_x::Int, L_y::Int, L_z::Int = size(ctx.spins)
    sum::Float64 = 0.0
    for x=1:L_x, y=1:L_y, z=1:L_z
        pos = Coord3D(x,y,z)
        for i = 1:3
            sum += -params.J[i] * (ctx.spins[pos,V()] ⋅ ctx.spins[pos+basis[i],V()])
        end
    end
    return sum
end

function local_field(params::SimulationParams, ctx::SimulationContext, pos::Coord3D)
    return sum(
        params.J[i] * (
            ctx.spins[pos+basis[i],V()] +
            ctx.spins[pos-basis[i],V()]
        )
        for i = 1:3
    )
end

function mag_layer(ctx::SimulationContext)
    if !ctx.spins.vectors_up_to_date[]
        error("Magnetization calculation requires vector representation of spins")
    end

    L_x::Int, L_y::Int, L_z::Int = size(ctx.spins)

    mag_x::Float64 = 0.0
    mag_y::Float64 = 0.0
    mag_sum::Float64 = 0.0
    for i = 1:L_z
        mag_x = mean(view(ctx.spins.vectors, :, :, i, 1))
        mag_y = mean(view(ctx.spins.vectors, :, :, i, 2))
        mag_sum += sqrt(mag_x^2 + mag_y^2)
    end

    return ( mag_sum / L_z)
end

function mag2_total(ctx::SimulationContext)
    if !ctx.spins.vectors_up_to_date[]
        error("Magnetization calculation requires vector representation of spins")
    end

    mag_x::Float64 = mean(view(ctx.spins.vectors, :, :, :, 1))
    mag_y::Float64 = mean(view(ctx.spins.vectors, :, :, :, 2))
    return mag_x^2 + mag_y^2
end

function mag_phase_diff(ctx::SimulationContext)
    if !ctx.spins.vectors_up_to_date[]
        error("Rel. phase calculation requires vector representation of spins")
    end

    L_x::Int, L_y::Int, L_z::Int = size(ctx.spins)

    scalar_prod_sum::Float64 = 0.0
    for i = 1:L_z
        scalar_prod_sum += (
            mean(view(ctx.spins.vectors, :, :, i, 1)) * mean(view(ctx.spins.vectors, :, :, mod1(i+1,L_z), 1)) +
            mean(view(ctx.spins.vectors, :, :, i, 2)) * mean(view(ctx.spins.vectors, :, :, mod1(i+1,L_z), 2))
        )
    end

    return ( scalar_prod_sum / L_z)
end

function J_d(params::SimulationParams, ctx::SimulationContext; axis::Int=1)
    L_x::Int, L_y::Int, L_z::Int = size(ctx.spins)
    sum::Float64 = 0.0
    for x=1:L_x, y=1:L_y, z=1:L_z
        pos = Coord3D(x,y,z)
        sum += ctx.spins[pos,V()] ⋅ ctx.spins[pos+basis[axis],V()]
    end
    return params.J[axis] * sum / (L_x*L_y*L_z)
end

# computes sin(θ_A - θ_B) for two Spin2D vectors A and B
sin_diff(A::Spin2D, B::Spin2D) = A[1]*B[2] - A[2]*B[1]

function J_p(params::SimulationParams, ctx::SimulationContext; axis::Int=1)
    L_x::Int, L_y::Int, L_z::Int = size(ctx.spins)
    sum::Float64 = 0.0
    for x=1:L_x, y=1:L_y, z=1:L_z
        pos = Coord3D(x,y,z)
        sum += sin_diff(ctx.spins[pos,V()],ctx.spins[pos+basis[axis],V()])
    end
    return params.J[axis]^2 * sum^2 * ctx.β[] / (L_x*L_y*L_z)
end
