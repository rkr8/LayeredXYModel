struct Lattice3D
    # Angle & cartesian representation of spins
    angles::Array{Float64, 3}
    vectors::Array{Float64, 4} # vectors[i,j,k,:] corresponds to vector at pos. i,j,k
    # Keep track of whether the angles or vectors are up to date
    angles_up_to_date::Base.RefValue{Bool}
    vectors_up_to_date::Base.RefValue{Bool}

    function Lattice3D(L_x::Int, L_y::Int, L_z::Int)
        angles = Array{Float64}(undef, L_x, L_y, L_z)
        vectors = Array{Float64}(undef, L_x, L_y, L_z, 2)
        angles_up_to_date = false
        vectors_up_to_date = false
        new(angles, vectors, Ref(angles_up_to_date), Ref(vectors_up_to_date))
    end
end

# randomly initialize lattice
function init_lattice!(lattice::Lattice3D)
    L_x, L_y, L_z = size(lattice.angles)
    lattice.angles .= 2π * rand(Float64, L_x, L_y, L_z)
    lattice.vectors[:, :, :, 1] .= cos.(lattice.angles)
    lattice.vectors[:, :, :, 2] .= sin.(lattice.angles)
    lattice.angles_up_to_date[] = true
    lattice.vectors_up_to_date[] = true
end

# initialize lattice with given angles
function init_lattice!(lattice::Lattice3D, angles::Array{Float64, 3})
    lattice.angles .= angles
    lattice.vectors[:, :, :, 1] .= cos.(lattice.angles)
    lattice.vectors[:, :, :, 2] .= sin.(lattice.angles)
    lattice.angles_up_to_date[] = true
    lattice.vectors_up_to_date[] = true
end

# Two singleton types for parametrization, to use multiple dispatch
struct A end # angle parametrization
struct V end # vector parametrization

const Coord3D = SVector{3, Int64}
const Couplings3D = SVector{3, Float64}
const Spin2D = SVector{2, Float64}
const x̂ = Coord3D(1,0,0)
const ŷ = Coord3D(0,1,0)
const ẑ = Coord3D(0,0,1)
const basis = (x̂,ŷ,ẑ)
const nn = (x̂,ŷ,ẑ,-x̂,-ŷ,-ẑ)

# Define size function for custom lattice
function Base.size(lattice::Lattice3D)
    return size(lattice.angles)
end

function Base.getindex(spins::Lattice3D, pos::Coord3D, ::A) # get angle
    L_x::Int, L_y::Int, L_z::Int = size(spins)
    x, y, z = pos
    return spins.angles[mod1(x,L_x), mod1(y,L_y), mod1(z,L_z)]
end

# TODO: maybe make this the default behavior
function Base.getindex(spins::Lattice3D, pos::Coord3D, ::V) # get vector
    L_x::Int, L_y::Int, L_z::Int = size(spins)
    x, y, z = pos
    #return view(spins.vectors, mod1(x,L_x), mod1(y,L_y), mod1(z,L_z), :)
    return Spin2D(
        spins.vectors[mod1(x,L_x), mod1(y,L_y), mod1(z,L_z), 1],
        spins.vectors[mod1(x,L_x), mod1(y,L_y), mod1(z,L_z), 2]
    )
end

function Base.setindex!(spins::Lattice3D, φ::Float64, pos::Coord3D, ::A) # set angle
    L_x::Int, L_y::Int, L_z::Int = size(spins)
    x, y, z = pos
    spins.angles[mod1(x,L_x), mod1(y,L_y), mod1(z,L_z)] = φ
    spins.vectors_up_to_date[] = false
end

# TODO: maybe make this the default behavior
function Base.setindex!(spins::Lattice3D, vec::Spin2D, pos::Coord3D, ::V) # set vector
    L_x::Int, L_y::Int, L_z::Int = size(spins)
    x, y, z = pos
    spins.vectors[mod1(x,L_x), mod1(y,L_y), mod1(z,L_z), :] .= vec
    spins.angles_up_to_date[] = false
end

# vectors -> angles
function update!(spins::Lattice3D, ::A)
    if !spins.angles_up_to_date[]
        spins.angles[:, :, :] .= atan.(
            view(spins.vectors, :, :, :, 2),
            view(spins.vectors, :, :, :, 1)
        )
        spins.angles_up_to_date[] = true
    end
    nothing
end

# angles -> vectors
function update!(spins::Lattice3D, ::V)
    if !spins.vectors_up_to_date[]
        spins.vectors[:, :, :, 1] .= cos.(spins.angles)
        spins.vectors[:, :, :, 2] .= sin.(spins.angles)
        spins.vectors_up_to_date[] = true
    end
    nothing
end

# Converts an sublattice index to a full lattice index
function to_lattice_ind(sublattice_pos::Coord3D, black::Bool)
    sl_x, sl_y, sl_z = sublattice_pos
    return Coord3D(
        sl_x,
        sl_y,
        2 * sl_z - (black ? (sl_x+sl_y+1) : (sl_x+sl_y))%2
    )
end

function shift_coordinates(pos::Coord3D, dim::Coord3D)
    x::Int, y::Int, z::Int = pos
    L_x::Int, L_y::Int, L_z::Int = dim
    return Coord3D(mod1(x,L_x), mod1(y,L_y), mod1(z,L_z))
end

# Static things
@kwdef struct SimulationParams
    relax_steps::Int = 100_000 # max. number of thermalization steps, corresponding to 15 bins
    sim_steps::Int = 200_000
    save_interval::Int = 1
    log_interval::Int = 8_000
    swap_interval::Int
    swap_measurements::Bool
    config_save_interval::Int
    J::Couplings3D = Couplings3D(1.0, 1.0, 1.0)
    βmin::Float64
    βmax::Float64
    nonlinear::Bool
    L::Int
    save_dir::String
    init_config_path::Union{String, Nothing} = nothing
    output_config_path::Union{String, Nothing} = nothing
    observables::Vector{Symbol}
end

# Dynamical things
struct SimulationContext
    β::Base.RefValue{Float64}
    last_log::Base.RefValue{Float64}
    cluster_size_sum::Base.RefValue{Int}
    swaps::Base.RefValue{Int}
    file::Base.RefValue{IOStream}
    spins::Lattice3D
    cluster::Set{Coord3D}
    to_check::Set{Coord3D}
    measurement::Vector{Float64}

    function SimulationContext(β::Float64, L::Int)
        new(
            Ref(β),
            Ref(time()),
            Ref(0),
            Ref(0),
            Ref{IOStream}(),
            Lattice3D(L, L, L),
            Set{Coord3D}(),
            Set{Coord3D}(),
            Float64[]
        )
    end
end