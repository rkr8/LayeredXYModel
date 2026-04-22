module xy_simulation

# external libraries
using MPI
using StaticArrays
using Random
using LinearAlgebra
using Statistics
using Serialization
using Dates
using UUIDs

export MPI

# source files
include("mpi.jl")
include("types.jl")
include("observables.jl")
include("von_mises.jl")
include("mc_algorithms.jl")
include("simulation_logic.jl")

# make internal components accessible
export Couplings3D, SimulationParams, run_simulation

end