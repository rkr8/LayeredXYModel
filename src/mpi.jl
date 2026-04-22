# MPI-related constants
const MPI_ROOT::Int = 0
const MPI_COMM::MPI.Comm = MPI.COMM_WORLD
const MPI_RANK() = MPI.Comm_rank(MPI_COMM)
const MPI_SIZE() = MPI.Comm_size(MPI_COMM)

println(args...) = Base.println("Rank $(MPI_RANK()): ", args...)
printroot(args...) = (MPI_RANK() == MPI_ROOT) && println(args...)

# partner rank for parallel tempering swaps
function get_partner_rank(iteration::Int)
    partner_rank = if isodd(iteration)
        (iseven(MPI_RANK()) ? MPI_RANK() + 1 : MPI_RANK() - 1)
    else
        (iseven(MPI_RANK()) ? MPI_RANK() - 1 : MPI_RANK() + 1)
    end
    #return mod(partner_rank, MPI_SIZE()) # periodic boundary conditions
    return partner_rank # open boundary conditions
end

# like MPI.Barrier, but only for 2 ranks
function sync_with_partner(partner_rank::Int)
    MPI.send(0, MPI_COMM; dest=partner_rank)
    MPI.recv(MPI_COMM; source=partner_rank)
    nothing
end