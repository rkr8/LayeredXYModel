# Layered XY Model

Monte Carlo simulation code for [arXiv:2603.19351](https://arxiv.org/abs/2603.19351).

## Installation

1. Install [Julia](https://julialang.org/downloads/).
2. Install OpenMPI on your system (typically preinstalled on HPC clusters).
3. Install Julia dependencies:
   ```bash
   julia --project=. -e "using Pkg; Pkg.instantiate();"
   ```

## Usage

View available options:
```bash
julia --project=. main.jl -h
```

Run a simulation:
```bash
mpiexec -n <num_ranks> julia --project=. main.jl [options]
```

See `./slurm_scripts` for example SLURM scripts and the parameters used in arXiv:2603.19351.