# Run the Monte Carlo simulation with parameters from the command line
using ArgParse
using xy_simulation

# For some reason console output isn’t line-buffered, when running
# the MC simulation on slurm. This is a work-around, to manually
# flush the stdout buffer after every line...
function Base.println(io::IO, xs...)
    print(io, xs..., '\n')
    flush(io)
end

# some example parameters:
# p = SimulationParams(
#     relax_steps = 0,
#     sim_steps = 30,
#     save_interval = 1,
#     log_interval = 1,
#     J = Couplings3D(1.0, 1.0, 1.0),
#     β = 0.45416,
#     L = 100,
#     save_dir = "./tmp",
#     #output_config_path = "./tmp/config_64.dat",
#     observables = [:mag2, :mag2_sl, :Jdx, :Jpx, :Jdz, :Jpz]
# )
# run_simulation(p)

function main()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--thermalization", "-t"
            help = "Number of thermalization steps"
            arg_type = Int
            default = 1_000
        "--simulation-steps", "-S"
            help = "Number of simulation steps"
            arg_type = Int
            default = 10_000
        "--save-interval", "-s"
            help = "Interval between saving observables"
            arg_type = Int
            default = 1
        "--log-interval", "-l"
            help = "Interval between logging progress"
            arg_type = Int
            default = 50
        "--swap-interval", "-w"
            help = "Interval between parallel tempering swaps"
            arg_type = Int
            default = 10
        "--swap-measurements", "-m"
            help = "Perform parallel tempering swaps also during measurements"
            arg_type = Bool
            default = false
        "--config-num", "-c"
            help = "Number of configurations to save (if output path is given)"
            arg_type = Int
            default = 10
        "--delta", "-D"
            help = "Anisotropy J_perp / J_parallel... 1 for isotropic"
            arg_type = Float64
            default = 1.0
        "--beta-min", "-b"
            help = "Minimum inverse temperature for parallel tempering"
            arg_type = Float64
            required = true
        "--beta-max", "-B"
            help = "Maximum inverse temperature for parallel tempering"
            arg_type = Float64
            required = true
        "--nonlinear", "-n"
            help = "If false, β is linearly spaced between βmin and βmax. If true, T is spaced linearly between Tmin=1/βmax and Tmax=1/βmin."
            arg_type = Bool
            default = false
        "--lattice-size", "-L"
            help = "Size of the lattice"
            arg_type = Int
            required = true
        "--save-dir", "-d"
            help = "Directory to save observables"
            arg_type = String
            required = true
        "--init-config", "-i"
            help = "Path to load the initial configuration"
            arg_type = String
            default = nothing
        "--output-config", "-o"
            help = "Path to save the final configuration"
            arg_type = String
            default = nothing
        "--observables", "-O"
            help = "Observables to save... select from: energy, mag2, mag_layer, mag_phase_diff, Jdx, Jdz, Jpx, Jpz"
            arg_type = Symbol
            nargs = '+'
            required = true
    end
    parsed_args = parse_args(s)

    # boilerplate... we could create the kwstruct directly from
    # the parsed_args dictionary, but the names have to match
    # TODO: include Δ
    p = SimulationParams(
        relax_steps = parsed_args["thermalization"],
        sim_steps = parsed_args["simulation-steps"],
        save_interval = parsed_args["save-interval"],
        log_interval = parsed_args["log-interval"],
        swap_interval = parsed_args["swap-interval"],
        swap_measurements = parsed_args["swap-measurements"],
        config_save_interval = parsed_args["simulation-steps"] ÷ parsed_args["config-num"],
        J = Couplings3D(1.0, 1.0, parsed_args["delta"]),
        βmin = parsed_args["beta-min"],
        βmax = parsed_args["beta-max"],
        nonlinear = parsed_args["nonlinear"],
        L = parsed_args["lattice-size"],
        save_dir = parsed_args["save-dir"],
        output_config_path = parsed_args["output-config"],
        init_config_path = parsed_args["init-config"],
        observables = parsed_args["observables"]
    )

    println("Initializing MPI...")
    MPI.Init()
    run_simulation(p)
    MPI.Finalize()
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

;