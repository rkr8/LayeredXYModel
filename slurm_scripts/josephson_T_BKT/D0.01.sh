#!/bin/bash

DELTA=0.01
L_MIN=2
L_MAX=50

# T_BKT=0.8929
# T_MAX=0.937545
# T_MIN=0.848255

BETA_MIN=1.066615
BETA_BKT=1.1199462426
BETA_MAX=1.178891
NUM_BETAS=8

THM_STEPS=400000
SIM_STEPS=600000
SWAP_INTERVAL=100
NUM_CONFIGS=500
OBSERVABLES="energy mag2 mag_layer mag_phase_diff Jdx Jpx Jdz Jpz"

SAVE_DIR="$SCRATCH/final_runs_measurements_josephson_T_BKT/$DELTA"
CONFIG_DIR="$SCRATCH/final_runs_configs_josephson_T_BKT/$DELTA"
PROJECT_DIR="$HOME/xy_simulation"
MAIN_SCRIPT="$PROJECT_DIR/main.jl"
LOG_DIR="$PROJECT_DIR/logs/$(date +"%Y-%m-%d_%H-%M-%S")"
mkdir -p $LOG_DIR

module load stack/2024-06 openmpi julia

julia --optimize=3 \
      --min-optlevel=3 \
      --check-bounds=no \
      --project=$PROJECT_DIR \
      -e "using Pkg; Pkg.instantiate(); Pkg.precompile(); using xy_simulation; println(\"Precompilation done.\")"

for L in $(seq $L_MIN 2 $L_MAX); do
    JOB_NAME="D${DELTA}_L${L}_B${BETA}_R${RUN}_PT"
    LOG_FILE="$LOG_DIR/$JOB_NAME.log"
    # Calculate max. expected runtime in seconds
    TOTAL_STEPS=$(($THM_STEPS + $SIM_STEPS))
    # 120 seconds minimum runtime for setup
    TIME_IN_SECONDS=$(echo "$TOTAL_STEPS*0.000001*$L*$L*$L + 120" | bc -l)
    TIME_IN_SECONDS=${TIME_IN_SECONDS%.*}
    # 4. Break into HH:MM:SS
    HOURS=$(( TIME_IN_SECONDS/3600 ))
    MINS=$(( (TIME_IN_SECONDS%3600)/60 ))
    SECS=$(( TIME_IN_SECONDS%60 ))
    TIME_STR=$(printf "%02d:%02d:%02d" "$HOURS" "$MINS" "$SECS")

    julia_command="julia --optimize=3 --min-optlevel=3 --check-bounds=no --project=$PROJECT_DIR $MAIN_SCRIPT \
                            --thermalization $THM_STEPS \
                            --simulation-steps $SIM_STEPS \
                            --swap-interval $SWAP_INTERVAL \
                            --swap-measurements false \
                            --config-num $NUM_CONFIGS \
                            --delta $DELTA \
                            --beta-min $BETA_MIN \
                            --beta-max $BETA_MAX \
                            --nonlinear true \
                            --lattice-size $L \
                            --save-dir "$SAVE_DIR/$L" \
                            --output-config $CONFIG_DIR \
                            --observables $OBSERVABLES"
    sbatch --ntasks=$NUM_BETAS \
            --cpus-per-task=1 \
            --time=$TIME_STR \
            --mem-per-cpu=2048 \
            --mail-type=NONE \
            --job-name=$JOB_NAME \
            --output=$LOG_FILE \
            --wrap="export JULIA_NUM_THREADS=1; export OPENBLAS_NUM_THREADS=1; mpirun $julia_command"
done