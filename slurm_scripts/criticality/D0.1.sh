#!/bin/bash

DELTA=0.1
L_VALS=(14 16 18 20 22 24 28 32 36 40 44 48 56 64 72)
BETA_C=0.7575
NU=0.67183
THM_STEPS=400000
SIM_STEPS=800000
SWAP_INTERVAL=100
OBSERVABLES="energy mag2 Jdx Jpx Jdz Jpz"
NUM_BETAS=48
NUM_CONFIGS=500

SAVE_DIR="$SCRATCH/final_runs_measurements/$DELTA"
CONFIG_DIR="$SCRATCH/final_runs_configs/$DELTA"
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

for L in ${L_VALS[@]}; do
    # Scale beta spacing for curve collapse
    BETA_MIN=$(echo "$BETA_C - 0.7 * e(l($L) * (-1.0 / $NU))" | bc -l)
    BETA_MAX=$(echo "$BETA_C + 0.7 * e(l($L) * (-1.0 / $NU))" | bc -l)

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