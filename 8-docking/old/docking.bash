#!/bin/bash
# =============================================================================
# AutoDock Vina - Batch Docking
# =============================================================================

# Paths
VINA_BIN="$HOME/Downloads/autodock_vina_1_1_2_linux_x86/bin/vina"
WORK_DIR="$HOME/Documents/vitoria/usp-masters/8-docking"
CONFIG_FILE="$2"
LIGANDS_DIR="$HOME/Documents/vitoria/usp-masters/6-preparation/ligands/pdbqt"

RECEPTOR="$1"
RECEPTOR_basename="${RECEPTOR%.*}"

# Outputs
OUT_DIR="$WORK_DIR/results_docking/$RECEPTOR_basename"
PROGRESS_LOG="$OUT_DIR/log.out"

# =============================================================================

# Aux functions
log_progress() {
    local ligand="$1"
    local status="$2"   # OK or FAIL
    local timestamp
    local seed="$3"

    timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "$timestamp | $status | $seed | $ligand" >> "$PROGRESS_LOG"
}
already_done() {
    local ligand="$1"
    grep -q "| OK | $ligand$" "$PROGRESS_LOG"
}

# Get all .pdbqt except for the protein
mapfile -t LIGANDS < <(find "$LIGANDS_DIR" -maxdepth 1 -name "*.pdbqt")

TOTAL=${#LIGANDS[@]}

# How many were done
DONE=0
for L in "${LIGANDS[@]}"; do
    already_done "$(basename "$L")" && ((DONE++))
done

echo "========================================"
echo " AutoDock Vina"
echo " Protein : $RECEPTOR_basename"
echo " Config   : $CONFIG_FILE"
echo " Ligands  : $TOTAL found | $DONE done"
echo " Progress : $PROGRESS_LOG"
echo " Output   : $OUT_DIR"
echo "========================================"

SUCCESS=0
FAIL=0
SKIP=0

N_REPS=3
for LIGAND_PATH in "${LIGANDS[@]}"; do
    LIGAND_FILE=$(basename "$LIGAND_PATH")
    LIGAND_NAME="${LIGAND_FILE%.pdbqt}"

    OUT_DIR_LIG="$OUT_DIR/$LIGAND_NAME"
    mkdir -p "$OUT_DIR_LIG"

    for REP in $(seq 1 $N_REPS); do

        SEED=$RANDOM

        TAG="${LIGAND_FILE}|rep${REP}"

        if already_done "$TAG"; then
            echo "Skipping: $LIGAND_NAME rep${REP}"
            ((SKIP++))
            continue
        fi

        OUT_FILE="$OUT_DIR_LIG/${LIGAND_NAME}_rep${REP}_out.pdbqt"
        LOG_FILE="$OUT_DIR_LIG/${LIGAND_NAME}_rep${REP}.log"

        echo ""
        echo ">>> Running: $LIGAND_NAME (rep ${REP})"

        "$VINA_BIN" \
            --config "$CONFIG_FILE" \
            --receptor "$RECEPTOR" \
            --ligand "$LIGAND_PATH" \
            --seed "$SEED" \
            --out "$OUT_FILE" \
            --log "$LOG_FILE"

        EXIT_CODE=$?

        if [[ $EXIT_CODE -eq 0 ]]; then
            log_progress "$TAG" "OK" "$SEED"
        else
            log_progress "$TAG" "FAIL" "$SEED"
        fi
    done
done

echo ""
echo "========================================"
echo " Done"
printf  "   Success: %d\n" "$SUCCESS"
printf  "   Fail   : %d\n" "$FAIL"
printf  "   SKipped: %d\n" "$SKIP"
echo " Results in: $OUT_DIR"
echo "========================================"

# Exit with error
[[ $FAIL -gt 0 ]] && exit 1
exit 0
