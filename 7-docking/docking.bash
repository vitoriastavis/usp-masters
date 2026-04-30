#!/bin/bash
# =============================================================================
# AutoDock Vina - Batch Docking 
# =============================================================================

# Paths
VINA_BIN="$HOME/Downloads/autodock_vina_1_1_2_linux_x86/bin/vina"
WORK_DIR="$HOME/Documents/vitoria/usp-masters/7-docking"   
CONFIG_FILE="$WORK_DIR/config.txt"

# Outputs
OUT_DIR="$WORK_DIR/results"
PROGRESS_LOG="$WORK_DIR/progress.log" 

# =============================================================================

# Aux functions
log_progress() {
    local ligand="$1"
    local status="$2"   # OK ou FALHA
    local timestamp
    timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "$timestamp | $status | $ligand" >> "$PROGRESS_LOG"
}
already_done() {
    local ligand="$1"
    grep -q "| OK | $ligand$" "$PROGRESS_LOG" 2>/dev/null
}

# Finds receptor file
RECEPTOR=$(grep -E '^\-\-receptor' "$CONFIG_FILE" | awk '{print $2}')

# Get all .pdbqt except for the protein
mapfile -t LIGANDS < <(find "$WORK_DIR" -maxdepth 1 -name "*.pdbqt" ! -name "$RECEPTOR" | sort)
TOTAL=${#LIGANDS[@]}

# How many were done
DONE=0
for L in "${LIGANDS[@]}"; do
    already_done "$(basename "$L")" && ((DONE++))
done

echo "========================================"
echo " AutoDock Vina"
echo " Protein : $RECEPTOR"
echo " Config   : $CONFIG_FILE"
echo " Ligands  : $TOTAL found | $DONE done"
echo " Progress : $PROGRESS_LOG"
echo " Output   : $OUT_DIR"
echo "========================================"

SUCESSO=0
FALHA=0
PULADO=0

for LIGAND_PATH in "${LIGANDS[@]}"; do
    LIGAND_FILE=$(basename "$LIGAND_PATH")
    LIGAND_NAME="${LIGAND_FILE%.pdbqt}"

    OUT_DIR_LIG="$OUT_DIR/$LIGAND_NAME/"
    mkdir -p OUT_DIR_LIG
    OUT_FILE="$OUTPUT_DIR/${LIGAND_NAME}_out.pdbqt"
    LOG_FILE="$OUTPUT_DIR/${LIGAND_NAME}.log"

    # Skip if it was already executed
    if already_done "$LIGAND_FILE"; then
        echo ""
        echo "Skipping: $LIGAND_NAME"
        ((PULADO++))
        continue
    fi

    echo ""
    echo ">>> [$((PULADO + SUCESSO + FALHA + 1))/$TOTAL] Running: $LIGAND_NAME"

    # --ligand, --out e --log sobrescrevem o que estiver no config.txt
    "$VINA_BIN" \
        --config  "$CONFIG_FILE" \
        --ligand  "$LIGAND_PATH" \
        --out     "$OUT_FILE" \
        --log     "$LOG_FILE" 

    EXIT_CODE=$?

    if [[ $EXIT_CODE -eq 0 ]]; then
        echo "    [OK] Save in: $OUT_FILE"
        log_progress "$LIGAND_FILE" "OK"
        ((SUCESSO++))
    else
        echo "    [FALHA] Exit code: $EXIT_CODE — check log: $LOG_FILE"
        log_progress "$LIGAND_FILE" "FALHA"
        ((FALHA++))
    fi
done

echo ""
echo "========================================"
echo " Done"
printf  "   Sucesso : %d\n" "$SUCESSO"
printf  "   Falha   : %d\n" "$FALHA"
printf  "   Pulados : %d\n" "$PULADO"
echo " Resultados em: $OUTPUT_DIR"
echo "========================================"

# Exit with error
[[ $FALHA -gt 0 ]] && exit 1
exit 0
