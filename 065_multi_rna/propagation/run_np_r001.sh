#!/bin/bash

PYTHON_SCRIPT="network_propagation.py"
NETWORK_FILE="string_network.txt"
RESTART_PROB=0.01

SEED_FILES=(
    "seed/Hepatocytes_Module1_Cirrhosis_Up_434genes.txt"
    "seed/T_Cells_Module1_Cirrhosis_Up_313genes.txt"
    "seed/Mesenchymal_Module1_Cirrhosis_Up_391genes.txt"
    "seed/Macrophages_Module1_Cirrhosis_Up_133genes.txt"
    "seed/NK_Cells_Module1_Cirrhosis_Up_116genes.txt"
    "seed/Endothelial_Cells_Module1_Cirrhosis_Up_511genes.txt"
    "seed/Plasma_Cells_Module1_Cirrhosis_Up_18genes.txt"
)

for SEED_FILE in "${SEED_FILES[@]}"; do
    CELLTYPE=$(basename "$SEED_FILE" | sed 's/_Module1_Cirrhosis_Up_.*//')
    OUTPUT_FILE="results_output/${CELLTYPE}_NP_Result_r001.txt"

    echo "=========================================="
    echo "Running: $CELLTYPE (restart_prob=$RESTART_PROB)"
    echo "  Seed: $SEED_FILE"
    echo "  Output: $OUTPUT_FILE"
    echo "=========================================="

    OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    python3 $PYTHON_SCRIPT \
        $NETWORK_FILE \
        "$SEED_FILE" \
        -o "$OUTPUT_FILE" \
        -e $RESTART_PROB \
        -constantWeight True \
        -normalize True

    echo "Top 10:"
    sort -k2 -nr "$OUTPUT_FILE" | head -n 10
    echo ""
done

echo "All done!"
