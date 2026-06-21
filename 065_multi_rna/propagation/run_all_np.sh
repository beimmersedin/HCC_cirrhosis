#!/bin/bash

# ==========================================
# Network Propagation (RWR) — 전체 seed 실행
# 0genes 파일 제외, 65개 seed
# ==========================================

cd /data1/project/yeonu/065_multi_rna/propagation

NETWORK_FILE="string_network.txt"
SEED_DIR="../WGCNA_80/results_output/seed"
OUT_DIR="results_output"
RESTART_PROB=0.1

mkdir -p "$OUT_DIR"

COUNT=0
TOTAL=$(ls "$SEED_DIR"/*.txt | grep -v "_0genes.txt" | wc -l)

for SEED_FILE in "$SEED_DIR"/*.txt; do
    # 0genes 파일 건너뛰기
    if [[ "$SEED_FILE" == *"_0genes.txt" ]]; then
        continue
    fi

    COUNT=$((COUNT + 1))
    BASENAME=$(basename "$SEED_FILE" .txt)
    OUTPUT_FILE="$OUT_DIR/${BASENAME}_NP_Result.txt"

    echo "=========================================="
    echo "[$COUNT/$TOTAL] $BASENAME"
    echo "  Seed: $SEED_FILE"
    echo "  Output: $OUTPUT_FILE"
    echo "=========================================="

    OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    python3 network_propagation.py \
        "$NETWORK_FILE" \
        "$SEED_FILE" \
        -o "$OUTPUT_FILE" \
        -e $RESTART_PROB \
        -constantWeight True \
        -normalize True

    echo "Top 5:"
    sort -k2 -nr "$OUTPUT_FILE" | head -n 5
    echo ""
done

echo "=========================================="
echo "All NP done! Total: $COUNT"
echo "=========================================="
