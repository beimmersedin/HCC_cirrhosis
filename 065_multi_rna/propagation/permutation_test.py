"""
Network Propagation Permutation Test

동일 크기의 랜덤 seed를 N번 반복하여 null distribution을 만들고,
실제 NP score의 empirical p-value를 계산한다.

Usage:
    python permutation_test.py \
        string_network.txt \
        seed/Hepatocytes_Module1_Cirrhosis_Up_434genes.txt \
        -o results_output/Hepatocytes_Permutation.csv \
        -n 1000 \
        -e 0.1

Input:
    1) network file   : string_network.txt (gene1 \t gene2 \t weight)
    2) seed file      : 유전자 이름이 한 줄에 하나씩

Output:
    CSV with columns:
    - Gene, Observed_Score, Empirical_Pvalue, BH_Adjusted_Pvalue, is_Seed, Rank
"""

import sys
import os
import argparse
import numpy as np
import pandas as pd

from network_propagation import Walker


def run_single_rwr(walker, seed_genes, restart_prob):
    """seed_genes(list) → {gene: score} dict 반환"""
    seed2weight = {g: 1.0 for g in seed_genes if g in walker.dic_node2idx}
    if not seed2weight:
        return {}
    result = walker.run_exp(seed2weight, restart_prob, normalize=True)
    return {node: score for node, score, _ in result}


def main():
    parser = argparse.ArgumentParser(description="Network Propagation Permutation Test")
    parser.add_argument("network", help="Network edge list (gene1 \\t gene2 \\t weight)")
    parser.add_argument("seed", help="Seed gene file (one gene per line)")
    parser.add_argument("-o", required=True, help="Output CSV path")
    parser.add_argument("-n", "--n_perm", type=int, default=1000, help="Number of permutations (default: 1000)")
    parser.add_argument("-e", "--restart_prob", type=float, default=0.1, help="Restart probability (default: 0.1)")
    args = parser.parse_args()

    # ----------------------------------------------------------
    # 1. 네트워크 한 번만 빌드
    # ----------------------------------------------------------
    print(f"[1/4] Building network: {args.network}")
    walker = Walker(args.network, constantWeight=True)
    all_nodes = sorted(walker.OG.nodes())
    node_list = list(walker.OG.nodes())
    print(f"       Network nodes: {len(all_nodes)}")

    # Seed 로드
    with open(args.seed) as f:
        seed_genes = [line.strip() for line in f if line.strip()]
    seeds_in_network = [g for g in seed_genes if g in walker.dic_node2idx]
    seed_size = len(seeds_in_network)
    print(f"       Seed genes: {len(seed_genes)} (in network: {seed_size})")

    # ----------------------------------------------------------
    # 2. 실제 RWR
    # ----------------------------------------------------------
    print(f"[2/4] Running observed RWR (restart_prob={args.restart_prob})...")
    observed = run_single_rwr(walker, seed_genes, args.restart_prob)

    # ----------------------------------------------------------
    # 3. Permutation (순차, Walker 재사용)
    # ----------------------------------------------------------
    print(f"[3/4] Running {args.n_perm} permutations (sequential, single Walker)...")
    count_ge = {gene: 0 for gene in all_nodes}
    rng = np.random.default_rng(42)

    for i in range(args.n_perm):
        random_seeds = rng.choice(node_list, size=seed_size, replace=False).tolist()
        perm_scores = run_single_rwr(walker, random_seeds, args.restart_prob)

        for gene in all_nodes:
            if perm_scores.get(gene, 0.0) >= observed.get(gene, 0.0):
                count_ge[gene] += 1

        if (i + 1) % 50 == 0 or (i + 1) == args.n_perm:
            print(f"       Completed: {i+1}/{args.n_perm}")

    # ----------------------------------------------------------
    # 4. Empirical p-value & BH correction
    # ----------------------------------------------------------
    print(f"[4/4] Computing p-values...")

    records = []
    seed_set = set(seed_genes)
    for gene in all_nodes:
        obs_score = observed.get(gene, 0.0)
        emp_pval = (count_ge[gene] + 1) / (args.n_perm + 1)
        records.append({
            "Gene": gene,
            "Observed_Score": obs_score,
            "Empirical_Pvalue": emp_pval,
            "is_Seed": "Yes" if gene in seed_set else "No",
        })

    df = pd.DataFrame(records)
    df = df.sort_values("Observed_Score", ascending=False).reset_index(drop=True)
    df["Rank"] = df.index + 1

    # Benjamini-Hochberg FDR correction
    n_genes = len(df)
    df = df.sort_values("Empirical_Pvalue").reset_index(drop=True)
    df["BH_rank"] = np.arange(1, n_genes + 1)
    df["BH_Adjusted_Pvalue"] = df["Empirical_Pvalue"] * n_genes / df["BH_rank"]
    df["BH_Adjusted_Pvalue"] = df["BH_Adjusted_Pvalue"][::-1].cummin()[::-1]
    df["BH_Adjusted_Pvalue"] = df["BH_Adjusted_Pvalue"].clip(upper=1.0)
    df = df.drop("BH_rank", axis=1)

    df = df.sort_values("Observed_Score", ascending=False).reset_index(drop=True)
    df["Rank"] = df.index + 1

    os.makedirs(os.path.dirname(args.o) or ".", exist_ok=True)
    df.to_csv(args.o, index=False)

    # Summary
    sig_005 = df[(df["BH_Adjusted_Pvalue"] < 0.05) & (df["is_Seed"] == "No")]
    sig_001 = df[(df["BH_Adjusted_Pvalue"] < 0.01) & (df["is_Seed"] == "No")]
    print(f"\n{'='*60}")
    print(f"Results saved: {args.o}")
    print(f"Total genes: {n_genes}")
    print(f"Significant non-seed (BH adj.p < 0.05): {len(sig_005)}")
    print(f"Significant non-seed (BH adj.p < 0.01): {len(sig_001)}")
    print(f"{'='*60}")

    nonseed = df[df["is_Seed"] == "No"].head(20)
    print("\nTop 20 non-seed genes:")
    print(nonseed[["Rank", "Gene", "Observed_Score", "Empirical_Pvalue", "BH_Adjusted_Pvalue"]].to_string(index=False))


if __name__ == "__main__":
    main()
