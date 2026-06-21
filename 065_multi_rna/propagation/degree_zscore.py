"""
Degree-aware Z-score for Network Propagation Results

NP score에서 degree(연결 수) bias를 보정하기 위해,
같은 degree 구간 내에서 z-score를 계산한다.

Usage:
    python degree_zscore.py

Input:
    - string_network.txt (PPI 네트워크)
    - results_output/{CellType}_NP_Final_Results.csv (기존 NP 결과)

Output:
    - results_output/{CellType}_Degree_Zscore.csv

Reference:
    Cowen et al. (2017) Nature Reviews Genetics 18(9):551-562
"""

import os
import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import norm

# ==========================================
# 설정
# ==========================================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
NETWORK_FILE = os.path.join(BASE_DIR, "string_network.txt")
RESULT_DIR = os.path.join(BASE_DIR, "results_output")

CELLTYPES = [
    "Hepatocytes", "T_Cells", "Mesenchymal", "Macrophages",
    "NK_Cells", "Endothelial_Cells", "Plasma_Cells"
]

# degree bin 설정: 5% quantile bins
N_BINS = 20  # 5% 간격 = 20개 bin
MIN_BIN_SIZE = 30  # bin 내 유전자 수가 이보다 적으면 인접 bin과 병합


def build_network(network_file):
    """네트워크 로드 → Graph 객체 반환"""
    G = nx.Graph()
    with open(network_file) as f:
        for line in f:
            g1, g2, w = line.strip().split("\t")
            G.add_edge(g1, g2)
    print(f"Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    return G


def calc_zscore_per_bin(group):
    """degree bin 내에서 z-score 계산"""
    mean = group["NP_Score"].mean()
    std = group["NP_Score"].std()
    if std == 0 or len(group) < 3:
        group["Z_score"] = 0.0
    else:
        group["Z_score"] = (group["NP_Score"] - mean) / std
    return group


def get_seed_neighbors(G, seed_file):
    """seed 유전자의 1차 이웃(direct neighbor) 집합 반환"""
    with open(seed_file) as f:
        seeds = set(line.strip() for line in f if line.strip())
    neighbors = set()
    for s in seeds:
        if s in G:
            neighbors.update(G.neighbors(s))
    neighbors -= seeds  # seed 자체는 제외
    return neighbors


def process_celltype(ct, G, degree_dict, bins):
    """단일 cell type에 대해 degree-aware z-score 계산"""
    np_file = os.path.join(RESULT_DIR, f"{ct}_NP_Final_Results.csv")
    df = pd.read_csv(np_file)

    # seed 1차 이웃 표기
    seed_file = None
    seed_dir = os.path.join(BASE_DIR, "seed")
    for fname in os.listdir(seed_dir):
        if fname.startswith(ct) and fname.endswith(".txt"):
            seed_file = os.path.join(seed_dir, fname)
            break
    if seed_file:
        seed_neighbors = get_seed_neighbors(G, seed_file)
        df["is_Seed_Neighbor"] = df["Gene"].isin(seed_neighbors).map(
            {True: "Yes", False: "No"}
        )
    else:
        df["is_Seed_Neighbor"] = "No"

    # degree 추가
    df["Degree"] = df["Gene"].map(degree_dict).fillna(0).astype(int)

    # degree bin 할당
    df["Degree_bin"] = pd.cut(
        df["Degree"], bins=bins, include_lowest=True, duplicates="drop"
    )

    # 작은 bin 병합
    df = merge_small_bins(df, MIN_BIN_SIZE)

    # bin별 z-score 계산
    df = df.groupby("Degree_bin", group_keys=False, observed=True).apply(
        calc_zscore_per_bin, include_groups=False
    )

    # z-score → p-value (one-tailed)
    df["Z_pvalue"] = 1 - norm.cdf(df["Z_score"])

    # Benjamini-Hochberg FDR correction
    df = df.sort_values("Z_pvalue").reset_index(drop=True)
    n = len(df)
    bh_rank = np.arange(1, n + 1)
    df["Z_adj_pvalue"] = (df["Z_pvalue"] * n / bh_rank).clip(upper=1.0)
    df["Z_adj_pvalue"] = df["Z_adj_pvalue"][::-1].cummin()[::-1]

    # 최종 정렬 (NP score 높은 순)
    df = df.sort_values("NP_Score", ascending=False).reset_index(drop=True)

    # 저장
    out_path = os.path.join(RESULT_DIR, f"{ct}_Degree_Zscore.csv")
    df.to_csv(out_path, index=False)

    # Summary
    nonseed = df[df["is_Seed"] == "No"]
    sig_005 = nonseed[nonseed["Z_adj_pvalue"] < 0.05]
    sig_001 = nonseed[nonseed["Z_adj_pvalue"] < 0.01]

    print(f"\n{ct}:")
    print(f"  Significant non-seed (adj.p < 0.05): {len(sig_005)}")
    print(f"  Significant non-seed (adj.p < 0.01): {len(sig_001)}")
    sig_novel = nonseed[
        (nonseed["Z_adj_pvalue"] < 0.05) & (nonseed["is_Seed_Neighbor"] == "No")
    ]
    print(f"  Novel genes (adj.p<0.05 & not seed neighbor): {len(sig_novel)}")
    print(f"  Top 10 non-seed:")
    top10 = nonseed.head(10)
    for _, row in top10.iterrows():
        nb = "*" if row["is_Seed_Neighbor"] == "Yes" else " "
        print(
            f"   {nb}{row['Gene']:>10}  score={row['NP_Score']:.6f}"
            f"  deg={int(row['Degree']):>4}  z={row['Z_score']:.2f}"
            f"  adj.p={row['Z_adj_pvalue']:.4f}"
        )

    return out_path


def merge_small_bins(df, min_size):
    """bin 내 유전자가 min_size 미만이면 인접 bin과 병합"""
    bin_col = "Degree_bin"
    bins_ordered = sorted(df[bin_col].dropna().unique())

    # 각 bin의 크기 확인 → 작은 bin을 다음 bin과 병합
    merge_map = {}
    i = 0
    while i < len(bins_ordered):
        current = bins_ordered[i]
        count = (df[bin_col] == current).sum()

        # 작은 bin이면 다음 bin과 합침
        if count < min_size and i + 1 < len(bins_ordered):
            next_bin = bins_ordered[i + 1]
            # current와 next_bin을 합친 새 구간
            merged = pd.Interval(current.left, next_bin.right, closed="right")
            merge_map[current] = merged
            merge_map[next_bin] = merged
            # 다음 bin도 소비했으므로 i+2
            bins_ordered[i + 1] = merged
            i += 2
        else:
            i += 1

    if merge_map:
        df[bin_col] = df[bin_col].map(lambda x: merge_map.get(x, x))

    return df


def main():
    # 1. 네트워크 로드 & degree 계산
    G = build_network(NETWORK_FILE)
    degree_dict = dict(G.degree())

    # 2. Degree bin 경계 설정 (전체 유전자의 degree 분포 기반)
    all_degrees = np.array(list(degree_dict.values()))
    bins = np.percentile(all_degrees, np.linspace(0, 100, N_BINS + 1))
    bins = np.unique(bins)
    print(f"Degree bins: {len(bins)-1} bins, range {int(bins[0])}~{int(bins[-1])}")

    # 3. Cell type별 처리
    for ct in CELLTYPES:
        process_celltype(ct, G, degree_dict, bins)

    print(f"\n{'='*60}")
    print("Done! All degree-aware z-score results saved.")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
