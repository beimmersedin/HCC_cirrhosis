# Differential Accessibility (DA) Analysis Summary

## Analysis: TC (Tumor Core) — Cirrhosis vs Hepatitis
- **Method**: Signac `FindMarkers` (Logistic Regression, latent.vars = nCount_ATAC, min.pct = 0.05)
- **Significance**: adjusted p-value < 0.05

---

## Cell Count per Cell Type

| Cell Type | Total | Cirrhosis | Hepatitis | Cirr % | Status |
|---|---:|---:|---:|---:|---|
| T_Cells | 3,126 | 1,213 | 1,913 | 38.8% | Analyzed |
| Hepatocytes | 2,734 | 175 | 2,559 | 6.4% | Analyzed |
| Macrophages | 2,170 | 470 | 1,700 | 21.7% | Analyzed |
| Mesenchymal | 1,778 | 512 | 1,266 | 28.8% | Analyzed |
| NK_Cells | 1,722 | 211 | 1,511 | 12.3% | Analyzed |
| DCs | 347 | 56 | 291 | 16.1% | Analyzed |
| Endothelial_Cells | 317 | 1 | 316 | 0.3% | **SKIP** (Cirrhosis < 10) |
| Plasma_Cells | 234 | 100 | 134 | 42.7% | Analyzed |
| B_Cells | 73 | 22 | 51 | 30.1% | Analyzed |

---

## DA Results Summary

| Cell Type | Peaks Tested | Significant DARs | Up in Cirrhosis (FC>0) | Up in Hepatitis (FC<0) | Elapsed |
|---|---:|---:|---:|---:|---|
| **Hepatocytes** | 17,938 | **689** | 687 | 2 | 8m 28s |
| **T_Cells** | 9,049 | **83** | 77 | 6 | 3m 33s |
| **Mesenchymal** | 12,969 | **38** | 37 | 1 | 3m 22s |
| **Macrophages** | 12,165 | **29** | 27 | 2 | 4m 18s |
| **NK_Cells** | 11,835 | **15** | 14 | 1 | 3m 21s |
| **DCs** | 17,195 | **1** | 1 | 0 | 1m 34s |
| **Plasma_Cells** | 13,847 | **0** | 0 | 0 | 1m 09s |
| **B_Cells** | 15,307 | **0** | 0 | 0 | 55s |
| Endothelial_Cells | — | — | — | — | SKIP |

---

## Key Observations

1. **Hepatocytes**가 689개로 압도적으로 많은 significant DAR을 보임. 간의 주요 실질세포로서 Cirrhosis 진행에 따른 chromatin remodeling이 가장 뚜렷함.
2. **방향성**: 전체 significant DAR 855개 중 843개(98.6%)가 Cirrhosis에서 accessibility 증가 (avg_log2FC > 0). Hepatitis → Cirrhosis 진행 시 전반적 chromatin opening이 관찰됨.
3. **chr8-8227872-8228698** peak이 T_Cells, Hepatocytes, Macrophages, NK_Cells 등 여러 세포유형에서 공통 top hit으로 반복 등장. 세포유형 공통 regulatory locus로 추정.
4. **B_Cells, Plasma_Cells**: significant DAR 0개. 세포 수 부족(B: 73개) 또는 chromatin landscape가 두 조건 간 안정적일 가능성.
5. **Endothelial_Cells**: Cirrhosis=1개로 비교 불가하여 skip됨.
