# HCC Cirrhosis Multiome Analysis

간세포암(HCC) 환자 유래 조직의 **single-cell Multiome (scRNA-seq + scATAC-seq)** 데이터를 이용해
**Cirrhosis vs Hepatitis** 배경에서의 세포타입별 발현·염색질 접근성 변화와 그 조절 메커니즘을 규명하는 프로젝트입니다.

- **샘플 축**: Cirrhosis vs Hepatitis
- **조직 구획**: NT (Non-Tumor) / PL (Premalignant) / TC (Tumor Core)
- **규모**: ~28,859 cells × 76,827 peaks (ATAC), 32,913 genes (RNA)
- **세포타입(9종)**: Hepatocyte, T_Cells, NK_Cells, Macrophages, DCs, B_Cells, Plasma_Cells, Mesenchymal, Endothelial_Cells

---

## 디렉토리 구조

```
hcc_cirrhosis/
├── 065_multi_rna/                  # scRNA-seq 분석
│   ├── data/                       # count matrix(.mtx), metadata, combined_final_annotated.h5ad
│   ├── Deg/                        # MAST DEG (세포타입 × Cirrhosis/Hepatitis, NT vs TC)
│   ├── Deg_80/                     # DEG (80% zero-rate 필터 버전) + boxplot
│   ├── WGCNA_80/                   # PyWGCNA: pseudobulk → 모듈 → hub/seed 유전자
│   └── propagation/                # STRING 네트워크 기반 network propagation (RWR) + permutation
│
├── 065_multi_atac/                 # scATAC-seq 분석
│   ├── analysis/                   # QC, overview, chromVAR downstream (.ipynb/.py) + figures
│   ├── DA/                         # Differential Accessibility (DAR): TC Cirrhosis vs Hepatitis
│   ├── TFmotif/                    # JASPAR motif enrichment + chromVAR TF activity
│   ├── Linkage/                    # Peak–Gene linkage (Signac LinkPeaks)
│   └── GRN/                        # TF–Gene 상관 기반 Gene Regulatory Network
│
├── 260420_최연우_cirrhosis.pptx    # 발표 자료
├── 260420_최연우_cirrhosis.pdf
├── environment.yml                 # conda 환경 (Python)
├── environment_crossplatform.yml
└── requirements.txt
```

각 분석 폴더의 `results_output/`에 산출 CSV/PNG/notebook이 저장됩니다.
R 환경은 `renv.lock`(065_multi_rna, 065_multi_atac)으로 관리합니다.

---

## 분석 파이프라인

### RNA (`065_multi_rna/`)

| 단계 | 도구 | 설명 |
|------|------|------|
| **DEG** | MAST (R) | 세포타입별 NT vs TC 차등발현, Cirrhosis/Hepatitis 각각. `volcano_plot.py`로 시각화 |
| **WGCNA** | PyWGCNA | pseudobulk → co-expression 모듈 → hub gene → seed gene 추출 |
| **Network Propagation** | RWR (STRING v12.0) | 모듈/DEG seed에서 random walk with restart, permutation test로 유의성 검정 |

### ATAC (`065_multi_atac/`)

| 단계 | 도구 | 설명 |
|------|------|------|
| **QC / Overview** | muon + scanpy | TF-IDF → LSI → UMAP → clustering, QC figures |
| **DA (DAR)** | Signac `FindMarkers` (LR) | TC에서 Cirrhosis vs Hepatitis 세포타입별 차등 접근성 peak |
| **TF Motif** | JASPAR + chromVAR | DAR의 motif enrichment, 세포타입별 TF activity |
| **Linkage** | Signac `LinkPeaks` | Peak–Gene correlation (TSS 거리 제약) |
| **GRN** | Spearman correlation | TF–target gene edge, master regulator 추출 |

---

## 주요 결과

- **DAR**: Hepatocyte가 **689개**로 압도적 (전체 855개 중). 전체 유의 DAR의 **98.6%가 Cirrhosis에서 accessibility 증가** → Hepatitis→Cirrhosis 진행 시 전반적 chromatin opening.
- **GRN**: Hepatocyte master regulator로 **HNF4A** (target 69개), 이어 CTCF, RXR 계열(NR1H2::RXRA, RXRG, RXRB).
- **공통 locus**: `chr8-8227872-8228698` peak이 여러 세포타입에서 공통 top hit.

세포타입별 DAR / cell count 상세는 [`065_multi_atac/DA/results_output/DA_summary.md`](065_multi_atac/DA/results_output/DA_summary.md) 참고.

---

## ⚠️ 알려진 한계 (해석 시 주의)

1. **Sample-level confounding (가장 중요)** — Cirrhosis가 사실상 환자 1명(C1_TC)에 가까움.
   현재 DAR/HNF4A 결과가 "Cirrhosis effect"가 아닌 "특정 환자 vs 나머지" 차이일 수 있음.
   → GeoMx 등 독립 데이터 검증 또는 limitation 명시 필요.
2. **Motif enrichment 필터 미적용** — 일부 TF 스크립트에서 enrichment p-value 필터가 주석 처리된 상태로 downstream 진행됨.
3. **"Causal" 용어** — 파이프라인 전반이 correlative evidence이므로 "causal" 네이밍은 재검토 권장.
4. **코드 중복** — 세포타입별 스크립트가 copy-paste 구조 → parameterized script로 모듈화 권장.
5. **통계적 rigor** — cell-level LR은 sample 내 비독립성 무시. pseudobulk(DESeq2/edgeR) 고려.

> 상세 평가는 [`065_multi_atac/참고.txt`](065_multi_atac/참고.txt)에 기록되어 있습니다.

---

## 환경 설정

```bash
# conda (Python: scanpy, muon, PyWGCNA 등)
conda env create -f environment.yml
conda activate hcc

# R (renv)
cd 065_multi_rna   # 또는 065_multi_atac
R -e 'renv::restore()'
```

주요 데이터(`.h5ad`, `.mtx`, STRING `.txt` 등)는 `.gitignore`로 제외되어 있을 수 있으니 별도 확보 필요.
