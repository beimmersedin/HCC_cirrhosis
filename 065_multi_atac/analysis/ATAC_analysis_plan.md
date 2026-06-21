# ATAC-seq Analysis Pipeline Plan

## Context
HCC Multiome 프로젝트의 ATAC 데이터(28,859 cells × 76,827 peaks)에 대한 분석 파이프라인 구축.
RNA 분석(QC, UMAP, clustering, annotation)은 완료 상태. ATAC은 QC만 완료.

## 현재 데이터로 각 단계 가능 여부

| 단계 | 가능 여부 | 필요한 것 | 현재 상태 |
|------|-----------|-----------|-----------|
| **1. QC → Dim Reduction → Clustering** | **즉시 가능** | ATAC count matrix + muon/scanpy | h5ad/h5mu 준비 완료 |
| **2. Differential Accessibility (DAR)** | **즉시 가능** | ATAC matrix + cell type annotation + group label | 모두 있음 |
| **3. TF Motif Enrichment** | **추가 작업 필요** | hg38 genome FASTA + motif DB | JASPAR motif는 있으나 genome FASTA 없음 |
| **4. Peak Annotation / Gene Linkage** | **추가 작업 필요** | GTF annotation file + RNA expression | RNA 있으나 GTF 파일 없음 |
| **5. RNA Integration (WNN)** | **즉시 가능** | RNA + ATAC 동일 barcode | h5mu에 둘 다 있음 |

### 즉시 가능한 것 (Step 1, 2, 5)
- ATAC h5ad에 28,859 × 76,827 peak matrix 있음
- Cell type annotation, group (NT/PL/TC), status (Hepatitis/Cirrhosis) 라벨 모두 있음
- RNA와 ATAC이 같은 cell barcode 공유 → 통합 바로 가능
- muon 0.1.7 + scanpy 1.11.5 설치 완료

### 추가 다운로드가 필요한 것 (Step 3, 4)
- **hg38 genome FASTA** (~3GB): motif scanning에 필요. peak 서열을 추출해야 motif 매칭 가능
- **GENCODE GTF** (~50MB): peak을 promoter/gene body/intergenic 등으로 분류하는 데 필요
- 다운로드 명령어 준비 가능 (wget으로 ~10분 소요)

---

## 구현 계획

### Step 1: QC → TF-IDF → LSI → UMAP → Clustering
**도구**: muon + scanpy (설치 완료)

```
1. muon.read_h5mu() → ATAC modality 추출
2. muon.atac.pp.tfidf(atac)  # TF-IDF normalization
3. muon.atac.tl.lsi(atac)    # LSI (PCA on TF-IDF)
4. sc.pp.neighbors(atac, use_rep='X_lsi')  # 1st component 제거 후
5. sc.tl.umap(atac)
6. sc.tl.leiden(atac)
7. RNA cell type annotation을 ATAC에 transfer (같은 barcode)
8. 시각화: UMAP by group/celltype/cluster
```

### Step 2: Differential Accessibility (DAR)
**도구**: scanpy (설치 완료)

```
1. sc.tl.rank_genes_groups(atac, groupby='celltype_assign')  # cell type marker peaks
2. TC vs NT per cell type (subset → Wilcoxon)
3. NT vs PL vs TC 전체 비교
4. Volcano plot, heatmap
```

### Step 3: TF Motif Enrichment
**필요**: hg38 FASTA 다운로드 (~3GB)

```
1. wget hg38.fa.gz (UCSC/Ensembl)
2. Peak 좌표에서 서열 추출 (pysam/pyfaidx)
3. JASPAR motif scanning (muon 내장 or MOODS)
4. DAR peaks에서 enriched motifs 계산
5. Cell type별 TF activity heatmap
```

### Step 4: Peak Annotation / Gene Linkage
**필요**: GENCODE GTF 다운로드 (~50MB)

```
1. wget gencode.v44.annotation.gtf.gz
2. Peak → genomic feature 매핑 (promoter/enhancer/intergenic)
3. Peak-Gene correlation (같은 cell에서 accessibility vs expression)
4. Distance constraint: TSS ±500kb
```

### Step 5: RNA + ATAC Integration (WNN)
**도구**: muon (설치 완료)

```
1. muon.pp.neighbors(mdata)  # multi-modal neighbors
2. muon.tl.umap(mdata)       # joint UMAP
3. RNA-only vs ATAC-only vs WNN UMAP 비교
4. Modality weight 분석
```

## 출력물
- `analysis/03_atac_analysis.ipynb` + `.py`
- `analysis/figures_atac_analysis/`
- 중간 결과 h5mu 저장

## 실행 순서
Step 1 → 2 → 5 (즉시 가능) → 3 → 4 (reference 다운로드 후)
