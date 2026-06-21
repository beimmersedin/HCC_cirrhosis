import os
import pickle
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist

BASE_DIR = '/data1/project/yeonu/065_multi_rna'
WGCNA_DIR = os.path.join(BASE_DIR, 'WGCNA_80/results_output')
DEG_DIR = os.path.join(BASE_DIR, 'Deg_80/results_output/overlap_genes')
SEED_DIR = os.path.join(BASE_DIR, 'WGCNA_80/results_output/seed')
os.makedirs(SEED_DIR, exist_ok=True)

CELLTYPE_MAP = {
    'Hepatocytes': 'Hepatocyte',
    'T_Cells': 'Tcell',
    'Mesenchymal': 'Mesenchymal',
    'Macrophages': 'Macrophage',
    'NK_Cells': 'NKcell',
    'Endothelial_Cells': 'Endothelial',
    'Plasma_Cells': 'Plasmacell',
    'DCs': 'DC',
    'B_Cells': 'Bcell',
}


def load_deg_genes(deg_name, condition, direction):
    prefix = f'overlap_{deg_name}_{condition}_{direction}_'
    for f in os.listdir(DEG_DIR):
        if f.startswith(prefix) and f.endswith('.txt'):
            with open(os.path.join(DEG_DIR, f)) as fh:
                return set(line.strip() for line in fh if line.strip())
    return set()


def main():
    for ct, deg_name in CELLTYPE_MAP.items():
        # WGCNA 로드
        p_path = os.path.join(WGCNA_DIR, ct, f'{ct}_Network.p')
        with open(p_path, 'rb') as f:
            pywgcna = pickle.load(f)
        MEs = pywgcna.datME
        gene_info = pd.read_csv(
            os.path.join(WGCNA_DIR, ct, 'hub_genes', 'all_genes_wgcna_info.csv'), index_col=0
        )

        module_colors = [col.replace('ME', '') for col in MEs.columns]

        # Eigengene clustering -> 2 clusters
        dist = pdist(MEs.T, metric='correlation')
        Z = linkage(dist, method='average')
        labels = fcluster(Z, t=2, criterion='maxclust')

        c1_colors = [module_colors[i] for i in range(len(module_colors)) if labels[i] == 1]
        c2_colors = [module_colors[i] for i in range(len(module_colors)) if labels[i] == 2]

        genes_c1 = set(gene_info[gene_info['moduleColors'].isin(c1_colors)].index)
        genes_c2 = set(gene_info[gene_info['moduleColors'].isin(c2_colors)].index)

        # DEG 로드
        cases = {
            'case1_Hep_Up': load_deg_genes(deg_name, 'Hepatitis_only', 'Up'),
            'case2_Cirr_Up': load_deg_genes(deg_name, 'Cirrhosis_only', 'Up'),
            'case3_Hep_Down': load_deg_genes(deg_name, 'Hepatitis_only', 'Down'),
            'case4_Cirr_Down': load_deg_genes(deg_name, 'Cirrhosis_only', 'Down'),
        }

        print(f'\n=== {ct} ===')

        for module_label, module_genes in [('Module1', genes_c1), ('Module2', genes_c2)]:
            for case_name, deg_genes in cases.items():
                overlap = sorted(module_genes & deg_genes)
                n = len(overlap)
                fname = f'{ct}_{module_label}_{case_name}_{n}genes.txt'
                fpath = os.path.join(SEED_DIR, fname)

                with open(fpath, 'w') as f:
                    f.write('\n'.join(overlap))

                print(f'  {fname}')

    print(f'\n\nDone! Total files: {len(os.listdir(SEED_DIR))}')
    print(f'Saved to: {SEED_DIR}')


if __name__ == '__main__':
    main()
