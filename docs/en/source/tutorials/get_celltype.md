# Celltype Annotation
## Installation 

install GPTBioinsightor using pip:

```shell
pip install gptbioinsightor
```

## Usage


### Demo

Here, we will use the classic 10x Genomics PBMC data to demonstrate how to use GPTBioinsightor. GPTBioinsightor is a Python program, and we use Scanpy for single-cell data analysis.

In Unix system, you can download pbmc data like:
```shell
mkdir data
wget http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz -O data/pbmc3k_filtered_gene_bc_matrices.tar.gz
cd data; tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
```

Then perform single-cell data processing in a Python environment:
```python

# For more detailed Scanpy data processing, please refer to  https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering-2017.html

import scanpy as sc

adata = sc.read_10x_mtx(
    "data/filtered_gene_bc_matrices/hg19/",  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file for faster subsequent reading
)

adata.var_names_make_unique()  

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# annotate the group of mitochondrial genes as "mt"
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)

adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()

sc.pp.normalize_total(adata, target_sum=1e4)

sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]

sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

sc.tl.leiden(
    adata,
    resolution=0.9,
    random_state=0,
    flavor="igraph",
    n_iterations=2,
    directed=False,
)
sc.tl.umap(adata)
sc.tl.rank_genes_groups(adata, "leiden", key_added="logreg_deg", method="logreg")
```

Performing cell type annotation using GPTBioinsightor:
```python
# set LLM API KEY
import os
os.environ['API_KEY'] = "sk-***"


import gptbioinsightor as gbi

# set background information of data
background = "Cells are PBMCs from a Healthy Donor" 

# here, I use Aliyun qwen2-72b-instruct
# you can set openai gpt-4o
res = gbi.get_celltype(adata, background=background, out="gbi.qwen.celltype.md", key="logreg_deg", topgenes=15,provider="aliyun", model="qwen2-72b-instruct")
res
# {'0': 'CD4+ T Helper Cells',
#  '1': 'B Cells',
#  '2': 'Monocytes/Macrophages',
#  '3': 'Natural Killer (NK) cells',
#  '4': 'Cytotoxic T Cells (CD8+)',
#  '5': 'Monocytes/Macrophages',
#  '6': 'Dendritic Cells',
#  '7': 'Platelets'}
```

Comparing the results with manual annotations based on classic gene markers
```python
cell_type_name = {
    "0": "CD4 T",
    "1": "B",
    "2": "FCGR3A+ Monocytes",
    "3": "NK",
    "4": "CD8 T",
    "5": "CD14+ Monocytes",
    "6": "Dendritic",
    "7": "Platelet",
}

adata.obs["celltype_manual"] = adata.obs["leiden"].map(
    cell_type_name
)
adata.obs["celltypes_gbi"] = adata.obs["leiden"].map(
    res
)
sc.pl.umap(adata, color=["leiden", "celltype_manual", "celltypes_gbi"], legend_loc="on data", frameon=False)

```
![cell cluster](../img/cell_cluster.png)


You can find more annotation information in `gbi.qwen.celltype.most.md`. The contents of `gbi.qwen.celltype.most.md` are as follows:
```markdown
# Most Possible celltypes
### Geneset 0: CD4+ T Helper Cells
**gene marker**: CD3D, CD3E, CCR7, CD27
**reason**: The presence of CD3D and CD3E, which are integral components of the T-cell receptor complex, along with CCR7 and CD27, which are characteristic of naïve and central memory CD4+ T helper cells, strongly supports this cell type.
**cell state/subtype**: Memory or naïve CD4+ T helper cells in a resting or surveillance state, ready to respond to antigenic challenges.

### Geneset 1: B Cells
**gene marker**: CD79A, MS4A1, CD79B, CD74, CD37
**reason**: These markers are highly specific to B lymphocytes, with CD79A and CD79B being components of the B-cell receptor complex, MS4A1 (CD20) being a well-known B-cell marker, and CD74 and CD37 also being commonly expressed in B cells.
**cell state/subtype**: Mature B cells, potentially activated and capable of antigen presentation, indicated by the presence of HLA-DRA.

### Geneset 2: Monocytes/Macrophages
**gene marker**: FCGR3A, FCER1G, AIF1, LILRA3, MT2A
**reason**: The combination of FCGR3A (CD16), FCER1G (part of Fc receptor complex), AIF1 (involved in macrophage activation), LILRA3 (implicated in immune regulation), and MT2A (a metal detoxification protein) strongly indicates monocytes/macrophages.
**cell state/subtype**: Activated monocytes/macrophages, possibly responding to inflammation or infection.

### Geneset 3: Natural Killer (NK) cells
**gene marker**: GNLY, GZMB, NKG7, PRF1, FCGR3A, TYROBP, XCL2, GZMA
**reason**: This set includes key markers of NK cell function, such as cytotoxic granule proteins (granzymes, perforin), signaling molecules (TYROBP), and the activating receptor CD16 (FCGR3A).
**cell state/subtype**: Activated NK cells, capable of cytotoxic activity against infected or transformed cells.

### Geneset 4: Cytotoxic T Cells (CD8+)
**gene marker**: CCL5, GZMK, NKG7, CST7, CD3D, GZMA, CTSW, CD8A, KLRG1, GZMH, NCR3
**reason**: The presence of CD8A, granzymes (GZMA, GZMK, GZMH), NKG7, and KLRG1 indicates cytotoxic T cells, which are known for their direct killing of infected or cancerous cells.
**cell state/subtype**: Activated or effector CD8+ T cells, potentially engaged in immune surveillance or responding to recent antigen exposure in a healthy individual.
### Geneset 5 : Monocytes/Macrophages
**gene marker**: S100A8, LYZ, S100A9, LGALS2, FCN1, CD14, GSTP1, FTL, TYROBP, GRN, APOBEC3A, GPX1
**reason**: The presence of a comprehensive set of markers, including S100A8, S100A9, CD14, and LYZ, strongly suggests monocytes/macrophages. These markers are indicative of both the cell lineage and the inflammatory state typical of these cells in response to stimuli.
**cell state/subtype**: Activated or inflammatory state due to the presence of alarmins and other inflammatory markers, indicating a response to infection or inflammation.

### Geneset 6 : Dendritic Cells
**gene marker**: HLA-DQA1, HLA-DPB1, HLA-DQB1, HLA-DRA, HLA-DPA1, HLA-DRB1, HLA-DRB5, CD74
**reason**: The high expression of MHC class II genes (HLA-DQA1, HLA-DPB1, etc.) and CD74, which is crucial for MHC class II antigen presentation, is characteristic of dendritic cells. These markers are essential for the function of antigen presentation to T cells.
**cell state/subtype**: Activated or mature dendritic cells, as indicated by the upregulation of MHC class II molecules, which occurs during the maturation process triggered by pathogen recognition.

### Geneset 7 : Platelets
**gene marker**: PPBP, PF4, GP9, GNG11
**reason**: The expression of PPBP, PF4, and GP9 is highly specific to platelets, which are crucial for hemostasis and thrombosis. GNG11, while not exclusive, supports the presence of platelet-related functions.
**cell state/subtype**: Activated or resting platelets. Given the presence of markers associated with platelet function and aggregation, these platelets might be in a state ready to respond to vascular damage or inflammation.

```