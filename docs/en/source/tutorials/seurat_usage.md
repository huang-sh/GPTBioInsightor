# Seurat usage

Currently, GPTBioinsightor does not have an R version. However, if you are a Seurat user, you can still use GPTBioinsightor conveniently.

## Usage


### Demo

Here, we will use the classic 10x Genomics PBMC data.

In Unix system, you can download pbmc data like:
```shell
mkdir data
wget http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz -O data/pbmc3k_filtered_gene_bc_matrices.tar.gz
cd data; tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
```

Then perform single-cell data processing in a R environment:
```R

# For more detailed Seurat data processing https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "data/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

pbmc <- RunUMAP(pbmc, dims = 1:10)

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers.fil <- pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.01)

##  save deg marker to a CSV file
write.csv(pbmc.markers.fil, "data/pbmc.markers.fil.csv")
```

Performing cell type annotation using GPTBioinsightor:
```python
# set LLM API KEY
import os
os.environ['API_KEY'] = "sk-***"

import gptbioinsightor as gbi

## read gene marker from Seurat output
result_dict = gbi.get_marker_from_seurat("data/pbmc.markers.fil.csv")

# set background information of data
background = "Cells are PBMCs from a Healthy Donor" 


# Here, use claude-3-5-sonnet-20241022 of anthropic, 
# but you also can use other supported LLM provider.
res = gbi.get_celltype(result_dict, background=background, out="Seurat.claude.celltype.md", topnumber=20, provider="anthropic", model="claude-3-5-sonnet-20241022")
res
# {0: 'T lymphocytes (T cells)',
#  1: 'Monocytes/Macrophages',
#  2: 'T lymphocytes (T cells)',
#  3: 'B Lymphocytes (B Cells)',
#  4: 'Cytotoxic T Cells (CD8+)',
#  5: 'Monocytes/Macrophages',
#  6: 'Natural Killer (NK) Cells',
#  7: 'Dendritic Cells',
#  8: 'Platelets'}
```

You can find more annotation information in `Seurat.qwen.celltype.likely.md` and `Seurat.qwen.celltype.most.md`. The contents of `Seurat.qwen.celltype.most.md` are as follows:
```markdown
# Most Possible celltypes
### Geneset 0: T lymphocytes (T cells)
**gene marker**: CD3D, CD3E, CCR7, LEF1, TCF7
**reason**: The presence of CD3D, CD3E, CCR7, LEF1, and TCF7 strongly supports T lymphocyte identity, with CD3D and CD3E being components of the T-cell receptor complex, CCR7 being involved in T-cell migration, and LEF1 and TCF7 being critical for T-cell development and function.
**cell state/subtype**: Naive or central memory T cells, as indicated by the expression of CCR7, which is typically associated with these subtypes.

### Geneset 1: Monocytes/Macrophages
**gene marker**: S100A8, LGALS2, FCN1, S100A9, CD14, TYROBP, MS4A6A, CST3, LYZ, TYMP, CFD, LST1, LGALS1, AIF1, GSTP1, GRN, GPX1, FTL, S100A6, FTH1
**reason**: Markers such as S100A8, S100A9, CD14, and LYZ are highly indicative of monocyte/macrophage lineage cells and their innate immune response functions.
**cell state/subtype**: Activated or inflammatory state due to the expression of pro-inflammatory markers, suggesting these cells are responding to infection or tissue damage.

### Geneset 2: T lymphocytes (T cells)
**gene marker**: CD3D, IL7R, CD2, TNFRSF4, CD3E, CD27, CD3G, LTB, LAT, MAL
**reason**: The expression of genes like CD3D, CD3E, CD3G, CD2, IL7R, and TNFRSF4, which are common in various subtypes of T cells, supports the identification of this set as T lymphocytes.
**cell state/subtype**: Resting or patrolling state, with the potential inclusion of memory T cells that are poised to respond rapidly to previously encountered pathogens.

### Geneset 3: B Lymphocytes (B Cells)
**gene marker**: CD79A, MS4A1 (CD20), CD79B, TCL1A, VPREB3, HLA-DQA1, HLA-DQB1, CD74, HLA-DRA, FCER2, BANK1, HLA-DRB1, HLA-DPA1, TSPAN13, HLA-DQA2, FCRLA, CD37, HLA-DRB5
**reason**: The presence of genes such as CD79A, CD79B, MS4A1 (CD20), and HLA class II molecules are highly specific to B cell lineage and antigen presentation capabilities.
**cell state/subtype**: Naive or memory B cells in a resting state, ready for activation upon encountering specific antigens.

### Geneset 4: Cytotoxic T Cells (CD8+)
**gene marker**: CD8A, CD8B, GZMA, GZMB, GZMK, GZMM, GZMH, PRF1, NKG7, KLRG1, LAG3, CST7
**reason**: Markers such as CD8A, CD8B, granzymes (GZMA, GZMK, GZMM, GZMH), perforin (PRF1), and NKG7 indicate the presence of cytotoxic T cells, which are capable of killing infected or cancerous cells.
**cell state/subtype**: Activated, effector cytotoxic T cells, likely terminally differentiated effector cells due to the expression of KLRG1.
### Geneset 5 : Monocytes/Macrophages
**gene marker**: CSF1R, CTSL, FCGR3A, HMOX1, SERPINA1
**reason**: The presence of CSF1R, a key receptor for monocyte/macrophage differentiation and survival, along with CTSL, HMOX1, and SERPINA1, which are associated with immune response and oxidative stress, strongly supports the classification of this geneset as characteristic of monocytes/macrophages.
**cell state/subtype**: Activated or inflammatory state due to the expression of genes related to immune response and oxidative stress.

### Geneset 6 : Natural Killer (NK) Cells
**gene marker**: GZMB, PRF1, GNLY, NKG7, XCL1, XCL2, CCL4, KLRD1, FCGR3A
**reason**: The geneset is dominated by markers of cytotoxic function (GZMB, PRF1, GNLY) and chemokines (XCL1, XCL2, CCL4) typical of NK cells, along with NK cell receptors (KLRD1, FCGR3A).
**cell state/subtype**: Activated, cytotoxic state capable of killing virus-infected and tumor cells.

### Geneset 7 : Dendritic Cells
**gene marker**: CLEC4C, CD1C
**reason**: The expression of CLEC4C (Mincle) and CD1C, which are well-established markers for dendritic cells, indicates that this geneset is characteristic of these antigen-presenting cells.
**cell state/subtype**: In a healthy donor, these dendritic cells could be in a resting or patrolling state, ready to present antigens and initiate immune responses.

### Geneset 8 : Platelets
**gene marker**: GP9, ITGA2B, GP1BA, PF4, ITGB3
**reason**: The geneset includes critical components of platelet function such as GP9, ITGA2B, GP1BA, PF4, and ITGB3, which are essential for platelet adhesion, aggregation, and signaling.
**cell state/subtype**: Activated or resting state, as indicated by the presence of genes involved in platelet activation and aggregation, although resting platelets can also express these genes.
```