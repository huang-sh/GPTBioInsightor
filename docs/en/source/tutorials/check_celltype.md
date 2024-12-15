# Celltype Checking

Sometimes, we want to check whether the cell annotation results from other software are correct, or which gene markers support this cell type, and we can use this `check_celltype`.

## Usage

Checking cell celltype annotation using GPTBioinsightor:
```python
## 
gene_dic = gbi.utils.get_gene_dict(adata, None, "logreg_deg", 15, True)

cg_dic = {
    "CD4 T": gene_dic["0"],
    "B": gene_dic["1"],
    "FCGR3A+ Monocytes": gene_dic["2"],
    "CTC": gene_dic["3"],  # right celltype is NK
    "CD8 T": gene_dic["4"],
    "CD14+ Monocytes": gene_dic["5"],
    "Dendritic": gene_dic["6"],
    "CSC": gene_dic["7"], # right celltype is Platelets
}


res_s = gbi.check_celltype(cg_dic, background=background, topnumber=20, provider="anthropic", model="claude-3-5-sonnet-20241022")
res_s

```

Follow is the output, we can see GPTBioInsightor found the wrong celltypes and give right celltypes:
```markdown
### CD4 T: YES
**reason**: The geneset includes markers such as CD3D, CD3E, and CCR7 which are known to be expressed in CD4+ T cells. Additionally, JUNB and IL32 have been associated with T cell activation and function.

### B: YES
**reason**: The geneset includes CD79A, CD79B, and MS4A1 (also known as CD20), which are well-established markers for B cells. HLA-DRA and other HLA class II molecules are also characteristic of antigen-presenting cells including B cells.

### FCGR3A+ Monocytes: YES
**reason**: The presence of FCGR3A, AIF1, and FCER1G are indicative of monocyte populations, particularly those expressing the Fc gamma receptor IIIa (FCGR3A). These receptors are involved in immune complex recognition and are specific to myeloid lineage cells like monocytes.

### CTC: NO
**reason**: CTC stands for Circulating Tumor Cells, but the geneset does not specifically point to tumor markers or properties typically associated with CTCs. Instead, genes such as GNLY, GZMB, and PRF1 suggest a cytotoxic lymphocyte population, possibly NK cells or CD8+ T cells.

### CD8 T: YES
**reason**: The geneset includes CD8A, GZMA, GZMB, and PRF1, all of which are hallmarks of CD8+ cytotoxic T lymphocytes. NKG7 and KLRG1 further support this cell type due to their roles in cytotoxic granule formation and T cell senescence.

### CD14+ Monocytes: YES
**reason**: The presence of CD14, LYZ, and S100A8/A9 confirms the identity of these cells as monocytes. These genes are specifically upregulated in monocytes and play key roles in innate immunity and inflammation.

### Dendritic: YES
**reason**: The geneset includes HLA-DRA, HLA-DPB1, and HLA-DQB1, which are characteristic of dendritic cells due to their role in antigen presentation. Additionally, FCER1A and CLEC10A are also known to be expressed by dendritic cells.

### CSC: NO
**reason**: CSC stands for Cancer Stem Cells, but the geneset does not contain markers that are typically associated with stemness or cancer stem cell properties. Instead, it contains genes such as PPBP and PF4 which are more commonly associated with platelets and megakaryocytes.

```