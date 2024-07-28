# GPTBioInsightor

![GPTBioInsightor](https://raw.githubusercontent.com/huang-sh/GPTBioInsightor/main/docs/en/source/img/logo.png)

GPTBioInsightor is a tool designed for single-cell data analysis, particularly beneficial for newcomers to a biological field or those in interdisciplinary areas who may lack sufficient biological background knowledge. GPTBioInsightor harnesses the powerful capabilities of large language models to help people quickly gain knowledge and insight, enhancing their work efficiency.

## Installation

Install GPTBioInsightor from PyPi:
```shell
pip install gptioinsightor
```

## Quick Start


```python
import gptbioinsightor as gbi 

### Set API KEY of LLM 
import os
os.environ['API_KEY'] = "sk-***"

# set background of your data
background = "Cells are PBMCs from a Healthy Donor" 

# make sure you have perform DEG analysis for adata,like: sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
# Here, use qwen2-72b-instruct of Aliyun, but you also can use openai gpt-4o
res = gbi.gptcelltype(adata, background=background, out="celltype.md", topgenes=15,provider="aliyun", model="qwen2-72b-instruct")
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

It will create two markdown files: 
- *.most.md: it describes the most likely cell type for each gene set, annotated using a LLM, including the supporting gene markers, the reason behind the annotation, and potential cellular states.
- *.like.md: it describes the all likely cell type for each gene set, annotated using a LLM, including the supporting gene markers, the reason behind the annotation, and potential cellular states.

like:
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

```

