# GPTBioInsightor

<table>
    <tr>
        <td><img src="https://raw.githubusercontent.com/huang-sh/GPTBioInsightor/main/docs/en/source/img/logo.png"></td><td>GPTBioInsightor is a tool designed for single-cell data analysis, particularly beneficial for newcomers to a biological field or those in interdisciplinary areas who may lack sufficient biological background knowledge. GPTBioInsightor utilizes the powerful capabilities of large language models to help people quickly gain knowledge and insight, enhancing their work efficiency.</td>
    </tr>
</table>

## Documention

Please checkout the documentations at:
- English: [GPTBioInsightor docs](https://gptbioinsightor.readthedocs.io/en/latest/)
- 中文: [GPTBioInsightor 文档](https://gptbioinsightor.readthedocs.io/zh-cn/latest/)

## Supported LLM provider
 - openai
 - anthropic
 - openrouter
 - groq
 - aliyun
 - zhipuai
 - siliconflow
 - deepseek
 - perplexity

## Get started
### Installation

Install GPTBioInsightor from PyPi:
```shell
pip install gptbioinsightor
```

### Usage


```python
import gptbioinsightor as gbi 

### Set API KEY of LLM 
import os

os.environ['API_KEY'] = "sk-***"
## or API KEY for anthropic
os.environ['ANTHROPIC_API_KEY'] = "sk-***"


# set background of your data
background = "Human Healthy Donor PBMCs" 

# make sure you have perform DEG analysis for adata,like: sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
# Here, use claude-3-5-sonnet-20241022 of anthropic, but you also can use other supported LLM provider.
res = gbi.get_celltype(adata, background=background, 
                       out="gbi.claude.celltype.md", 
                       topnumber=15,provider="anthropic", 
                       n_jobs=4,model="claude-3-5-sonnet-20241022")
# {'0': 'CD4+ T Helper Cells',
#  '1': 'B Cells',
#  '2': 'Monocytes/Macrophages',
#  '3': 'Natural Killer (NK) cells',
#  '4': 'Cytotoxic T Cells (CD8+)',
#  '5': 'Monocytes/Macrophages',
#  '6': 'Dendritic Cells',
#  '7': 'Platelets'}
```

It will output a markdown file, like:
```markdown
# CellType Analysis
## cluster geneset 0

### Gene List
\```
LDHB, LTB, RGCC, IL32, NOSIP, CD3D, CD3E, TMEM123, VIM, TMEM66, FYB, JUNB, CCR7, CD27, MYL12A
\```

### Celltype Prediction
#### Optimal Celltype: T Cells (likely CD4+ T cells)
**Key Markers**:
- Cell-specific: CD3D, CD3E, CCR7, CD27, FYB
- Context-specific: LTB, JUNB, VIM

**Evidence and Reasoning**
- **PRIMARY EVIDENCE**: The presence of CD3D and CD3E, which are essential components of the T cell receptor complex, strongly indicates a T cell population. These markers are highly specific to T cells.
- **SECONDARY EVIDENCE**: CCR7 is a chemokine receptor that is typically expressed on naïve and central memory T cells, suggesting that these cells may be in a non-activated or memory state.
- **ADDITIONAL EVIDENCE**: CD27 and FYB are also known to be expressed in T cells, particularly in activated T cells. LTB (lymphotoxin beta) and JUNB (a transcription factor) are involved in T cell activation and function.

**Validation**: Other gold standard markers for T cells (not in Geneset 0) include CD4, CD8, and TCRα/β. For CD4+ T cells, additional markers like CD45RA and CD45RO can be used to distinguish between naïve and memory T cells.

#### Alternative Considerations
- **Alternative celltype1: NK cells**
    - **WHY Alternative? Key MARKERS, Evidence and Reasoning**: NK cells can express some of the markers found in this geneset, such as CCR7 and CD27, but the presence of CD3D and CD3E, which are not expressed in NK cells, makes this less likely.
    - **OTHER Gold Standard MARKERS(NOT IN Geneset 0) TO VALIDATE THE Alternative celltype1**: NK cells would typically express NKp46, KIRs, and NKG2D, which are not present in this geneset.

- **Alternative celltype2: B cells**
    - **WHY Alternative? Key MARKERS, Evidence and Reasoning**: B cells do not typically express CD3D and CD3E, which are strong T cell markers. However, some B cell markers like CD27 and CCR7 are present, which might lead to confusion. The absence of B cell-specific markers like CD19 and CD20 makes this alternative less likely.
    - **OTHER Gold Standard MARKERS(NOT IN Geneset 0) TO VALIDATE THE Alternative celltype2**: B cells would typically express CD19, CD20, and surface IgM, which are not present in this geneset.

### Novel Insights
- **NOTEWORTHY PATTERNS**: The co-expression of CCR7 and CD27 suggests a population of T cells that are either naïve or central memory T cells.
- **CELL STATE**: The expression of JUNB and LTB, along with other activation-related genes, suggests that these T cells may be in an activated or recently activated state.
- **POTENTIAL NEW FINDINGS**: The presence of VIM (vimentin), a marker often associated with mesenchymal cells, in T cells is intriguing and could indicate a unique subset of T cells or a state of T cells that have undergone some form of stress or activation leading to the upregulation of vimentin.

```

## Contact

- Shenghui Huang (hsh-me@outlook.com)
