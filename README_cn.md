# GPTBioInsightor

<table>
    <tr>
        <td><img src="https://raw.githubusercontent.com/huang-sh/GPTBioInsightor/main/docs/en/source/img/logo.png"></td><td>GPTBioInsightor 是一款专为单细胞数据分析设计的工具，尤其适用于那些刚刚涉足某一生物学领域或从事跨学科研究的研究者，在初期阶段可能不具备充足的生物学背景知识。通过利用大型语言模型的强大能力，GPTBioInsightor 能够帮助用户迅速获取知识和深入见解，从而有效提升其工作效率。</td>
    </tr>
</table>

## 文档

详细文档见：[GPTBioInsightor中文文档](https://gptbioinsightor.readthedocs.io/zh-cn/latest/)

## 支持 LLM 提供者
 - openai
 - anthropic
 - openrouter
 - groq
 - aliyun
 - zhipuai
 - siliconflow
 - deepseek
 - perplexity


## 快速开始

### 安装

使用pip安装：
```shell
pip install gptbioinsightor
```

### 用法


```python
import gptbioinsightor as gbi 

### 设置大语言模型的API KEY
import os
os.environ['API_KEY'] = "sk-***"

# 设置数据的背景信息
background = "Cells are PBMCs from a Healthy Donor" 

# 确保adata已经进行了差异基因分析，例如: sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
# 我们使用Aliyun 的 qwen2-72b-instruct 进行演示, 但你也可以使用 openai gpt-4o
res = gbi.get_celltype(adata, background=background, out="celltype.md", topnumber=15,provider="aliyun", model="qwen2-72b-instruct")
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

它将创建两个 Markdown 文件：
- *.most.md: 描述每个基因集合中最可能的细胞类型，使用大语言模型进行注释，包括支持的gene markers、注释的理由以及潜在的细胞状态
- *.like.md: 描述每个基因集合中所有可能的细胞类型，使用大语言模型进行注释，包括支持的gene markers、注释的理由以及潜在的细胞状态

例如:
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

## 联系

- Shenghui Huang (hsh-me@outlook.com)
