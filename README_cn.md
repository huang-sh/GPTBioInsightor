# GPTBioInsightor

<table>
    <tr>
        <td><img src="https://raw.githubusercontent.com/huang-sh/GPTBioInsightor/main/docs/en/source/img/logo.png"></td><td>GPTBioInsightor 是一款专为单细胞数据分析设计的工具，尤其适用于那些刚刚涉足某一生物学领域或从事跨学科研究的研究者，在初期阶段可能不具备充足的生物学背景知识。通过利用大型语言模型的强大能力，GPTBioInsightor 能够帮助用户迅速获取知识和深入见解，从而有效提升其工作效率。</td>
    </tr>
</table>

## 安装

使用pip安装：
```shell
pip install gptbioinsightor
```

## 用法


```python
import gptbioinsightor as gbi 

### 设置大语言模型的API KEY
import os
os.environ['API_KEY'] = "sk-***"

# 设置数据的背景信息
background = "Cells are PBMCs from a Healthy Donor" 

# 确保adata已经进行了差异基因分析，例如: sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
# 我们使用Aliyun 的 qwen2-72b-instruct 进行演示, 但你也可以使用 openai gpt-4o
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

它将创建两个 Markdown 文件：
- *.most.md: 描述每个基因集合中最可能的细胞类型，使用大语言模型进行注释，包括支持的gene markers、注释的理由以及潜在的细胞状态
- *.like.md: 描述每个基因集合中所有可能的细胞类型，使用大语言模型进行注释，包括支持的gene markers、注释的理由以及潜在的细胞状态

例如:
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

