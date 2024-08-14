# GPTBioInsightor

<table>
    <tr>
        <td><img src="https://raw.githubusercontent.com/huang-sh/GPTBioInsightor/main/docs/en/source/img/logo.png"></td><td>GPTBioInsightor 是一款专为单细胞数据分析设计的工具，尤其适用于那些刚刚涉足某一生物学领域或从事跨学科研究的研究者，在初期阶段可能不具备充足的生物学背景知识。通过利用大型语言模型的强大能力，GPTBioInsightor 能够帮助用户迅速获取知识和深入见解，从而有效提升其工作效率。</td>
    </tr>
</table>

## 文档

详细文档见：[GPTBioInsightor中文文档](https://gptbioinsightor.readthedocs.io/zh-cn/latest/)

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
res = gbi.get_celltype(adata, background=background, out="celltype.md", topgenes=15,provider="aliyun", model="qwen2-72b-instruct")
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

## 联系

- Shenghui Huang (hsh-me@outlook.com)
