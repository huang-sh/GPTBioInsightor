# 自定义细胞类型评分体系


默认评分体系如下：

| **证据维度** | **评分标准** | **理由** | **分值†** |
|----------------|---------------------------|------------------------|--------------|
| **标记物表达** | 识别到匹配的细胞类型/状态标记物 | 强有力地支持该聚类的身份判断 | +45（上限） |
|                | 识别到高特异性的窄范围标记物   | 指向细分亚型或激活状态         | +15（上限） |
|                | 与其他候选共享标记物           | 说明存在歧义，降低置信度       | −10 |
|                | 检测到不应表达的负向标记物     | 与目标身份矛盾                 | −30 |
| **通路富集**   | 富集通路符合细胞状态           | 捕捉功能程序（例：“干扰素反应”） | +15（上限） |
|                | 富集通路符合细胞类型           | 描述总体生物学匹配（例：“TCR 信号”） | +5（上限） |
|                | 通路同时出现在其他候选中       | 降低特异性                     | −10 |
|                | 富集结果相互矛盾               | 指向功能不匹配                 | −20 |
| **生物学背景** | 在组织/条件下合理的细胞类型     | 与先验知识一致                 | +10（上限） |
|                | 在组织/条件下合理的细胞状态     | 支持激活或分化状态             | +10（上限） |
|                | 生物学上不合理的类型/状态       | 与实验背景矛盾                 | −30 |


## 自定义评分提示

当项目需要强调不同的证据维度时，可以覆盖默认权重。构造一个新的 `score_prompt` 字符串，保持与默认模板相同的结构——逐项说明各证据维度的含义与分值区间，模型会据此标准化评分。

### 编写自定义评分方案时的几点建议
- 保留 `<Scoring_Criteria>...</Scoring_Criteria>` 包裹，方便 LLM 正确解析结构。
- 明确写出设计意图（例如“罕见发育标记权重更高”），减少模型理解偏差。
- 为正向得分与冲突处罚设定清晰上限，保持同一维度内的尺度一致。
- 若需要原始模板，可以调用 `gptbioinsightor.get_score_prompt()` 获取内置版本后再修改。

### 示例：聚焦免疫细胞的评分方案

```python
score_prompt = """
<Scoring_Criteria>
Marker Profile (60 pts)
- Matching immune lineage markers present: max 40
- Activation markers for effector/memory states present: max 20
- Shared markers with non-immune candidates: -15
- Negative markers for the lineage present: -30

Pathway Profile (30 pts)
- Interferon/inflammatory pathways enriched: 15
- Cytotoxic or antigen-presentation pathways enriched: 15
- Pathway overlap with alternative candidates: -10
- Conflicting metabolic pathways: -20

Biological Context (10 pts)
- Plausible immune cell type in the sampled tissue: 5
- Plausible activation state for the condition: 5
- Implausible immune cell in this context: -25
</Scoring_Criteria>
"""

res = gbi.get_celltype(
    adata,
    background=background,
    out="gbi.celltype.md",
    key="deg_key",
    pathway=pathway_dic,
    topnumber=15,
    provider="openai",
    model="gpt-4o",
    n_jobs=4,
    score_prompt=score_prompt,
)
```

API 会接受任意结构清晰的提示文本，建议在复审输出时逐步调整权重，观察候选细胞类型排名的变化，以确保评分准则符合实验预期。
