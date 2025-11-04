# User-defined Cell-type scoring system


The default scoring system is:

| **Evidence Stream** | **Scoring Criterion** | **Rationale** | **Points†** |
|----------------------|-----------------------|----------------|--------------|
| **Marker Profile** | Matching cell-type/state markers detected | Strong evidence the cluster represents the proposed identity | +45 (max) |
|  | Narrow (high-specificity) markers detected | Supports a fine-grained subtype or activation state | +15 (max) |
|  | Shares markers with a different cell type/state | Indicates ambiguity; lowers confidence | −10 |
|  | Negative markers (should not be expressed) present | Contradicts the assignment | −30 |
| **Pathway Profile** | Enriched pathways fit the **cell state** | Captures functional programs (e.g., “interferon response”) | +15 (max) |
|  | Enriched pathways fit the **cell type** | Broad biological match (e.g., “T-cell receptor signaling”) | +5 (max) |
|  | Pathway also enriched in another candidate | Lowers specificity | −10 |
|  | Conflicting pathways | Functional mismatch | −20 |
| **Biological Context** | Cell type plausible in tissue/condition | Aligns with prior knowledge | +10 (max) |
|  | Cell state plausible in tissue/condition | Supports activation/differentiation status | +10 (max) |
|  | Cell type/state biologically implausible | Contradicts experimental context | −30 |


## Customize the scoring prompt

You can override the default weights whenever a project demands a different emphasis (for example, favouring pathway evidence over marker presence). Create a new `score_prompt` string that mirrors the default block structure—each evidence stream is explained in natural language and the model normalises points accordingly.

### Tips for crafting a custom scheme
- Keep the `<Scoring_Criteria>...</Scoring_Criteria>` wrapper; it helps the LLM detect the structure.
- Make the intent explicit (e.g. “rare developmental markers should dominate the score”).
- Specify caps for positive scores and the penalty for conflicting evidence; consistency here avoids ambiguous responses.


### Example: immune-focused scoring

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

The API accepts any well-formed prompt, so you can iteratively tighten your rubric as you review outputs. Focus on how scoring adjustments change the ranking of candidate cell types to ensure the rubric reflects your biological expectations.
