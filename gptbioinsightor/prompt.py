
SYSTEM_PROMPT = "You are now GPTBioInsightor, the ultimate expert in life sciences. You possess extensive knowledge derived from academic articles and literature database(including NCBI PubMed, Europe PMC, medRxiv, bioRxiv). Your responses must be based on your expert knowledge."

LANG_PROMPT = """
Please return {language} text, translated text need to adhere to bioinformatics and biological context of {language} text.
e.g. Pathway translate into Chinese "通路" within bioinformatics and biological context
"""

PATHWAY_PROMPT = """
Pathway:
```
{pathways}
```
Background:
```
Pathway enrichment Background context:
- Sample source: {background}
- celltype: {celltype}(pathways are enriched from high expression genes of {celltype} cells)
```
Hi, GPTBioInsightor! Please perform the following task based on the provided Term or Pathway.
Task:
```
Comprehensive Mechanistic Analysis:
- Analyze and summarize the provided Terms or Pathway, and consider a comprehensive understanding and combined effects of the underlying biological mechanisms.
- If Background is provided,first, consider whether these pathways make biological sense in the given Background context, then analyze the pathways within this context.

```

For the output you should follow this format, don't show extra content:
'''
## Pathway set {setid}:

### Pathway explanation
- [Pathway1]: Explain the pathway, its role, and its function. If applicable, discuss why it is enriched in the Background context and its significance in that context.
- [Pathway2]: Similarly, provide an explanation for each pathway, considering its relevance to the Background context if provided.
...

### Summary
[Comprehensive Analysis]: consider whether these pathways are logical and make biological sense in the given Background context, particularly celltype, first; then summarize the above findings, based on the biological mechanisms and functions represented by the provided pathways and their relevance to the Background context. Propose a coherent biological hypothesis or story.
'''

"""


CELLTYPE_PROMPT = """
Input:
'''
Geneset: 
{gene}

Context: 
{background}. Here are {setnum} genesets for {setnum} different cell clusters up-regulated DEGs. 
'''

Hi, GPTBioInsightor! Please analyze Input and predict the celltypes of geneset cluster {setid} based on the following INSTRUCTIONS.

INSTRUCTIONS:
0. Analyze each gene in cluster {setid}, check cell-specific and context-specific markers
1. prioritize single Gold Standard marker or gene marker combinations for celltype prediction
2. Consider the context of Input for celltype prediction; e.g. tissue, disease, etc.
3. Integrate Context of Input to speculate on some novel insights 
4. Focus on positive evidence, avoid using marker absence as primary reasoning.
5. Exclude celltypes with clear negative markers in the geneset.
6. Consider each cluster has different celltype prediction in most time, exclude celltypes represented by ohter cluster gene markers. 
7. Consider one Optimal celltypes and two alternative celltypes.

Output Format:, without any additional prompt or string:
'''
## cluster geneset {setid}

### Gene List
```
[cluster {setid} gene list]
```

### Celltype Prediction
#### Optimal Celltype: [OPTIMAL CELLTYPE NAME]
**Key Markers**:
- Cell-specific: [CELL-SPECIFIC MARKERS]
- Context-specific: [CONTEXT-SPECIFIC MARKERS]

**Evidence and Reasoning**
- [PRIMARY EVIDENCE]
- [SECONDARY EVIDENCE]
- [ADDITIONAL EVIDENCE AS NEEDED]

**Validation**: [OTHER Gold Standard MARKERS(NOT IN Geneset {setid}) TO VALIDATE THE OPTIMAL CELLTYPE]

#### Alternative Considerations
- Alternative celltype1
    - [WHY Alternative? Key MARKERS, Evidence and Reasoning]
    - [OTHER Gold Standard MARKERS(NOT IN Geneset {setid}) TO VALIDATE THE Alternative celltype1]

- Alternative celltype2
    - [WHY Alternative? Key MARKERS, Evidence and Reasoning]
    - [OTHER Gold Standard MARKERS(NOT IN Geneset {setid}) TO VALIDATE THE Alternative celltype2]

### Novel Insights
- [NOTEWORTHY PATTERNS]
- [CELL STATE]
- [POTENTIAL NEW FINDINGS]
'''
"""


SUBTYPE_PROMPT = """
Hi, GPTBioInsightor! Please determine cell subtypes of {celltype} for each geneset.Your reasoning process must be based on INSTRUCTION.

GENESET:
'''
{genesets}
'''

INSTRUCTION:
1. determine cell subtype according to specific gene markers, please provide evidence and reason
2. give full consideration to context of cell: BACKGROUND, determine the most logical subtype
3. consider the context: BACKGROUND; speculate the cell state, such as stress, invasive, proliferative, developmental stages, or other transient or dynamically responsive properties


BACKGROUND:
{background}

For the output you should follow this format:
'''
### [geneset id] : [ SUBTYPE NAME ] 
** gene marker**: [ALL GENE MARKER SUPPORTED THE CELLTYPE]
** subtype gene marker**: [SPECIFIC GENE MARKER FOR CELL SUBTYPE]
**reason**: [REASON]
**cell state**: [POTENTIAL CELL STATE]

### [geneset id] : [ SUBTYPE NAME ] 
** gene marker**: [ALL GENE MARKER SUPPORTED THE CELLTYPE]
** subtype gene marker**: [SPECIFIC GENE MARKER FOR CELL SUBTYPE]
**reason**: [REASON]
**cell state**: [POTENTIAL CELL STATE]
...

'''
"""



CHECK_TYPE_PROMPT = """
Hi, GPTBioInsightor! Please check celltype of corresponding to each geneset.Your reasoning process must be based on INSTRUCTION.

GENESET:
'''
{genesets}
'''

INSTRUCTION:
1. review and check if the celltype is reasonable based on the provided genes.
2. give full consideration to context of cell: BACKGROUND
3. provide new celltype in reason if the celltype is not reasonable 


BACKGROUND:
{background}

For the output you should follow this format:
'''
### [geneset celltype]: [ YES or NO ]  // Is this type reasonable?
**reason**: [REASON] // Why do you think it is or isn't reasonable?

### [geneset celltype]: [ YES or NO ]  // Is this type reasonable?
**reason**: [REASON] // Why do you think it is or isn't reasonable?
...

'''
"""
