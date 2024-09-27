
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


LIKELY_CELLTYPE_PROMPT = """
Geneset {setid}:
```gene list
{gene}
```

Hi, GPTBioInsightor! Please analyze the above geneset and determine three most likely celltypes based on the following INSTRUCTIONS.

INSTRUCTIONS:
0. Evaluate each gene individually, providing evidence and rationale for every potential cell type.
1. Prioritize cell-specific gene markers
2. Consider context-specific gene markers and samples/cells source within context (BACKGROUND)
3. Analyze the celltype context (BACKGROUND) to speculate on cell states, such as stress responses, invasiveness, proliferation rates, developmental stages, or other transient/dynamic properties.
4. Focus on positive evidence; avoid using the absence of markers as primary reasoning.
5. Evaluate the possibility of mixed cell populations or transitional states.
6. If applicable, note any unexpected gene combinations that might suggest novel cell states or types.


BACKGROUND:
{background}

Please format your output as follows, without any additional content:

'''
## Geneset {setid}:

### Gene List
```
[gene list]
```

### Potential Cell Types

#### [CELLTYPE1]
**Gene Markers**:
- cell-specific: [CELL-SPECIFIC GENE MARKERS]
- context-specific: [CONTEXT-SPECIFIC GENE MARKERS]

**Evidence**: [DETAILED EVIDENCE SUPPORTING THIS CELL TYPE]

**Rationale**: [COMPREHENSIVE REASONING]

**Potential Cell State**: [SPECULATED STATE BASED ON BACKGROUND]

#### [CELLTYPE2]
**Gene Markers**:
- cell-specific: [CELL-SPECIFIC GENE MARKERS]
- context-specific: [CONTEXT-SPECIFIC GENE MARKERS]

**Evidence**: [DETAILED EVIDENCE SUPPORTING THIS CELL TYPE]

**Rationale**: [COMPREHENSIVE REASONING]

**Potential Cell State**: [SPECULATED STATE BASED ON BACKGROUND]

#### [CELLTYPE3]
**Gene Markers**:
- cell-specific: [CELL-SPECIFIC GENE MARKERS]
- context-specific: [CONTEXT-SPECIFIC GENE MARKERS]

**Evidence**: [DETAILED EVIDENCE SUPPORTING THIS CELL TYPE]

**Rationale**: [COMPREHENSIVE REASONING]

**Potential Cell State**: [SPECULATED STATE BASED ON BACKGROUND]

### Additional Observations
[ANY NOTEWORTHY PATTERNS, UNUSUAL GENE COMBINATIONS, OR POTENTIAL NOVEL INSIGHTS]
'''
"""


FINAL_CELLTYPE_PROMPT = """
Hi, GPTBioInsightor! Please determine the most likely cell type for each gene set from the provided potential cell types. Your analysis should be based on the following INSTRUCTIONS and context(BACKGROUND) information.

INSTRUCTIONS:
1. Provide comprehensive evidence and reasoning for the most likely cell type of each gene set wtih context(BACKGROUND).
2. Prioritize cell-specific and context-specific gene markers in your analysis.
3. Fully integrate the BACKGROUND information into your analysis to determine the most logical cell type.
4. Speculate on the cell state, considering factors such as stress response, invasiveness, proliferation rate, developmental stage, or other transient/dynamic properties.
5. Do not use the absence of markers as primary evidence; focus on positive evidence.
6. Exclude cell types with clear negative markers present in the gene set.
7. Evaluate the possibility of mixed cell populations or transitional states if the gene set suggests this.
8. Note any unusual gene combinations or expression patterns that might indicate novel cell states or types.

BACKGROUND:
{background}. Above are {geneset_num} genesets and their potential celltypes, geneseach geneset is highly expressed relative to other gene sets..

For the output you should follow this format:
'''
### [geneset id] : [ CELLTYPE NAME]

**Gene Markers**:
- cell-specific: [CELL-SPECIFIC GENE MARKERS]
- context-specific: [CONTEXT-SPECIFIC GENE MARKERS]

**Evidence and Reasoning**:
1. [MAIN EVIDENCE POINT, like cell-specific marker]
2. [SECONDARY EVIDENCE POINT, like context-specific marker ]
3. [ADDITIONAL EVIDENCE POINTS AS NEEDED]

**Cell State/Subtype**: [SPECULATED CELL STATE OR SUBTYPE BASED ON BACKGROUND]
**Alternative Considerations**: [BRIEFLY MENTION OTHER CELL TYPES CONSIDERED AND WHY THEY WERE RULED OUT]


### [geneset id] : [ CELLTYPE NAME ] 

**Gene Markers**:
- cell-specific: [CELL-SPECIFIC GENE MARKERS]
- context-specific: [CONTEXT-SPECIFIC GENE MARKERS]

**Evidence and Reasoning**:
1. [MAIN EVIDENCE POINT, like cell-specific marker]
2. [SECONDARY EVIDENCE POINT, like context-specific marker ]
3. [ADDITIONAL EVIDENCE POINTS AS NEEDED]

**Cell State/Subtype**: [SPECULATED CELL STATE OR SUBTYPE BASED ON BACKGROUND]
**Alternative Considerations**: [BRIEFLY MENTION OTHER CELL TYPES CONSIDERED AND WHY THEY WERE RULED OUT]

[REPEAT FOR EACH GENE SET]

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
