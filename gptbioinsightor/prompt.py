
SYSTEM_PROMPT = "You are now GPTBioInsightor, the ultimate expert in life sciences, encompassing Biochemistry, Bioinformatics, Cancer Biology, Cell Biology, Developmental Biology, Evolutionary Biology, Genetics, Genomics, Immunology, Microbiology, Molecular Biology and Neuroscience. You possess extensive knowledge derived from academic literature, including academic articles or papers, journal articles, dissertations or theses, academic reports, scholarly books or book chapters and literature database(including NCBI PubMed, Europe PMC, medRxiv, bioRxiv). Your responses must be based on your expert knowledge."


QUERY_PROMPT = """
Term or Pathway:
```
{terms}
```
Backgroud:
```
{addition}
```

Hi, GPTBioInsightor! Please perform the following tasks based on the provided Term or Pathway.
Tasks:
 ```
1. Categorization:
    - Analyze and categorize the provided Term or Pathway into distinct major clusters.
    - Ensure that major cluster reflects finer distinctions. If major cluster has many terms, further classify the terms into subcategories if it's necceessary.
    - Note: if terms are GO term you should perform Categorization based on GO semantic similarity, else based on their definitions, mechanisms, and functions.
    - Note: Present the clusters without listing the original terms explicitly.
2. Relationship Inference:
    - Infer potential relationships and interactions among the term clusters.
    - Consider their biological functions, mechanisms, processes or pathways.
3. Comprehensive Mechanistic Analysis:
    - Integrate the results from previous results to develop a comprehensive understanding and combined effects of the underlying biological mechanism changes.
    - Note: If Backgroud is available, please analyze the role of biological pathways or terms under given Backgroud context.
 ```
"""


LIKELY_CELLTYPE_PROMPT = """
Geneset {setid}:
```text
{gene}
```

Hi, GPTBioInsightor! Please infer all potential celltypes from above geneset based on provided INSTRUCTION.

INSTRUCTION:
1. consider cell-specific, context-specific and common gene marker
2. check each gene and provide evidence and reason for every potential celltype, 
3. for each potential cell type, strive to provide a more comprehensive list of gene markers
4. give full consideration to context of cell: BACKGROUND
5. consider the context of celltype: BACKGROUND; speculate the cell state or subtype of celltype, such as stress, invasive, proliferative, developmental stages, or other transient or dynamically responsive properties
6. do not use the lacking of marker as your reason or evidence

BACKGROUND:
{background}

For the output you should follow this format, don't show extra content:
'''
## Geneset {setid}:
[Gene]

### Potential cell types

#### [CELLTYPE1]
**gene marker**: [ALL POSSIBLE GENE MARKER SUPPORTED THIS CELL TYPE]
**reason**: [REASON]
**cell state**: [POTENTIAL CELL STATE UNDER BACKGROUND]

#### [CELLTYPE2]
**gene marker**: [ALL POSSIBLE GENE MARKER SUPPORTED THIS CELL TYPE]
**reason**: [REASON]
**cell state**: [POTENTIAL CELL STATE UNDER BACKGROUND]

#### [CELLTYPE3]
**gene marker**: [ALL POSSIBLE GENE MARKER SUPPORTED THIS CELL TYPE]
**reason**: [REASON]
**cell state**: [POTENTIAL CELL STATE UNDER BACKGROUND]
...

'''
"""


FINAL_CELLTYPE_PROMPT = """
Hi, GPTBioInsightor! Please determine the most likely celltype from multiple potential celltypes of each geneset.Your reasoning process must be based on INSTRUCTION.

INSTRUCTION:
1. you should provide evidence and reason for most likely celltype of each geneset
2. give full consideration to context of cell: BACKGROUND, determine the most logical celltype
4. consider the context of celltype: BACKGROUND; speculate the cell state or subtype of celltype, such as stress, invasive, proliferative, developmental stages, or other transient or dynamically responsive properties
5. consider cell-specific and context-specific gene markers first, then supported gene marker number
6. don't use the lacking of marker as reason or evidence
7. you can also propose a more specific cell type if you are sure it is better than provided celltypes

BACKGROUND:
{background}. Above are {geneset_num} genesets and their potential celltypes, geneseach geneset is highly expressed relative to other gene sets..

For the output you should follow this format:
'''
### [geneset id] : [ CELLTYPE ]  // just celltype name
**gene marker**: [ALL GENE MARKER SUPPORTED THE CELLTYPE]
**reason**: [REASON] // provide detailed evidence and reason for this CELLTYPE
**cell state/subtype**: [POTENTIAL CELL STATE/SUBTYPE UNDER BACKGROUND]

### [geneset id] : [ CELLTYPE ] // just celltype name
**gene marker**: [ALL GENE MARKER SUPPORTED THE CELLTYPE]
**reason**: [REASON] // provide detailed evidence and reason for this CELLTYPE
**cell state/subtype**: [POTENTIAL CELL STATE/SUBTYPE UNDER BACKGROUND]
...

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
### [geneset id] : [ SUBTYPE ] 
** gene marker**: [ALL GENE MARKER SUPPORTED THE CELLTYPE]
** subtype gene marker**: [SPECIFIC GENE MARKER FOR CELL SUBTYPE]
**reason**: [REASON]
**cell state**: [POTENTIAL CELL STATE]

### [geneset id] : [ SUBTYPE ] 
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
