
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
1. you should check each gene and provide evidence and reason for every potential cell type, 
2. for each potential cell type, strive to provide a more comprehensive list of gene markers
3. give full consideration to BACKGROUND information
4. give full consideration the cell states that the cell type is in, such as the stress-like state, invasive state, proliferation or other transient or dynamically responsive property of a cell to a background context—rather than a cell type
5. consider context-specific gene marker and common gene marker
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
Hi, GPTBioInsightor! Please determine the most likely celltype from multiple potential celltypes of each geneset.Your reasoning process must be based on following INSTRUCTION.

INSTRUCTION:
1. you should provide evidence and reason for most likely celltype of each geneset!
2. give full consideration to BACKGROUND, determine the most logical celltype! 
3. give full consideration to cell states that the cell type are in, such as the stress-like state, invasive state or other transient or dynamically responsive property of a cell to a background context—rather than a cell type!
5. check if there are more genes within the respective gene set that support most likely celltype, besides the aforementioned gene marker!
6. do not use the lacking of marker as your reason or evidence!
 

BACKGROUND:
{background}. Above are {geneset_num} genesets and their potential celltypes, genes of each geneset are upregulated in a corresponding cell type compared to other cell types.

For the output you should follow this format:
'''
### [geneset id] : [ MOST LIKELY CELLTYPE ]
**gene marker**: [ALL GENE MARKER SUPPORTED THE CELLTYPE]
**reason**: [REASON]
**cell state**: [POTENTIAL CELL STATE UNDER BACKGROUND]

### [geneset id] : [ MOST LIKELY CELLTYPE ]
**gene marker**: [ALL GENE MARKER SUPPORTED THE CELLTYPE]
**reason**: [REASON]
**cell state**: [POTENTIAL CELL STATE UNDER BACKGROUND]

### [geneset id] : [ MOST LIKELY CELLTYPE ]
**gene marker**: [ALL GENE MARKER SUPPORTED THE CELLTYPE]
**reason**: [REASON]
**cell state**: [POTENTIAL CELL STATE UNDER BACKGROUND]

...

'''
"""
