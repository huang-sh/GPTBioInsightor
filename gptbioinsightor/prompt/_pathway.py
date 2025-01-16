
## source: https://idekerlab.ucsd.edu/gsai/
PATHWAY_NAMING = """
Write a critical analysis of the biological processes performed by this system of interacting proteins.

Base your analysis on prior knowledge available in your training data. After completing your analysis, propose a brief and detailed name for the most prominent biological process performed by the system.
    
After completing your analysis, please also assign a confidence level to the process name you selected. This confidence level should follow the name in parentheses and be one of the following: Low, Medium, or High. A low confidence level indicates the least confidence, while a high confidence level reflects the most confidence. This confidence level helps gauge how accurately the chosen name represents the functions and activities within the system of interacting proteins. When determining your confidence level, consider the proportion of genes in the protein system that participate in the identified biological process. For instance, if you select “Ribosome biogenesis” as the process name but only a few genes in the system contribute to this process, the confidence level should be low compared to a scenario where a majority of the genes are involved in “Ribosome biogenesis.”
     
Put your chosen name at the top of the analysis as 'Name: <name>'. 
Put your confidence in a new line after the name as ‘LLM self-assessed confidence: <confidence level>’

Be concise: Avoid unnecessary words.
Be factual: Do not editorialize.
Be specific: Avoid overly general statements such as ‘the proteins are involved in various cellular processes’.
Group proteins: If group of  proteins has similar functions then discuss their interplay, synergistic, or antagonistic effects, and functional integration within the system, instead of listing facts about individual proteins
Avoid generic process names: Choose detailed and specific names, avoiding generic names like ‘Cellular Signaling and Regulation’.
Contextual relevance: Ensure the analysis is relevant to the specific biological context of the system of interacting proteins if provided.

If you cannot identify a prominent biological process for the proteins in the system, I want you to communicate this in you analysis and name the process: "System of unrelated proteins". Provide a confidence of ‘None’ for a "System of unrelated proteins".
    
To help you in your work, I am providing an example system of interacting proteins and the corresponding example analysis output.

The example system of interacting proteins is:
PDX1, SLC2A2, NKX6-1, GLP1, GCG.

The example analysis output is:
Name: Pancreatic development and glucose homeostasis 
LLM self-assessed confidence: High

1. PDX1 is a homeodomain transcription factor involved in the specification of the early pancreatic epithelium and 
its subsequent differentiation. It activates the transcription of several genes including insulin, somatostatin, glucokinase 
and glucose transporter type 2. It is essential for maintenance of the normal hormone-producing phenotype in the 
pancreatic beta-cell. In pancreatic acinar cells, forms a complex with PBX1b and MEIS2b and mediates the activation 
of the ELA1 enhancer.

2. NKX6-1 is also a transcription factor involved in the development of pancreatic beta-cells during the secondary transition. 
Together with NKX2-2 and IRX3, controls the generation of motor neurons in the neural tube and belongs to the neural progenitor 
factors induced by Sonic Hedgehog (SHH) signals.

3.GCG and GLP1, respectively glucagon and glucagon-like peptide 1, are involved in glucose metabolism and homeostasis. 
GCG raises blood glucose levels by promoting gluconeogenesis and is the counter regulatory hormone of Insulin. 
GLP1 is a potent stimulator of Glucose-Induced Insulin Secretion (GSIS). Plays roles in gastric motility and 
suppresses blood glucagon levels. Promotes growth of the intestinal epithelium and pancreatic islet mass both by islet 
neogenesis and islet cell proliferation.

4. SLC2A2, also known as GLUT2, is a facilitative hexose transporter. In hepatocytes, it mediates bi-directional 
transport of glucose accross the plasma membranes, while in the pancreatic beta-cell, it is the main transporter responsible 
for glucose uptake and part of the cell's glucose-sensing mechanism. It is involved in glucose transport in the small intestine 
and kidney too.

To summarize, the genes in this set are involved in the specification, differentiation, growth and functionality of the pancreas, 
with a particular emphasis on the pancreatic beta-cell. Particularly, the architecture of the pancreatic islet ensures proper 
glucose sensing and homeostasis via a number of different hormones and receptors that can elicit both synergistic and antagonistic 
effects in the pancreas itself and other peripheral tissues.
    

Here are the interacting proteins: 

Proteins:
{geneset}
"""

## source: https://idekerlab.ucsd.edu/gsai/
BIO_PROCESS_PROMPT = """
<task>
Analyze the following set of proteins that have been identified as potentially acting together in a particular biological context. Do their known roles and relationships suggest that there a specific process (or processes) that may be occurring? Do you see any strong relationships to specific disease processes or other phenotypes? Its OK if you don’t see a significant pattern, don’t try to find one if it isn’t there. In that case, just briefly describe the roles and interactions.

Be concise and factual. Only speculate about novel protein functions or interactions when you think there is very strong evidence. Be very clear when you are making a conjecture, a speculation. Use language such as “it might be the case that (speculation) based on (supporting reasoning)”.

Be specific: Avoid overly general statements such as ‘the proteins are involved in various cellular processes’.
</task>

<process>
Identify groups of  proteins with similar functions or which you think interact in this context. For each group, briefly discuss their interplay, synergistic, or antagonistic effects, and functional integration within the system, instead of listing facts about individual proteins.

Then write a brief summary paragraph.

</process>
    
<proteins>
{geneset}
</proteins>

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
Hi, BioInsightor! Please perform the following task based on the provided Term or Pathway.
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
