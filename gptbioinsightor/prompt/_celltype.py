CELLTYPE_PROMPT = """
<Candidate_celltype>
{candidate}
</Candidate_celltype>
<INSTRUCTIONS>
1. Context Analysis
   Cell type:
      - eg. Common Cell Types in Tissues, like Glial cells are not in blood
   Cell state:
      - experimental condition, dieseas and more
      - eg. Endothelial Cells within tumor tissue: active in neovascularization supporting tumor growth
   Sample realism:
      - consider preparation-related contaminants, residual populations, or technical artifacts (e.g. platelets frequently persist in PBMC isolations)

2. Marker-Based Analysis: High Priority for Cell type prediction

   Cell-type-specific markers: 
      - Broad Category of celltype markers: eg. PTPRC for immune cell
      - Narrow Category of celltype markers: eg. CD79A for B cell; ACTA2 for smooth muscle cell; 
      - marker combinations: eg. PPBP,PF4,CD9,KRT8 for CTC
   Cell-state-specific Markers: Phenotypic Features Acquired with Cell State Changes
      - eg. Fibroblasts/Stellate cell will get ACTA2(α-SMA) when activated

3. Pathway-Based Analysis: High Priority for Cell state prediction
   Cell-state-specific Pathways:
      - eg. Stress response: p53 signaling, HIF-1 signaling
   Cell-type-specific Function: 
      - eg. Immune cells: TCR signaling, BCR signaling

4. Quality Control
   Inclusion Criteria:
      - Strong cell type-specific markers
      - Coherent pathways enrichment
      - Context-appropriate cell type predictions
   
   Exclusion Criteria:
      - Presence of definitive negative markers
      - Conflicting pathway patterns
      - Biologically implausible combinations
      - Fully idealized interpretations that ignore realistic contaminant signatures

   Distinguishing Cell Types and Cell States:
      - Do not confuse cell state features(Marker, pathway) with cell type features.
      - e.g. if (celltypeA with activation state) == (celltypeB), then celltypeA and celltypeB are both candidate cell types

   Potential celltype and state:
      - MUST review each cell type in Candidate_celltype.
      - New celltype could be provided if there are NO matched in Candidate_celltype.

   Validation
      - Provide Gold standard markers(not in cluster {setid} geneset) for cell type validation
</INSTRUCTIONS>
<Reminder>
Focus on existing evidence, Top DEGs are selected for Input Geneset, so Geneset contain limited genes and do not use the lack of classic markers as the basis for your reasoning.
But make reasonable inferential extensions, such as :
- transcription factor regulation
- gene interactions
- metabolic characteristics
- and other plausible speculations.
Prioritize identifying the correct broad lineage for cluster {setid}. Scan <Candidate_celltype> for references to other clusters; only add refined subtypes or states when those clusters share the same lineage, and explicitly state the distinguishing evidence that separates cluster {setid} from its peers.
Always tie your reasoning back to the specified cluster identifier and echo the final label you selected for cluster {setid}.
It should be noted that possible contaminants, residual populations, or technical artifacts may remain due to limitations in isolation, processing, or measurement. Treat these non-ideal components as potential minor confounders when forming conclusions rather than assuming a perfectly pure sample.
Always express every predicted cell type using the standardized Cell Ontology term (only show official CL label, not CL ID).
</Reminder>
<Input>
  <biological_context>
    {background}
  </biological_context>
  <Geneset>
    cluster {setid} markers: {gene} 
  </Geneset>
  <Pathway>
    {pathway}
  </Pathway>
</Input>
<Task>
As BioInsightor, please think about three most likely celltypes for Input based on INSTRUCTIONS. // If you provide a comprehensive and professional analysis, you will publish articles in Cell/Nature/Science
</Task>
"""


CELLTYPE_REPORT = """
Output Format (high socre celltype appear first), without any additional prompt or string. Use standardized Cell Ontology term names for every [CELLTYPE NAME]:
'''

### Celltype Prediction
#### Celltype1: [CELLTYPE NAME] (score: [SCORE])

[Cell state inferring]

**Key Markers**:
- Cell-type-specific: [CELL-SPECIFIC MARKERS]                   

**Evidence and Reasoning** // combine gene marker, enrichment pathway and Cell Origin, >120 words
- [PRIMARY EVIDENCE]
- [SECONDARY EVIDENCE]
- [ADDITIONAL EVIDENCE AS NEEDED]

**Validation**: [Gold Standard MARKERS TO VALIDATE Celltype1]

#### Celltype2: [CELLTYPE NAME] (score: [SCORE])

[Cell state inferring]

**Key Markers**:
- Cell-type-specific: [CELL-SPECIFIC MARKERS]
- Cell-state-specific: [STATE-SPECIFIC MARKERS] 

**Evidence and Reasoning** // combine gene marker, enrichment pathway and Cell Origin, >120 words
- [PRIMARY EVIDENCE]
- [SECONDARY EVIDENCE]
- [ADDITIONAL EVIDENCE AS NEEDED]

**Validation**: [Gold Standard MARKERS TO VALIDATE Celltype2]

#### Celltype3: [CELLTYPE NAME] (score: [SCORE])

[Cell state inferring]

**Key Markers**:
- Cell-type-specific: [CELL-SPECIFIC MARKERS]
- Cell-state-specific: [STATE-SPECIFIC MARKERS] 

**Evidence and Reasoning** // combine gene marker, enrichment pathway and Cell Origin, >120 words
- [PRIMARY EVIDENCE]
- [SECONDARY EVIDENCE]
- [ADDITIONAL EVIDENCE AS NEEDED]

**Validation**: [Gold Standard MARKERS TO VALIDATE Celltype3]

### Novel Insights // combine gene marker, enriched pathway and biological context, >100 words
- [NOTEWORTHY PATTERNS]
- [POTENTIAL NEW FINDINGS]
'''
"""


CELLTYPE_SCORE = """
<Scoring_Criteria_with_pathway>
Marker Profile (60 pts) // common and classical marker will get high score
- Matching cell type or state markers present: max 45  
- Narrow markers of cell type or state present: max 15 
- Share common markers with other cell type or state: -10 // refer to above Candidate_celltype
- Negative markers present: -30

Pathway Profile (20 pts)
- Enriched pathways match cell state: 15
- Enriched pathways match cell type: 5
- Shared Pathway with other cell type or state: -10
- conflicting pathways: -20

Biological Context (20 pts)
- Plausible cell type in Context: 10
- Plausible cell state in Context: 10
- implausible cell type or state in Context: -30
</Scoring_Criteria_with_pathway>

<Scoring_Criteria_without_pathway>
Marker Profile (70 pts)  // commson and classical marker will get high score
- Matching cell type or state markers present: max 50  
- Narrow markers of cell type or state present: max 20 
- Share common markers with other cell type or state: -15 // refer to above Candidate_celltype
- Negative markers present: -30

Biological Context (30 pts)
- Plausible cell type in Context: 15
- Plausible cell state in Context: 15
- implausible cell type or state in Context: -30
</Scoring_Criteria_without_pathway>
<Reminder>
- Score each cell type independently according to Scoring Criteria, without being influenced by the scores of other cell types.
- Do not use the lack of classic markers as the basis for your scoring.
- Negative markers are markers which are impossible to be present in the cell type.
- For matching cell type or state markers, e.g. if scoring B cell, PTPRC is pan-leukocyte marker, match B cell
- For narrow markers of cell type or state, e.g. if scoring B cell, PTPRC is pan-leukocyte marker for immnue cell(broad Category for B cell), not Narrow marker for B cell
- Give full points if there are classical or well-known markers to support cell type or state
- Align your scoring roster with the latest `Candidate Roster` provided earlier in the conversation if `Candidate Roster` exist; treat it as the authoritative list of cell types to evaluate.
</Reminder>
<Task>
Please review and correct above content if anything wrong, then score the each Cell Type Prediction of cluster geneset using a scoring criteria (100 points),
if there is no pathway provided, use Scoring_Criteria_without_pathway, otherwise use Scoring_Criteria_with_pathway, and please notice Reminder
In addition to your thinking process, the final result should be returned with format: response_format,  do not include tag
</Task>
<response_format>
{CELLTYPE1}: {SCORE}
{CELLTYPE2}: {SCORE}
{CELLTYPE3}: {SCORE}
</response_format>
"""


USER_DEFINED_CELLTYPE_SCORE = """
Scoring_Criteria_PROMPT

<Reminder>
- Score each cell type independently according to Scoring Criteria, without being influenced by the scores of other cell types.
- Do not use the lack of classic markers as the basis for your scoring.
- Negative markers are markers which are impossible to be present in the cell type.
- For matching cell type or state markers, e.g. if scoring B cell, PTPRC is pan-leukocyte marker, match B cell
- For narrow markers of cell type or state, e.g. if scoring B cell, PTPRC is pan-leukocyte marker for immnue cell(broad Category for B cell), not Narrow marker for B cell
- Give full points if there are classical or well-known markers to support cell type or state
- Align your scoring roster with the latest `Candidate Roster` provided earlier in the conversation if `Candidate Roster` exist; treat it as the authoritative list of cell types to evaluate.
</Reminder>
<Task>
Please review and correct above content if anything wrong, then score the each Cell Type Prediction of cluster geneset using Scoring_Criteria (100 points),
In addition to your thinking process, the final result should be returned with format: response_format,  do not include tag
</Task>
<response_format>
{CELLTYPE1}: {SCORE}
{CELLTYPE2}: {SCORE}
{CELLTYPE3}: {SCORE}
</response_format>
"""


CELLTYPE_QC_PROMPT = """
## cluster geneset {setid}

### Gene List
Top genes
```
[cluster {setid} gene list]
```
enrichment pathway
```
[cluster {setid} pathway]
```
"""


SUBTYPE_PROMPT = """
<Input>
  <biological_context>
    {background}
  </biological_context>
  <Geneset>
    {genesets}
  </Geneset>
</Input>
<Task>
Hi, BioInsightor! Please determine cell subtypes or state of {celltype} for each geneset of Input Geneset.
Your reasoning process must be based on INSTRUCTION.the final result should be returned with format: response_format,  do not include tag
</Task>

<INSTRUCTION>
1. determine cell subtype according to specific gene markers, please provide evidence and reason
2. give full consideration to Input biological_context, determine the most logical subtype
3. consider the biological_context; speculate the cell state, such as stress, invasive, proliferative, developmental stages, or other transient or dynamically responsive properties
</INSTRUCTION>
<response_format>
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
<response_format>
"""


CHECK_TYPE_PROMPT = """

<Task>
Please determine whether Gene_marker can characterize this Celltype within Biological_context.
Please output in the response_format
</Task>
<INSTRUCTIONS>
  1. review and check if the Celltype is reasonable based on the provided gene in Geneset.
  2. give full consideration to context of cell: pathway and Biological_context if provided
  3. provide new celltype in reason if the Celltype is not reasonable 
</INSTRUCTIONS>
<Input>
  <Celltype>
    {celltype}
  </Celltype>
  <Gene_marker>
    {genest}
  </Gene_marker>
  <Pathway>
    {pathway}
  </Pathway>>
  <Biological_context>
    {background}
  </Biological_context>
</Input>

<response_format>
### [Input celltype]: [ YES or NO ]  // Is this type reasonable?
**reason**: [REASON] // Why do you think it is or isn't reasonable?


### [Corrected celltype] // If the celltype is reasonable, provide original celltype
[Cell state inferring]

**Key Markers**:
- Cell-type-specific: [CELL-SPECIFIC MARKERS]
- Cell-state-specific: [STATE-SPECIFIC MARKERS] 

**Evidence and Reasoning**
- [PRIMARY EVIDENCE]
- [SECONDARY EVIDENCE]
- [ADDITIONAL EVIDENCE AS NEEDED]

**Validation**: [Gold Standard MARKERS TO VALIDATE Celltype1]
</response_format>
"""


PRE_CELLTYPE_PROMPT1 = """
<Input>
  <Biological_context>
    {background}
  </Biological_context>
</Input>
<Task>
Within the <Biological_context>, please list the {num} most likely broad cell types 
(e.g., T cell, B cell, fibroblast, epithelial cell, endothelial cell) that are commonly detected in scRNA-Seq data analysis. 
For each cell type, provide up to two potential cell states relevant to the given biological context.
Only include biologically reasonable cell types for the context, and exclude irrelevant ones.
If relevant, prioritize stromal (fibroblasts, stellate, endothelial, pericytes, vSMCs), 
epithelial/parenchymal, and major immune lineages (T, B, NK, macrophages, dendritic cells, neutrophils, mast cells).
</Task>
"""

PRE_CELLTYPE_PROMPT2 = """
In a scRNA-Seq study derived from the biological context of '{background}', 
list the most likely and common broad cell types expected to be identified. 
For each, also include potential cell states that are biologically reasonable in this context.
"""

PRE_CELLTYPE_MERGE_PROMPT = """
<Task>
Based on the previous three conversation, integrate all above three answers(cell types and states) according to Instruction.
Output result using response_format(without any additional words, string or tag):
</Task>
<Instruction>
1. include previous all cell types and states, if there are any common cell types or states that have NOT been mentioned, please ADD them!
2. correct and exclude unreasonable cell types and states within Biological_context
3. provide classical gene markers for each cell type and state in your output 
</Instruction>
<response_format>
Candidate celltype:
[celltype1]: [classical marker]
   - [cell state1]: [specifc gene markers], [phenotype characteristic]
   - [cell state2]: [specifc gene markers], [phenotype characteristic]
   - [cell state3]: [specifc gene markers], [phenotype characteristic]
   - ...
celltype2: classical marker
   - [cell state1]: [specifc gene markers], [phenotype characteristic]
   - [cell state2]: [specifc gene markers], [phenotype characteristic]
   - ...
......
</response_format>
"""


CELLSTATE_PROMPT = """
<Input>
  <Celltype>
    {celltype}
  </Celltype>
  <biological_context>
    {background}
  </biological_context>    
  <Geneset>
    {gene} 
  </Geneset>
  <Pathway>
    {pathway}
  </Pathway>
</Input>
<Task>
As BioInsightor, please predict the cell state of Input Celltype according to biological_context, Geneset and Pathway  // If you provide a comprehensive and professional analysis, you will publish articles in Cell/Nature/Science
</Task>
"""
