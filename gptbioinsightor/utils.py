from __future__ import annotations

import os
import sys
from pathlib import Path
from functools import lru_cache

import scanpy as sc
import pandas as pd
from anndata import AnnData

from .exception import ApiKeyMissingError


def get_marker_from_seurat(path: str | Path) -> dict:
    """\
    generate a gene dict from Seurat FindAllMarkers output csv file

    Parameters
    ----------
    path : str | Path
        gene marker csv path

    Returns
    -------
    dict
        gene marker dict
    """
    df = pd.read_csv(path)
    marker_dict = df.groupby('cluster', observed=True)['gene'].agg(list).to_dict()
    return marker_dict


def get_gene_dict(input, group, key, topnumber, rm_genes):
    if isinstance(input, AnnData):
        deg_df = sc.get.rank_genes_groups_df(input, group=group, key=key)
        gene_dic = {}
        for gid, sdf in deg_df.groupby("group", observed=True):
            gene_dic[gid] = sdf["names"].tolist()
    elif isinstance(input, dict):
        gene_dic = input.copy()
    if rm_genes:
        for k in gene_dic.keys():
            gene_dic[k] = [g for g in gene_dic[k] if not g.startswith(('MT-', 'RPL', 'RPS', "ENSG"))]
    for k in gene_dic.keys():
        gene_dic[k] = gene_dic[k][:topnumber]
    return gene_dic


def parse_model(provider, model):
    if provider is None:
        items = model.split(":")
        provider = items[0]
        model = ":".join(items[1:])
    return provider, model


def get_api_key(provider=None):
    if provider is not None:
        API_KEY = os.getenv(f"{provider.upper()}_API_KEY")
        if API_KEY is None:
            API_KEY = os.getenv("API_KEY")
    else:
        API_KEY = os.getenv("API_KEY")
    if API_KEY is None:
        raise ApiKeyMissingError(f"Note: API key not found, please set {provider.upper()}_API_KEY or API_KEY")
    return API_KEY


def set_api_key(api_key: str, provider: str | None = None):
    """\
    set api key for different providers

    Parameters
    ----------
    api_key : str
        api key of the LLM provider 
    provider : str | None, optional
        LLM provider, by default None
    """
    if provider is None:
        os.environ["API_KEY"] = api_key
    else:
        os.environ[f"{provider.upper()}_API_KEY"] = api_key


def get_celltype_name(text):
    for line in text.split("\n"):
        if line.startswith("####"):
            try:
                return line.split(":")[1]
            except IndexError:
                print("LLM doesn't output result accroding predefined format")
         

class Outputor:
    def __init__(self, path: str | Path | None = None) -> None:
        self.path = path
        if self.path is None:
            self.handle = sys.stdout
        else:
            self.handle = open(self.path, "w", encoding="utf-8")
            
    def write(self, text):
        print(text, file=self.handle)
    
    def close(self):
        if self.path is not None: 
            self.handle.close()

# @lru_cache(maxsize=500)
def list_celltypes(num, background, provider, model, base_url, sys_prompt):
    from .core import Agent
    from .prompt import PRE_CELLTYPE_PROMPT1, PRE_CELLTYPE_PROMPT2, PRE_CELLTYPE_MERGE_PROMPT
    from concurrent.futures import ThreadPoolExecutor
    from functools import partial

    query_num = 3
    text = PRE_CELLTYPE_PROMPT1.format(num=num, background=background)
    agent = Agent(model=model, provider=provider, sys_prompt=sys_prompt, base_url=base_url)
    agent.repeat_query(text, n=3, use_context=False)
    chat_msg = [
        {"role": "user", "content": PRE_CELLTYPE_PROMPT2.format(background=background)}, 
        {"role": "assistant", "content": agent.query(PRE_CELLTYPE_MERGE_PROMPT)}
    ]
    return chat_msg


def agent_pipe(agent, pct_txt):
    from .prompt import CELLTYPE_SCORE, CELLTYPE_REPORT
    pre_res = agent.query(pct_txt, use_context=True, add_context=True, use_cache=True)
    scores = agent.query(CELLTYPE_SCORE, use_context=True, add_context=True, use_cache=False)
    report_prompt = CELLTYPE_REPORT.format(score=scores)
    agent.query(report_prompt, use_context=True, add_context=True, use_cache=False)
    return agent.get_history(role="assistant")


def score_heatmap(score_dic, cutoff=0, figsize=(10, 6), cmap='viridis'):
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    df = pd.DataFrame(score_dic).T.apply(pd.to_numeric)
    df = df[df >= cutoff].dropna(axis=1, how='all') 
    plt.figure(figsize=figsize)
    base_size = min(figsize) * 2   
    font_size = max(base_size / max(df.shape), 8) 
    heatmap = sns.heatmap(df, annot=True, cmap=cmap, fmt='g', linewidths=0.5, annot_kws={"size": font_size})
    heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=12)
    plt.title('CellType Score Heatmap')
    plt.xlabel('CellTypes')
    plt.ylabel('Cluster')
    return heatmap


def unify_name(dic, model, provider=None, base_url=None):
    from .core import Agent
    agent = Agent(model=model, provider=provider, sys_prompt=None, base_url=base_url)

    fmt_demo = """
    {
        '0': {'Plasma': '85', 'Memory B Cells': '70', 'Activated B Cells': '60'},
        '1': {'T Cells': '85', 'NK Cells': '25', 'Dendritic Cells': '5'}
    }
    """
    correct_txt = "If there is something format error,like {'1': 'Platelets: 44'}, also correct it as  {'1': {'Platelets': '44'}}"
    text = f"""
    ```JSON
    {str(dic)}
    ```
    Unify the cell type names in this JSON data, using the same format and term to represent the same cell type.
    {correct_txt}
    Only return the corrected JSON format data,without any additional characters or text, such as "", ```or ', like:

    {fmt_demo}

    """
    new_dic_str = agent.query(text, use_context=False, add_context=False, use_cache=True)
    try:
        new_dic = eval(new_dic_str)
    except:
        print("Failed to unify the cell type names")
        new_dic = dic
    return new_dic
