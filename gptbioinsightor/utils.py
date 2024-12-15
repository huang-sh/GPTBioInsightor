from __future__ import annotations

import os
import sys
from pathlib import Path

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



def get_pre_celltype_chat(cluster_num, background, provider, model, base_url, sys_prompt):
    from .core import query_model
    from .prompt import PRE_CELLTYPE_PROMPT
    from concurrent.futures import ThreadPoolExecutor
    text = PRE_CELLTYPE_PROMPT.format(number=cluster_num, background=background)
    msg = [{"role": "user", "content": text}]
    from functools import partial
    _aux_query = partial(query_model, provider=provider, model=model, base_url=base_url, sys_prompt=sys_prompt)
    
    with ThreadPoolExecutor(max_workers=3) as executor:
        results = list(executor.map(_aux_query , [msg, msg, msg]))
    new_text = f"""Please integrate the three list into one list,
    list1:
    '''
    {results[0]}
    '''
    list2:
    '''
    {results[1]}
    '''
    list3:
    '''
    {results[2]}
    '''
    please provide result using follow format:
    ```
    celltype1: classical marker
        - cell state1: specifc gene markers, characteristic
        - cell state2: specifc gene markers, characteristic
    celltype2: classical marker
        - cell state1: specifc gene markers, characteristic
        - cell state2: specifc gene markers, characteristic
    celltype3: classical marker
        - cell state1: specifc gene markers, characteristic
        - cell state2: specifc gene markers, characteristic
    ......
    ```
    """
    msg = [{"role": "user", "content": new_text}]
    response = query_model(msg, provider=provider, model=model, base_url=base_url, sys_prompt=sys_prompt)
    chat_msg = [{"role": "user", "content": text}, {"role": "assistant", "content": response}]
    return chat_msg
