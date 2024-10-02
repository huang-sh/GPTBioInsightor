from __future__ import annotations

import os
from pathlib import Path

import scanpy as sc
import pandas as pd
from anndata import AnnData


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


class ApiKeyMissingError(Exception):
    """Exception raised for missing API_KEY."""
    def __init__(self, message="API_KEY is missing"):
        self.message = message
        super().__init__(self.message)


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
                
            