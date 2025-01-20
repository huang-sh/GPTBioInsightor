from __future__ import annotations

import os
from functools import partial
from typing import Union, List, Dict, Iterable
from concurrent.futures import ThreadPoolExecutor
import gseapy as gp
from anndata import AnnData
import scanpy as sc


def enrich(
    input: AnnData | dict, 
    organism: str = 'human',
    gene_sets: Union[List[str], str, Dict[str, str]] = "GO_Biological_Process_2023",
    group: str | Iterable[str] | None = None, 
    key: str = "rank_genes_groups", 
    logfc: float = 0.5, 
    pval: float = 0.2, 
    qval: float = 0.5,
    n_jobs: int = 2
):
    """\
    Enrichment analysis using GSEApy.

    Parameters
    ----------
    input : AnnData | dict
        _description_
    organism : str, optional
        Enrichr supported organism. Select from (human, mouse, yeast, fly, fish, worm), by default 'human'
    gene_sets : Union[List[str], str, Dict[str, str]], optional
        Input Enrichr Libraries, by default "GO_Biological_Process_2023", please check https://maayanlab.cloud/Enrichr/#libraries
    group : str | Iterable[str] | None, optional
        which groups for enrichment, by default None for all
    key : str, optional
        rank_genes_groups key, by default "rank_genes_groups"
    logfc : float, optional
        logfc, by default 0.5
    pval : float, optional
        pval, by default 0.2
    qval : float, optional
        qval, by default 0.5
    n_jobs : int | None, optional
        set multiple jobs for querying, by default 2

    Returns
    -------
    _type_
        _description_
    """
    if isinstance(input, AnnData):
        deg_df = sc.get.rank_genes_groups_df(input, group=group, key=key)
        keep = (deg_df["pvals"]<pval) & (deg_df["pvals_adj"]<qval) & (deg_df["logfoldchanges"]>logfc)
        deg_df = deg_df.loc[keep, ]
        gene_dic = {}
        for gid, sdf in deg_df.groupby("group", observed=True):
            gene_dic[gid] = sdf["names"].tolist()
    elif isinstance(input, dict):
        gene_dic = input.copy()

    if n_jobs is None:
        n_jobs = min(os.cpu_count()//2, len(gene_dic))

    _aux_enrichr = partial(gp.enrichr, gene_sets=gene_sets, organism=organism)
    terms_dic = {}
    with ThreadPoolExecutor(max_workers=n_jobs) as executor:
        results = executor.map(_aux_enrichr, gene_dic.values())
        for k, res in zip(gene_dic.keys(), results):
            edf = res.res2d
            terms_dic[k] = {}
            for gs, sdf in res.res2d.groupby("Gene_set", observed=True):
                if sdf.shape[0]==0:
                    terms_dic[k][gs] = []
                    continue
                sdf = sdf.head(5)
                if gs.startswith("GO_"):
                    terms = sdf.Term.str.split(" \(GO").str[0].tolist()
                elif gs.startswith("WikiPathways"):
                    terms = sdf.Term.str.split("WP").str[0].tolist()
                else:
                    terms = sdf.Term.tolist()
                terms_dic[k][gs] = terms
    return terms_dic
