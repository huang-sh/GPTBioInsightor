from __future__ import annotations

import os, sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from collections.abc import Iterable

from anndata import AnnData

from .core import query_model
from . import utils as ul
from .prompt import *


def _query_celltype(queryid, gene_txt, cluster_num, background, provider, model, base_url, sys_prompt, history=None):
    text = CELLTYPE_PROMPT.format(setid=queryid, gene=gene_txt, setnum=cluster_num,background=background)
    msg = [{"role": "user", "content": text}]
    if history is not None: 
        msg = history + msg
    response = query_model(msg, provider=provider, model=model, base_url=base_url, sys_prompt=sys_prompt)
    return response


def get_celltype(
    input: AnnData | dict, 
    out: Path| str = None, 
    background: str = None, 
    key: str = "rank_genes_groups", 
    topnumber: int = 15, 
    n_jobs: int | None = None,
    provider: str | None = None,
    model: str | None = None,
    group: str | Iterable[str] | None = None,  
    base_url: str | None = None, 
    rm_genes=True
) -> dict:
    """\
    Annotating genesets using LLM, providing cell types, supporting gene markers, reasons, and potential cell state annotations.

    Parameters
    ----------
    input : AnnData | dict
        An AnnData object or geneset dict
    out : Path | str, optional
        output path, by default None
    background : str, optional
        background information of input data, by default None
    key : str, optional
        rank_genes_groups key, by default "rank_genes_groups"
    topnumber : int, optional
        select top gene for analysis, by default 15
    n_jobs : int | None, optional
        set multiple jobs for querying LLM, by default None
    provider : str| None, optional
        LLM provider, by default None
        "openai" for chatgpt
        "aliyun" for qwen
        "moonshot" for kimi
    model : str | None, optional
        set a model based on LLM provider, by default None
    group : str | Iterable, optional
         Which group, by default None
    base_url : str | None, optional
        customized LLM API url, by default None
    rm_genes : bool, optional
        remove rb and mt genes, by default True

    Returns
    -------
    dict
        a celltypes dict
    """
    sys_prompt = SYSTEM_PROMPT
    gene_dic = ul.get_gene_dict(input, group, key, topnumber, rm_genes)
    ot = ul.Outputor(out)
    ot.write("# CellType Analysis")
    chat_msg = ul.get_pre_celltype_chat(len(gene_dic), background, provider, model, base_url, sys_prompt)
    ot.write("## Potential CellType")
    ot.write(f"In scRNA-Seq data background of '{background}', the following Potential CellType to be identified:")
    ot.write(chat_msg[-1]["content"])
    for i in chat_msg:
        print(i["content"])

    if n_jobs is None:
        n_jobs = min(os.cpu_count()//2, len(gene_dic))

    def _aux_func(args):
        gene_txt = "\n".join([f"cluster {k}: {','.join(genes[:topnumber])}" for k,genes in args[1].items()])
        return _query_celltype(args[0], gene_txt, len(args[1]), background, provider, model, base_url, sys_prompt, history=chat_msg)

    celltype_ls = []
    with ThreadPoolExecutor(max_workers=n_jobs) as executor:
        results = executor.map(_aux_func, [(k, gene_dic) for k in gene_dic.keys()])
        for gsid, res in zip(gene_dic.values(), results):
            res = res.strip("```").strip("'''").strip()
            ot.write(res)
            ctn = ul.get_celltype_name(res)
            celltype_ls.append(ctn)
    ot.close()
    if len(gene_dic.keys()) == len(celltype_ls):
        celltype_dic = {k:celltype_ls[idx] for idx, k in enumerate(gene_dic.keys())}
    else: # low-capability model may not get correct output
        celltype_dic = {}
        print("The model may not be producing correct outputs; please try using a better model")
    return celltype_dic


def get_subtype(
    input, 
    out: Path | str | None = None, 
    celltype: str = None,
    background: str = None, 
    group: Iterable[str] | None = None,  
    key: str = "rank_genes_groups", 
    topnumber: int = 15, 
    n_jobs: int | None = None, 
    provider: str | None = None,
    model: str | None = None,
    base_url: str | None = None, 
    rm_genes=True
) -> dict:
    """\
    Annotating cell subtypes using LLM, providing cell types, supporting gene markers, reasons, and potential cell state annotations.
    
    Parameters
    ----------
    input : _type_
        An AnnData object or geneset dict
    out : Path | str | None, optional
        output path, by default None
    celltype : str, optional
        major cell type, by default None
    background : str, optional
        background information of input data, by default None
    group : Iterable[str] | None, optional
        which group, by default None
    key : str, optional
        deg group key, by default "rank_genes_groups"
    topnumber : int, optional
        select top gene for analysis, by default 15
    n_jobs : int | None, optional
        set multiple jobs for querying LLM, by default None
    provider : str | None, optional
        LLM provider, by default None
        "openai" for chatgpt
        "aliyun" for qwen
        "moonshot" for kimi
    model : str | None, optional
        set a model based on LLM provider, by default None
    base_url : str | None, optional
        customized LLM API url, by default None
    rm_genes : bool, optional
        remove rb and mt genes, by default True

    Returns
    -------
    dict
        a cell subtypes dict
    """
    gene_dic = ul.get_gene_dict(input, group, key, topnumber, rm_genes)
    ot = ul.Outputor(out)
    ot.write("# Cell Subtype")

    genesets = [] 
    for k in gene_dic.keys():
        genestr = ",".join(gene_dic[k])
        genesetstr = f"geneset {k}: {genestr}"
        genesets.append(genesetstr)
    genesets_txt = "\n".join(genesets)
    msgs = [
        {"role": "user", "content": SUBTYPE_PROMPT.format(celltype=celltype,genesets=genesets_txt, background=background)}
    ]
    
    response = query_model(msgs, provider=provider, model=model, base_url=base_url, sys_prompt=SYSTEM_PROMPT)
    res_content = response.strip("```").strip("'''")
    ot.write(res_content)
    ot.close()
    subtype_ls = [line.split(":")[1].strip() for line in res_content.split("\n") if line.startswith("###")]

    if len(gene_dic.keys()) == len(subtype_ls):
        subtype_dic = {k:subtype_ls[idx] for idx, k in enumerate(gene_dic.keys())}
    else: # low-capability model may not get correct output
        subtype_dic = {}
        print("The model may not be producing correct outputs; please try using a better model")
    return subtype_dic


def check_celltype(
    input: AnnData | dict, 
    out: Path| str = None, 
    background: str = None, 
    key: str = "rank_genes_groups", 
    topnumber: int = 15, 
    n_jobs: int | None = None, 
    provider: str | None = None,
    model: str | None = None,
    group: str | Iterable[str] | None = None,  
    base_url: str | None = None, 
    rm_genes=True
):
    """\
    Check the reason why genesets are annotated as these celltypes.

    Parameters
    ----------
    input : AnnData | dict
        An AnnData object or geneset dict
    out : Path | str, optional
        output path, by default None
    background : str, optional
        background information of input data, by default None
    key : str, optional
        deg group key, by default "rank_genes_groups"
    topnumber : int, optional
        select top gene for analysis, by default 15
    n_jobs : int | None, optional
        set multiple jobs for querying LLM, by default None
    provider : str | None, optional
        LLM provider, by default None
        "openai" for chatgpt
        "aliyun" for qwen
        "moonshot" for kimi
    model : str | None, optional
        set a model based on LLM provider, by default None
    group : str | Iterable[str] | None, optional
        _description_, by default None
    base_url : str | None, optional
        customized LLM API url, by default None
    rm_genes : bool, optional
        remove rb and mt genes, by default True

    Returns
    -------
    None
    """
    sys_prompt = SYSTEM_PROMPT
    gene_dic = ul.get_gene_dict(input, group, key, topnumber, rm_genes)

    ot = ul.Outputor(out)
    ot.write("# CellType checking")

    genesets = [] 
    for k in gene_dic.keys():
        genestr = ",".join(gene_dic[k])
        genesetstr = f"Celltype: {k}; geneset: {genestr}"
        genesets.append(genesetstr)
    genesets_txt = "\n".join(genesets)
    msgs = [
        {"role": "user", 
        "content": CHECK_TYPE_PROMPT.format(genesets=genesets_txt, background=background)}
    ]
    
    response = query_model(msgs, provider=provider, model=model, base_url=base_url, sys_prompt=sys_prompt)
    res_content = response.strip("```").strip("'''")
    ot.write(res_content)
    ot.close()
