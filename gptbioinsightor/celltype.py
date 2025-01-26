from __future__ import annotations

import os, sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from collections.abc import Iterable

from anndata import AnnData

from .core import query_model, Agent
from . import utils as ul
from .prompt import *
from .structure import extract_score


def get_celltype(
    input: AnnData | dict, 
    out: Path| str = None, 
    background: str = None, 
    pathway: dict | None = None,
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
        "deepseek" for DeepSeek
        "anthropic" for claude
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
    genes_ls = []
    for ck, genes in gene_dic.items():
        genes_ls.append(f"   - cluster {ck}: {','.join(genes[:topnumber])}")
    all_gene_txt = "\n".join(genes_ls)
    chat_msg = ul.list_celltype(len(gene_dic), background, provider, model, base_url, sys_prompt)
    ot = ul.Outputor(out)
    ot.write("# CellType Analysis")
    ot.write("GPTBioInsightor is powered by AI, so mistakes are possible. Review output carefully before use")
    ot.write("## Potential CellType")
    ot.write(f"In scRNA-Seq data background of '{background}', the following Potential CellType to be identified:")
    ot.write(chat_msg[-1]["content"])

    with ThreadPoolExecutor(max_workers=n_jobs) as executor:
        futures = {}
        pathway = {} if pathway is None else pathway
        pathway_txt_dic = {}
        for k, genes in gene_dic.items():
            agent = Agent(model=model, provider=provider, sys_prompt=sys_prompt, base_url=base_url)
            gene_txt = f"   - cluster {k}: {','.join(genes[:topnumber])}"
            cluster_pathway = pathway.get(k, {})
            pw_txt = ""
            for db, pw in cluster_pathway.items():
                pw_txt += f"    - {db}: {','.join(pw)}\n"
            pathway_txt_dic[k] = pw_txt
            pct_txt = CELLTYPE_PROMPT.format(
                candidate=chat_msg[-1]["content"], 
                setid=k, gene=gene_txt, 
                setnum=len(gene_dic),
                background=background, 
                pathway=pw_txt
            )
            future = executor.submit(ul.agent_pipe, agent, pct_txt)
            futures[k] = future
        score_dic = {}
        for k, future in futures.items():
            reps = future.result()
            ot.write(f"## cluster geneset {k}\n")
            ot.write(f"### Gene List\n")
            ot.write(f"Top genes\n```\n{','.join(gene_dic[k])}\n```\n")
            ot.write(f"enrichment pathway\n```\n{pathway_txt_dic[k]}```\n")
            ot.write("### celltype thinking\n")
            ot.write(reps[0])
            ot.write("### Score\n")
            ot.write(reps[1])
            ot.write("### Report\n")
            ot.write(reps[2])
            score_ls = extract_score(reps[1], provider, model, base_url).score_ls
            score_dic[k] = {i[0]: i[1] for i in score_ls}
    score_dic = ul.unify_name(score_dic, model, provider, base_url)
    return score_dic


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
    subtype_lines = [line for line in res_content.split("\n") if line.startswith("###")][-len(gene_dic):]
    subtype_ls = [line.split(":")[1].strip() for line in subtype_lines]

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


def get_cellstate(
    input: AnnData | dict, 
    out: Path| str = None, 
    background: str = None, 
    pathway: dict | None = None,
    deg_key: str = "rank_genes_groups", 
    topnumber: int = 15, 
    n_jobs: int | None = None,
    provider: str | None = None,
    model: str | None = None,
    group: str | Iterable[str] | None = None,  
    base_url: str | None = None, 
    rm_genes=True
) -> dict:
    """\
    Annotating cell type state using LLM.

    Parameters
    ----------
    input : AnnData | dict
        An AnnData object or geneset dict
    out : Path | str, optional
        output path, by default None
    background : str, optional
        background information of input data, by default None
    deg_key : str, optional
        deg key, by default "rank_genes_groups"
    topnumber : int, optional
        select top gene for analysis, by default 15
    n_jobs : int | None, optional
        set multiple jobs for querying LLM, by default None
    provider : str| None, optional
        LLM provider, by default None
        "openai" for chatgpt
        "aliyun" for qwen
        "deepseek" for DeepSeek
        "anthropic" for claude
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
    None
    """
    sys_prompt = SYSTEM_PROMPT
    gene_dic = ul.get_gene_dict(input, group, deg_key, topnumber, rm_genes)
    ot = ul.Outputor(out)
    ot.write("# Cell State Analysis")
    ot.write("GPTBioInsightor is powered by AI, so mistakes are possible. Review output carefully before use")

    with ThreadPoolExecutor(max_workers=n_jobs) as executor:
        futures = {}
        pathway = {} if pathway is None else pathway
        pathway_txt_dic = {}
        for k, genes in gene_dic.items():
            agent = Agent(model=model, provider=provider, sys_prompt=sys_prompt, base_url=base_url)
            gene_txt = f"   - cluster {k}: {','.join(genes[:topnumber])}"
            cluster_pathway = pathway.get(k, {})
            pw_txt = ""
            for db, pw in cluster_pathway.items():
                pw_txt += f"    - {db}: {','.join(pw)}\n"
            pathway_txt_dic[k] = pw_txt
            pct_txt = CELLSTATE_PROMPT.format(
                celltype=k,
                gene=gene_txt, 
                background=background, 
                pathway=pw_txt
            )
            future = executor.submit(agent.query, pct_txt)
            futures[k] = future
        score_dic = {}
        for k, future in futures.items():
            reps = future.result()
            ot.write(f"## cluster geneset {k}\n")
            ot.write(f"### Gene List\n")
            ot.write(f"Top genes\n```\n{','.join(gene_dic[k])}\n```\n")
            ot.write(f"enrichment pathway\n```\n{pathway_txt_dic[k]}```\n")
            ot.write("### cell state thinking\n")
            ot.write(reps)
