from __future__ import annotations

import os
from pathlib import Path
from collections.abc import Iterable
from concurrent.futures import ThreadPoolExecutor

from .prompt import *
from . import utils as ul
from .core import query_model, Agent
from .constant import LANG_DIC
from anndata import AnnData

def _query_pathway(queryid, pathways, celltype, background, provider, model, base_url, sys_prompt, lang):
    language = lang if lang not in LANG_DIC else LANG_DIC[lang]
    lang_prompt = LANG_PROMPT.format(language=language)
    prompt_template = PATHWAY_PROMPT + lang_prompt
    pathway_txt = "\n".join(pathways)
    text = prompt_template.format(setid=queryid, pathways=pathway_txt, background=background, celltype=celltype)

    msg = [{"role": "user", "content": text}]
    response = query_model(msg, provider=provider, model=model, base_url=base_url, sys_prompt=sys_prompt)
    return response

def depict_pathway(
    input: dict, 
    out: Path| str = None, 
    celltype_dic: dict = None,
    background: str = None, 
    database: str = None, 
    topnumber: int = 15, 
    n_jobs: int | None = None, 
    provider: str | None = None,
    model: str | None = None,
    base_url=None,
    lang: str = "en"
):
    """
    interpret GO term or pathway with actual biological context using LLM.

    Parameters
    ----------
    input : dict
        pathway input
    out : Path | str, optional
        output path, by default None
    celltype_dic : dict, optional
        celltype of pathway, by default None
    background : str, optional
        background information of pathway input, by default None
    database : str, optional
        pathway database, by default None
    topnumber : int, optional
        select top number  for analysis, by default 15
    n_jobs : int | None, optional
        set multiple jobs for querying LLM, by default None
    provider : str | None, optional
        LLM provider, by default None
        "openai" for chatgpt
        "aliyun" for qwen
        "moonshot" for kimi
    model : str | None, optional
        set a model based on LLM provider, by default None
    base_url : _type_, optional
        customized LLM API url, by default None
    sys_prompt : bool, optional
        use system prompt, by default True
    lang : str, optional
        language setting, by default "en"

    Returns
    -------
    _type_
        None
    """
    ot = ul.Outputor(out)
    ot.write("# Pathway summary")

    if n_jobs is None:
        n_jobs = min(os.cpu_count()//2, len(input))

    def _aux_func(item):
        return _query_pathway(*item, background, provider, model, base_url, SYSTEM_PROMPT, lang)

    with ThreadPoolExecutor(max_workers=n_jobs) as executor:
        iterables = [(k, v[:topnumber], celltype_dic[k]) for k, v in input.items()]
        results = executor.map(_aux_func, iterables)
        for res in results:
            res = res.strip("```").strip("'''")
            ot.write(res)
    ot.close()


def name_pathway(
    input: AnnData | dict, 
    out: Path| str = None, 
    background: str = None, 
    key: str = "rank_genes_groups", 
    topnumber: int = 100, 
    n_jobs: int | None = None,
    provider: str | None = None,
    model: str | None = None,
    group: str | Iterable[str] | None = None,  
    base_url: str | None = None, 
    rm_genes=True
):
    """\
    Process naming and analysis using LLM.

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
    None
    """
    gene_dic = ul.get_gene_dict(input, group, key, topnumber, rm_genes)
    ot = ul.Outputor(out)
    ot.write("# Process Naming and Analysis")
    ot.write("GPTBioInsightor is powered by AI, so mistakes are possible. Review output carefully before use")

    with ThreadPoolExecutor(max_workers=n_jobs) as executor:
        futures = {}
        for k, genes in gene_dic.items():
            agent = Agent(model=model, provider=provider, sys_prompt=None, base_url=base_url)
            query_txt = PATHWAY_NAMING.format(geneset=" ".join(genes))
            future = executor.submit(agent.query, query_txt)
            futures[k] = future
        for k, future in futures.items():
            reps = future.result()
            ot.write(f"## cluster geneset {k}\n")
            ot.write(f"{reps}\n")


def analyse_pathway(
    input: AnnData | dict, 
    out: Path| str = None, 
    background: str = None, 
    key: str = "rank_genes_groups", 
    topnumber: int = 100, 
    n_jobs: int | None = None,
    provider: str | None = None,
    model: str | None = None,
    group: str | Iterable[str] | None = None,  
    base_url: str | None = None, 
    rm_genes=True
):
    """\
    Biological Process Analysis Using LLM.

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
    None
    """
    gene_dic = ul.get_gene_dict(input, group, key, topnumber, rm_genes)
    ot = ul.Outputor(out)
    ot.write("# Biological Process Analysis")
    ot.write("GPTBioInsightor is powered by AI, so mistakes are possible. Review output carefully before use")

    with ThreadPoolExecutor(max_workers=n_jobs) as executor:
        futures = {}
        for k, genes in gene_dic.items():
            agent = Agent(model=model, provider=provider, sys_prompt=None, base_url=base_url)
            query_txt = BIO_PROCESS_PROMPT.format(geneset=" ".join(genes))
            future = executor.submit(agent.query, query_txt)
            futures[k] = future
        for k, future in futures.items():
            reps = future.result()
            ot.write(f"## cluster geneset {k}\n")
            ot.write(f"{reps}\n")
