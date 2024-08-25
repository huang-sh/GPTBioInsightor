from __future__ import annotations

import os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

from .prompt import *
from .core import query_model
from .constant import LANG_DIC


def _query_pathway(queryid, pathways, celltype, background, provider, model, base_url, sys_prompt, lang):
    language = lang if lang not in LANG_DIC else LANG_DIC[lang]
    lang_prompt = LANG_PROMPT.format(language=language)
    prompt_template = PATHWAY_PROMPT + lang_prompt
    pathway_txt = "\n".join(pathways)
    text = prompt_template.format(setid=queryid, pathways=pathway_txt, background=background, celltype=celltype)

    msg = [{"role": "user", "content": text}]
    response = query_model(msg, provider=provider, model=model, base_url=base_url, sys_prompt=sys_prompt)
    return response.choices[0].message.content

def depict_pathway(
    input: dict, 
    out: Path| str = None, 
    celltype_dic: dict = None,
    background: str = None, 
    database: str = None, 
    topnumber: int = 15, 
    n_jobs: int | None = None, 
    provider: str = "openai", 
    model: str | None = None,
    base_url=None,
    sys_prompt=True,
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
    provider : str, optional
        LLM provider, by default "openai"
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
    if out is None:
        out_handle = sys.stdout
    else:
        out_handle = open(out, "w")
        print("# Pathway summary", file=out_handle)
    if n_jobs is None:
        n_jobs = min(os.cpu_count()//2, len(input))

    def _aux_func(item):
        return _query_pathway(*item, background, provider, model, base_url, sys_prompt, lang)

    with ThreadPoolExecutor(max_workers=n_jobs) as executor:
        iterables = [(k, v[:topnumber], celltype_dic[k]) for k, v in input.items()]
        results = executor.map(_aux_func, iterables)
        for res in results:
            res = res.strip("```").strip("'''")
            print(res, file=out_handle)
    if out is not None: 
        out_handle.close()
