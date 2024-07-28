from __future__ import annotations

import os, sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from collections.abc import Iterable

import scanpy as sc
from anndata import AnnData

from .core import query_model
from .prompt import LIKELY_CELLTYPE_PROMPT, FINAL_CELLTYPE_PROMPT


def _query_celltype(genes, queryid, background, provider, model, base_url, sys_prompt):
    text = LIKELY_CELLTYPE_PROMPT.format(setid=queryid, gene=",".join(genes), background=background)
    msg = [{"role": "user", "content": text}]
    response = query_model(msg, provider=provider, model=model, base_url=base_url, sys_prompt=sys_prompt)
    return response.choices[0].message.content


def gptcelltype(
    input: AnnData | dict, 
    out: Path| str = None, 
    background: str = None, 
    key: str = "rank_genes_groups", 
    topgenes: int = 15, 
    n_jobs: int | None = None, 
    provider: str = "openai", 
    model: str | None = None,
    group: str | Iterable[str] | None = None,  
    base_url: str | None = None, 
    rm_genes=True, 
    sys_prompt=True
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
    topgenes : int, optional
        select top gene for analysis, by default 15
    n_jobs : int | None, optional
        set multiple jobs for querying LLM, by default None
    provider : str, optional
        LLM provider, by default "openai"
        "openai" for chatgpt
        "aliyun" for qwen
        "moonshot" for kimi
    model : str | None, optional
        set a model based on LLM provider, by default None
    group : str | Iterable, optional
         Which group, by default None
    base_url : str | None, optional
        customized LLM API url by default None
    rm_genes : bool, optional
        rm rb and mt genes, by default True
    sys_prompt : bool, optional
        use system prompt, by default True

    Returns
    -------
    dict
        a celltypes dict
    """
    if isinstance(input, AnnData):
        deg_df = sc.get.rank_genes_groups_df(input, group=group, key=key)
        gene_dic = {}
        for gid, sdf in deg_df.groupby("group"):
            gene_dic[gid] = sdf["names"].tolist()
    elif isinstance(input, dict):
        gene_dic = input
    if rm_genes:
        for k in gene_dic.keys():
            gene_dic[k] = [g for g in gene_dic[k] if not g.startswith(('MT-', 'RPL', 'RPS'))]
    for k in gene_dic.keys():
            gene_dic[k] = gene_dic[k][:topgenes]
    if out is None:
        likely_handle, most_handle = sys.stdout, sys.stdout
    else:
        out = Path(out)
        likely_path = out.with_suffix(f".likely{out.suffix}")
        most_path = out.with_suffix(f".most{out.suffix}")
        likely_handle = open(likely_path, "w")
        most_handle = open(most_path, "w")
        print("# All possible celltypes", file=likely_handle)
        print("\n\n# Most Possible celltypes", file=most_handle)
        
    if n_jobs is None:
        n_jobs = min(os.cpu_count()//2, len(gene_dic))

    def _aux_func(item):
        return _query_celltype(item[1][:topgenes], item[0], background, provider, model, base_url, sys_prompt)
        
    likely_res_ls = []
    with ThreadPoolExecutor(max_workers=n_jobs) as executor:
        results = executor.map(_aux_func, gene_dic.items())
        for (gsid, genes), res in zip(gene_dic.items(), results):
            res = res.strip("```").strip("'''")
            print(res, file=likely_handle)
            likely_res_ls.append(res)

    celltype_ls = []
    for i in range(0, len(likely_res_ls), 5):
        part_likely_res = "".join(likely_res_ls[i:i + 5])
        msgs = [
            {"role": "user", "content": part_likely_res},
            {"role": "user", "content": FINAL_CELLTYPE_PROMPT.format(geneset_num=len(likely_res_ls[i:i + 5]), background=background)}
        ]

        response = query_model(msgs, provider=provider, model=model, base_url=base_url, sys_prompt=sys_prompt)
        
        res_content = response.choices[0].message.content.strip("```").strip("'''")
        print(res_content, file=most_handle)
        celltype_ls.extend([line.split(":")[1].strip() for line in res_content.split("\n") if line.startswith("###")])

    if out is not None: 
        likely_handle.close()
        most_handle.close()
        
    celltype_dic = {k:celltype_ls[idx] for idx, k in enumerate(gene_dic.keys())}
    return celltype_dic

