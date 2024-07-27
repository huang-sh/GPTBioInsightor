import os, sys
from concurrent.futures import ThreadPoolExecutor

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
    input, 
    out=None, 
    background=None, 
    group=None, 
    key="rank_genes_groups", 
    topgenes=15, 
    n_jobs=None, 
    rm_genes=True, 
    provider="qwen", 
    model=None, 
    base_url=None, 
    sys_prompt=True
):
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
        outhandle = sys.stdout
    else:
        outhandle = open(out, "w")
        print("# All possible celltypes", file=outhandle)
        
    if n_jobs is None:
        n_jobs = min(os.cpu_count()//2, len(gene_dic))

    def _aux_func(item):
        return _query_celltype(item[1][:topgenes], item[0], background, provider, model, base_url, sys_prompt)
        
    text = ""
    with ThreadPoolExecutor(max_workers=n_jobs) as executor:
        results = executor.map(_aux_func, gene_dic.items())
        for (gsid, genes), res in zip(gene_dic.items(), results):
            res = res.strip("```").strip("'''")
            print(res, file=outhandle)
            text += res

    msgs = [{"role": "user", "content": text},
            {"role": "user", "content": FINAL_CELLTYPE_PROMPT.format(geneset_num=len(gene_dic), background=background)}]
    response = query_model(msgs, provider=provider, model=model, base_url=base_url, sys_prompt=sys_prompt)
    print("\n\n# Most Possible celltypes", file=outhandle)
    res_content = response.choices[0].message.content.strip("```").strip("'''")
    print(res_content, file=outhandle)
    if out is not None: 
        outhandle.close()
    celltypes = [line.split(":")[1].strip() for line in res_content.split("\n") if line.startswith("###")]
    celltype_dic = {k:celltypes[idx] for idx, k in enumerate(gene_dic.keys())}
    return celltype_dic

