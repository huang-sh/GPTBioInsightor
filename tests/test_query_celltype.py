from gptbioinsightor.celltype import _query_celltype, get_celltype, get_subtype
from gptbioinsightor.prompt import SYSTEM_PROMPT
import os
import pytest


base_url=None
MODEL_LIST = pytest.MODEL_LIST

def test_query_celltype(provider, model):
    genes = ["CD19", "MS4A1", "CD79A", "CD79", "CCR7"]
    queryid = "TEST_GENE_SET"
    background = "Human blood"
    content = _query_celltype(genes, queryid, background, provider, model, base_url, SYSTEM_PROMPT)
    if model in MODEL_LIST:
        assert "B cell" in content


def test_get_celltype(provider, model):
    gene_dic = {
        "gs1": ["CD19", "MS4A1", "CD79A", "CD79", "CCR7"],
        "gs2": ["CD3", "CD3E", "CD8", "CD3D", "25CD" ]}
    background = "Human blood"
    celltype_dic = get_celltype(gene_dic, background=background,
                    provider=provider, model=model,base_url=base_url)
    
    celltype_ls = list(celltype_dic.values())
    print(celltype_dic)
    if model in MODEL_LIST:
        assert "B" in celltype_ls[0]
        assert "T" in celltype_ls[1]


def test_get_subtype(provider, model):
    gene_dic = {
        "gs1": ["CD14", "S100A8", "S100A9", "LYZ", "FCN1", " TYROBP"],
        "gs2": ["FCGR3A", "MS4A7", "CDKN1C", "CKB", "LILRA3", "IFITM3"]
        }
    background = "Human blood"
    subtype_dic = get_subtype(gene_dic, background=background, provider=provider, model=model,base_url=base_url)
    if model in MODEL_LIST:
        assert ["gs1", "gs2"] == list(subtype_dic.keys())
