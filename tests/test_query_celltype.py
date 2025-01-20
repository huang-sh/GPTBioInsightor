from gptbioinsightor.celltype import get_celltype, get_subtype, get_cellstate
from gptbioinsightor.prompt import SYSTEM_PROMPT
import os
import pytest


base_url=None
MODEL_LIST = pytest.MODEL_LIST

def test_get_celltype(provider, model):
    gene_dic = {
        "gs1": ["CD19", "MS4A1", "CD79A", "CD79", "CCR7"],
        "gs2": ["CD3", "CD3E", "CD8", "CD3D", "25CD" ]}
    background = "Human blood"
    celltype_dic = get_celltype(gene_dic, background=background,
                    provider=provider, model=model,base_url=base_url)
    
    if model in MODEL_LIST:
        assert "B" in list(celltype_dic["gs1"].keys())[0]
        assert "T" in list(celltype_dic["gs2"].keys())[0]


def test_get_subtype(provider, model):
    gene_dic = {
        "gs1": ["CD14", "S100A8", "S100A9", "LYZ", "FCN1", " TYROBP"],
        "gs2": ["FCGR3A", "MS4A7", "CDKN1C", "CKB", "LILRA3", "IFITM3"]
        }
    background = "Human blood"
    subtype_dic = get_subtype(gene_dic, celltype="monocyte",background=background, provider=provider, model=model,base_url=base_url)
    if model in MODEL_LIST:
        assert ["gs1", "gs2"] == list(subtype_dic.keys())


def test_get_cellstate(provider, model):
    gene_dic = {
        "gs1": ["CD14", "S100A8", "S100A9", "LYZ", "FCN1", " TYROBP"],
        "gs2": ["FCGR3A", "MS4A7", "CDKN1C", "CKB", "LILRA3", "IFITM3"]
        }
    background = "Human blood"
    get_cellstate(gene_dic, background=background, provider=provider, model=model)
 
