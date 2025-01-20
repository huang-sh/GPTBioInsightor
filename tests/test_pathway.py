import pytest
from gptbioinsightor import name_pathway, analyse_pathway

MODEL_LIST = pytest.MODEL_LIST

def test_name_pathway(provider, model):
    gene_dic = {
        "gs1": ["CD14", "S100A8", "S100A9", "LYZ", "FCN1", " TYROBP"],
        "gs2": ["FCGR3A", "MS4A7", "CDKN1C", "CKB", "LILRA3", "IFITM3"]
        }
    name_pathway(gene_dic,provider=provider, model=model)
    analyse_pathway(gene_dic,provider=provider, model=model)
