import os
from gptbioinsightor import utils
import scanpy as sc


def test_get_marker_from_seurat():
    gene_marker = utils.get_marker_from_seurat("tests/data/pbmc.markers.fil.csv")
    assert list(gene_marker.keys()) == [0, 1, 2, 3, 4, 5, 6, 7, 8]


def test_parse_model():
    assert ("openai", "gpt-4o", "https://api.openai.com/v1/") == utils.parse_api(
        "openai", "gpt-4o", "https://api.openai.com/v1/"
    )
    assert ("openai", "gpt-4o", "https://api.openai.com/v1/") == utils.parse_api(
        None, "openai:gpt-4o", "https://api.openai.com/v1/"
    )
    assert (
        "openrouter",
        "openai/gpt-4o:beta",
        "https://openrouter.ai/api/v1",
    ) == utils.parse_api(
        None, "openrouter:openai/gpt-4o:beta", "https://openrouter.ai/api/v1"
    )


def test_get_gene_dict():
    topnumber = 15
    rm_gene = True
    key = "rank_genes_groups"
    groups = None
    adata = sc.read_h5ad("tests/data/test.h5ad")
    gene_dic = utils.get_gene_dict(adata, groups, key, topnumber, rm_gene)
    assert list(gene_dic.keys()) == ["0", "1", "2", "3", "4", "5", "6", "7"]
    assert len(gene_dic["0"]) == topnumber

    groups = ["0", "1"]
    gene_dic = utils.get_gene_dict(adata, groups, key, topnumber, rm_gene)
    assert list(gene_dic.keys()) == ["0", "1"]
    assert len(gene_dic["0"]) == topnumber

    groups = ["0", "1"]
    rm_gene = False
    gene_dic = utils.get_gene_dict(adata, groups, key, topnumber, rm_gene)
    bool_ls = [i.startswith(("MT-", "RPL", "RPS")) for i in gene_dic["0"]]
    assert any(bool_ls)


def test_unify_name(provider, model):
    dic = {"0": {"Natural Killer Cells ": "85"}, "1": "NK Cell: 80"}
    new_dic = utils.unify_name(dic, model, provider=provider, base_url=None)
    assert list(new_dic["0"].keys())[0] == list(new_dic["1"].keys())[0]


def test_set_api_key():
    utils.set_api_key("test_key")
    test_key = os.getenv("API_KEY")
    assert test_key == "test_key"

    utils.set_api_key("openai_api_key", "openai")
    openai_key = os.getenv("OPENAI_API_KEY")
    assert openai_key == "openai_api_key"

    utils.set_api_key("meta_api_key", provider="meta")
    meta_key = os.getenv("META_API_KEY")
    assert meta_key == "meta_api_key"
