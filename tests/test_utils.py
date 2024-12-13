from gptbioinsightor import utils
import scanpy as sc

def test_get_marker_from_seurat():
    gene_marker = utils.get_marker_from_seurat("tests/data/pbmc.markers.fil.csv")
    assert list(gene_marker.keys()) == [0,1,2,3,4,5,6,7,8]


def test_parse_model():
    assert ("openai", "gpt-4o") == utils.parse_model("openai", "gpt-4o")
    assert ("openai", "gpt-4o") == utils.parse_model(None, "openai:gpt-4o")
    assert ("openrouter", "openai/gpt-4o:beta") == utils.parse_model(None, "openrouter:openai/gpt-4o:beta")

def test_get_gene_dict():
    topnumber = 15
    rm_gene = True 
    key = "rank_genes_groups"
    groups = None
    adata = sc.read_h5ad("tests/data/test.h5ad")
    gene_dic = utils.get_gene_dict(adata, groups, key, topnumber, rm_gene)
    assert list(gene_dic.keys()) == ["0","1","2","3","4","5","6","7"]
    assert len(gene_dic["0"]) == topnumber

    groups = ["0", "1"]
    gene_dic = utils.get_gene_dict(adata, groups, key, topnumber, rm_gene)
    assert list(gene_dic.keys()) ==  ["0", "1"]
    assert len(gene_dic["0"]) == topnumber

    groups = ["0", "1"]
    rm_gene = False
    gene_dic = utils.get_gene_dict(adata, groups, key, topnumber, rm_gene)
    bool_ls = [i.startswith(('MT-', 'RPL', 'RPS')) for i in gene_dic["0"]]
    assert any(bool_ls)
