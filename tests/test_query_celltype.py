from gptbioinsightor.celltype import _query_celltype, gptcelltype


def test_query_celltype():
    genes = ["CD19", "MS4A1", "CD79A", "CD79", "CCR7"]
    queryid = "TEST_GENE_SET"
    background = "Human blood"
    provider = "aliyun"
    model = "qwen2-72b-instruct"
    base_url= None
    sys_prompt = True
    content = _query_celltype(genes, queryid, background, provider, model, base_url, sys_prompt)
    assert "B cell" in content


def test_gptcelltype():
    gene_dic = {
        "gs1": ["CD19", "MS4A1", "CD79A", "CD79", "CCR7"],
        "gs2": ["CD3", "CD3E", "CD8", "CD3D", "25CD" ]}
    background = "Human blood"
    provider = "aliyun"
    model = "qwen2-72b-instruct"
    celltype_dic = gptcelltype(gene_dic, background=background, provider=provider, model=model)
    celltype_ls = list(celltype_dic.values())
    assert "B" in celltype_ls[0]
    assert "T" in celltype_ls[1]
