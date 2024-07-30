from gptbioinsightor import utils


def test_get_marker_from_seurat():
    gene_marker = utils.get_marker_from_seurat("tests/data/pbmc.markers.fil.csv")
    assert list(gene_marker.keys()) == [0,1,2,3,4,5,6,7,8]
