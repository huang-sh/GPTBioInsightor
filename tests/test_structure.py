from gptbioinsightor import structure


def test_extract_score(provider, model):
    text_txt = """
    ```thinking
    here is thinking content 
    ```
    B cell: 100  
    T cell: 25  
    cancer cell: 0
    """
    score_dic = structure.extract_score(text_txt, provider, model, None)
    assert score_dic.score_ls[0].score == 100
