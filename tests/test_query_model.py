import pytest
from gptbioinsightor.core import query_model

MODEL_LIST = pytest.MODEL_LIST

def test_aliyun_query_model(provider, model):
    msg = [{"role": "user", "content": "Are you a LLM? Just answer yes or no."}]
    resp = query_model(msg, provider=provider, model=model, sys_prompt=None)
    print(model, resp)
    if model in MODEL_LIST:
        assert "YES" in resp.upper()
