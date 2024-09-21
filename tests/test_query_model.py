from gptbioinsightor.core import query_model


def test_aliyun_query_model(request):
    model = request.config.getoption("--model", default="qwen2-1.5b-instruct")
    msg = [{"role": "user", "content": "Are you a LLM? Just answer yes or no."}]
    resp = query_model(msg, provider='aliyun', model=model, sys_prompt=False)
    if model in ["qwen2-72b-instruct", "gpt-4o"]:
        assert "YES" in resp.upper()
