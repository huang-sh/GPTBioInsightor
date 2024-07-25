from gptbioinsightor.core import query_model

def test_aliyun_query_model():
    msg = [{"role": "user", "content": "Are you a LLM of Aliyun? Just answer yes or no."}]
    resp = query_model(msg, provider='aliyun', model='qwen2-72b-instruct', sys_prompt=False)
    assert "YES" in resp.choices[0].message.content.upper()