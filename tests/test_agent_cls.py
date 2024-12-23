from gptbioinsightor.core import Agent 
import pytest
from gptbioinsightor.core import query_model

MODEL_LIST = pytest.MODEL_LIST


def test_agent(provider, model):
    agent = Agent(model=model, provider=provider, sys_prompt=None, base_url=None)
    
    response_any1 = agent.query("tell me something", add_context=False)
    response_any2 = agent.query("tell me something", add_context=False)
    assert response_any1 == response_any2

    response_m1 = agent.query("what is your models", add_context=False,use_cache=False)
    response_m2 = agent.query("what is your models", add_context=False,use_cache=False)
    assert response_m1 != response_m2

    response1 = agent.query("Please remember I am Bob")
    response2 = agent.query("what is my name?")
    assert len(agent.history) == 4
    assert "Bob" in response2

    agent.update_context([{'role': 'user', 'content': 'test'}, {'role': 'assistant', 'content': "test1"}])
    assert len(agent.history) == 6
    assert agent.history[-1]["content"] == "test1"

    agent.update_context({'role': 'user', 'content': 'test2'})
    assert len(agent.history) == 7
    assert agent.history[-1]["content"] == "test2"
