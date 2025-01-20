import os
import pytest


def pytest_addoption(parser):
    parser.addoption("--provider",  default="deepseek", help="model provider")
    parser.addoption("--model",  default="deepseek-chat", help="model selection")


@pytest.fixture()
def provider(request):
    return request.config.getoption("--provider")


@pytest.fixture()
def model(request):
    return request.config.getoption("--model")


def pytest_configure(config):
    pytest.MODEL_LIST = [
        "Qwen/Qwen2.5-72B-Instruct", "qwen2-72b-instruct", "gpt-4o", "deepseek-chat"
    ]
    
# @pytest.fixture(autouse=True, scope="session")
# def set_global_variable():
#     import builtins
#     builtins.MODEL_LIST = ["qwen2-72b-instruct", "gpt-4o", "Qwen/Qwen2.5-72B-Instruct"]
