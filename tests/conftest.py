import os
import pytest


def pytest_addoption(parser):
    parser.addoption("--provider",  default="siliconflow", help="model provider")
    parser.addoption("--model",  default="Qwen/Qwen2.5-72B-Instruct", help="model selection")


@pytest.fixture()
def provider(request):
    return request.config.getoption("--provider")


@pytest.fixture()
def model(request):
    return request.config.getoption("--model")


@pytest.fixture(autouse=True, scope="session")
def set_global_variable():
    import builtins
    builtins.MODEL_LIST = ["qwen2-72b-instruct", "gpt-4o", "Qwen/Qwen2.5-72B-Instruct"]
