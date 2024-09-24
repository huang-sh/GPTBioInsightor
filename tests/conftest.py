import os
import pytest


def pytest_addoption(parser):
    parser.addoption("--provider",  default="siliconflow", help="model provider")
    parser.addoption("--model",  default="Qwen/Qwen2.5-72B-Instruct", help="model selection")


@pytest.fixture()
def provider(request):
    value = request.config.getoption("--provider")    
    API_KEY = os.getenv(f"{value.upper()}_API_KEY")
    if API_KEY is None:
        API_KEY = os.getenv("API_KEY")
    os.environ['API_KEY'] = API_KEY
    return value


@pytest.fixture()
def model(request):
    return request.config.getoption("--model")


@pytest.fixture(autouse=True, scope="session")
def set_global_variable():
    import builtins
    builtins.MODEL_LS = ["qwen2-72b-instruct", "gpt-4o", "Qwen/Qwen2.5-72B-Instruct"]
