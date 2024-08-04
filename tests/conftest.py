
def pytest_addoption(parser):
    parser.addoption("--model",  default="qwen2-1.5b-instruct", help="model selection")
