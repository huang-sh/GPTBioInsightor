import os
from openai import OpenAI
from .prompt import SYSTEM_PROMPT
from .constant import API_SOURCE

class ApiKeyMissingError(Exception):
    """Exception raised for missing API_KEY."""
    def __init__(self, message="API_KEY is missing"):
        self.message = message
        super().__init__(self.message)


def query_model(msgs, provider='openai', model='gpt-4o', base_url=None, sys_prompt=True):
    if base_url is None:
        base_url = API_SOURCE[provider]
        
    API_KEY = os.getenv("API_KEY")
    if API_KEY is None:
        raise ApiKeyMissingError("Note: API key not found, please set API_KEY")
        
    client = OpenAI(
        api_key= API_KEY, 
        base_url=base_url,
    )
    sys_msg = [{"role": "system", "content": SYSTEM_PROMPT}] if sys_prompt else []
    response = client.chat.completions.create(
        model=model,
        messages=sys_msg + msgs
    )
    return response

