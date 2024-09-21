import os
from openai import OpenAI
from .prompt import SYSTEM_PROMPT
from .constant import API_SOURCE


class ApiKeyMissingError(Exception):
    """Exception raised for missing API_KEY."""
    def __init__(self, message="API_KEY is missing"):
        self.message = message
        super().__init__(self.message)


def openai_client(msgs, apikey, model, provider, base_url=None, sys_prompt=True):
    if base_url is None:
        base_url = API_SOURCE[provider]

    client = OpenAI(
        api_key= apikey, 
        base_url=base_url,
    )
    sys_msg = [{"role": "system", "content": SYSTEM_PROMPT}] if sys_prompt else []
    response = client.chat.completions.create(
        model=model,
        messages=sys_msg + msgs
    )
    return response.choices[0].message.content


def anthropic_client(msgs, model, apikey, sys_prompt=True):
    import anthropic

    client = anthropic.Anthropic(
        api_key=apikey
    )
    sys_msg = [{"role": "system", "content": SYSTEM_PROMPT}] if sys_prompt else []
    response = client.messages.create(
        model=model,
        max_tokens=8192, # 
        messages=msgs
    )
    content = response.content[0].text
    return content


def query_model(msgs, provider, model, base_url=None, sys_prompt=True):
    API_KEY = os.getenv("API_KEY")
    if API_KEY is None:
        raise ApiKeyMissingError("Note: API key not found, please set API_KEY")
    
    if provider == "anthropic":
        content = anthropic_client(msgs, model, API_KEY, sys_prompt=sys_prompt)
    else:
        content = openai_client(msgs, API_KEY, model, provider, base_url=base_url, sys_prompt=sys_prompt)
    return content
