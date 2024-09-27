from openai import OpenAI
from .prompt import SYSTEM_PROMPT
from .constant import API_SOURCE
from .utils import get_api_key


def openai_client(msgs, apikey, model, provider, base_url=None, sys_prompt=None, tools=None):
    if base_url is None:
        base_url = API_SOURCE[provider]

    client = OpenAI(
        api_key= apikey, 
        base_url=base_url,
    )
    if sys_prompt is None:
        sys_msg = []
    else:
        sys_msg = [{"role": "system", "content": sys_prompt}] 
    response = client.chat.completions.create(
        model=model,
        messages=sys_msg + msgs,
        tools = tools,
        top_p= 0.7,
        temperature= 0.95
    )
    return response.choices[0].message.content


def anthropic_client(msgs, model, apikey, sys_prompt=''):
    import anthropic

    client = anthropic.Anthropic(
        api_key=apikey
    )
    sys_prompt = '' if sys_prompt is None else sys_prompt
    response = client.messages.create(
        model=model,
        system=sys_prompt,
        max_tokens=8192, 
        messages=msgs
    )
    content = response.content[0].text
    return content


def query_model(msgs, provider, model, base_url=None, sys_prompt=None, tools=None):
    API_KEY = get_api_key(provider)
    if provider == "anthropic":
        content = anthropic_client(msgs, model, API_KEY, sys_prompt=sys_prompt)
    else:
        content = openai_client(msgs, API_KEY, model, provider, base_url=base_url, sys_prompt=sys_prompt, tools=tools)
    return content
