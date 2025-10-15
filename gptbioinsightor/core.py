import os
from functools import lru_cache, partial
from concurrent.futures import ThreadPoolExecutor


from .utils import get_api_key, parse_api
from .exception import APIStatusError, ApiBalanceLow


def openai_client(msgs, apikey, model, base_url=None, sys_prompt=None, temperature=0.5):
    from openai import OpenAI

    client = OpenAI(api_key=apikey, base_url=base_url)
    if model.startswith("o1"):
        sys_prompt = None
        kwargs = {}
    else:
        kwargs = {"top_p": 0.5, "temperature": temperature}
    if sys_prompt is None:
        query_msgs = msgs
    else:
        query_msgs = [{"role": "system", "content": sys_prompt}] + list(msgs)
    try:
        response = client.chat.completions.create(
            model=model, messages=query_msgs, **kwargs
        )
        if response.choices[0].finish_reason != "stop":
            print("finish_reason: ", response.choices[0].finish_reason)
    except APIStatusError as e:
        raise ApiBalanceLow(base_url, e.message, e.response, e.body)

    return response.choices[0].message.content


def anthropic_client(msgs, model, apikey, sys_prompt=None, temperature=0.5):
    import anthropic

    client = anthropic.Anthropic(api_key=apikey)
    if sys_prompt is None:
        sys_prompt = ""
    else:
        sys_prompt = [
            {"type": "text", "text": sys_prompt, "cache_control": {"type": "ephemeral"}}
        ]

    if model in [
        "claude-3-5-sonnet-20241022",
        "claude-3-5-sonnet-20240620",
        "claude-3-5-haiku-20241022",
    ]:
        max_tokens = 8192
    else:
        max_tokens = 4096
    response = client.messages.create(
        model=model,
        system=sys_prompt,
        max_tokens=max_tokens,
        temperature=temperature,
        messages=msgs,
    )
    if response.stop_reason != "end_turn":
        print("stop_reason: ", response.stop_reason)
    content = response.content[0].text
    return content


def query_model(
    msgs, provider, model, base_url=None, sys_prompt=None, temperature=None
):
    provider, model, base_url = parse_api(provider, model, base_url)
    API_KEY = get_api_key(provider)
    if provider == "anthropic":
        content = anthropic_client(msgs, model, API_KEY, sys_prompt=sys_prompt)
    else:
        content = openai_client(
            msgs,
            API_KEY,
            model,
            base_url=base_url,
            sys_prompt=sys_prompt,
            temperature=temperature,
        )
    return content


## source: https://stackoverflow.com/a/66729248
def deep_freeze(thing):
    from collections.abc import Collection, Mapping, Hashable
    from frozendict import frozendict

    if thing is None or isinstance(thing, str):
        return thing
    elif isinstance(thing, Mapping):
        return frozendict({k: deep_freeze(v) for k, v in thing.items()})
    elif isinstance(thing, Collection):
        return tuple(deep_freeze(i) for i in thing)
    elif not isinstance(thing, Hashable):
        raise TypeError(f"unfreezable type: '{type(thing)}'")
    else:
        return thing


def deep_freeze_args(func):
    import functools

    @functools.wraps(func)
    def wrapped(*args, **kwargs):
        return func(*deep_freeze(args), **deep_freeze(kwargs))

    return wrapped


class Agent:
    def __init__(self, model=None, provider=None, sys_prompt=None, base_url=None):
        self.provider = provider
        self.model = model
        self.sys_prompt = sys_prompt
        self.base_url = base_url
        self.history = []
        self._query_with_cache = deep_freeze_args(lru_cache(maxsize=500)(query_model))
        self._query_without_cache = query_model

    def configure(self, provider=None, model=None, base_url=None):
        """Update runtime configuration for downstream API calls."""
        if provider is not None:
            self.provider = provider
        if model is not None:
            self.model = model
        if base_url is not None:
            self.base_url = base_url

    def query(
        self, text, use_context=True, add_context=True, use_cache=True, temperature=1
    ):
        if use_context:
            query_msg = self.history.copy()
        else:
            query_msg = []
        query_msg.append({"role": "user", "content": text})
        if use_cache:
            _query_func = self._query_with_cache
        else:
            _query_func = self._query_without_cache

        response = _query_func(
            query_msg,
            self.provider,
            self.model,
            base_url=self.base_url,
            sys_prompt=self.sys_prompt,
            temperature=temperature,
        )
        if add_context:
            self.history.extend(
                [
                    {"role": "user", "content": text},
                    {"role": "assistant", "content": response},
                ]
            )

        return response

    def _multi_query(
        self,
        texts,
        n_jobs=None,
        use_context=True,
        add_context=True,
        use_cache=True,
        temperature=1,
    ):
        if n_jobs is None:
            n_jobs = min(max(1, os.cpu_count() // 2), len(texts))

        query_func = partial(
            self.query,
            use_context=use_context,
            add_context=False,
            use_cache=use_cache,
            temperature=temperature,
        )
        with ThreadPoolExecutor(max_workers=n_jobs) as executor:
            results = list(executor.map(query_func, texts))
        if add_context:
            for text, res in zip(texts, results):
                self.create_conversation(text, res)
        return results

    def repeat_query(
        self, text, n=1, n_jobs=None, use_context=True, add_context=True, temperature=1
    ):
        if n_jobs is None:
            n_jobs = min(max(1, os.cpu_count() // 2), n)
        texts = [text] * n
        return self._multi_query(
            texts,
            n_jobs=n_jobs,
            use_context=use_context,
            add_context=add_context,
            use_cache=False,  # We don't want to get same results when using repeated query
            temperature=temperature,
        )

    def multi_query(
        self,
        texts,
        n_jobs=None,
        use_context=True,
        add_context=True,
        use_cache=True,
        temperature=1,
    ):
        return self._multi_query(
            texts,
            n_jobs=n_jobs,
            use_context=use_context,
            add_context=add_context,
            use_cache=use_cache,
            temperature=temperature,
        )

    def update_context(self, message):
        if isinstance(message, dict):
            self.history.append(message)
        elif isinstance(message, list):
            self.history.extend(message)

    def get_history(self, role=None):
        if role is None:
            return self.history
        elif role in ["user", "assistant"]:
            return [msg["content"] for msg in self.history if msg["role"] == role]

    def create_conversation(self, user_txt, assistant_txt):
        if self.provider == "anthropic":
            conversation = [
                {
                    "role": "user",
                    "content": [
                        {
                            "type": "text",
                            "text": user_txt,
                            "cache_control": {"type": "ephemeral"},
                        }
                    ],
                },
                {
                    "role": "assistant",
                    "content": [
                        {
                            "type": "text",
                            "text": assistant_txt,
                            "cache_control": {"type": "ephemeral"},
                        }
                    ],
                },
            ]
        else:
            conversation = [
                {"role": "user", "content": user_txt},
                {"role": "assistant", "content": assistant_txt},
            ]
        self.history.extend(conversation)
        return conversation
