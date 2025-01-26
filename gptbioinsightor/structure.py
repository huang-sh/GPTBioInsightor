from pydantic import BaseModel, Field
import instructor
from typing import List, Tuple
from litellm import completion

from .utils import get_api_key, parse_api


class ScoreInfo(BaseModel):
 
    score_thinking: str = Field(
        description="the pure thinking content celltype score, it may be within ```thinking``` or tag <thinking></thinking>, or under the 'thinking' word or header"
    )
 
    score_ls: List[Tuple[str, float]] = Field(
        ...,
        description="A list where each element is a Tuple. Tuple is like (celltype, score)"
    )

def extract_score(content, provider, model, base_url):
    provider, model, base_url = parse_api(provider, model, base_url)
    API_KEY = get_api_key(provider)

    client = instructor.from_litellm(completion )
    
    resp = client.chat.completions.create(
        model=model,
        api_key=API_KEY,
        api_base=base_url,
        max_tokens=4000,
        messages=[
            {
                "role": "user",
                "content": content, 
            }
        ],
        response_model=ScoreInfo,
    )
    return resp
