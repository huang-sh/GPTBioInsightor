# Model selection

Currently, GPTBioinsightor supports three large language models:
- ChatGPT by OpenAI
- Qwen by Aliyun
- Moonshot

GPTBioinsightor provides three parameters to select model :：
- `provider`: set the model provider, which can be: "openai", "aliyun", "moonshot"
- `model`: the specific model to use
- `base_url`: for customizing the model's API

## ChatGPT

demo：
```python
### set ChatGPT API KEY
import os
os.environ['API_KEY'] = "sk-***"

import gptbioinsightor as gbi

background = "Cells are PBMCs from a Healthy Donor" 

res = gbi.get_celltype(adata, background=background, out="gbi.md", provider="openai", model="gpt-4o-mini")

```

For more models, please refer to: https://platform.openai.com/docs/models/


## Qwen

demo:
```python
### Qwen API KEY
import os
os.environ['API_KEY'] = "sk-***"

import gptbioinsightor as gbi

background = "Cells are PBMCs from a Healthy Donor" 

res = gbi.get_celltype(adata, background=background, out="gbi.md", provider="aliyun", model="qwen2-72b-instruct")

```

For more models, please refer to: https://help.aliyun.com/zh/dashscope/developer-reference/tongyi-qianwen-7b-14b-72b-api-detailes


## Moonshot

demo：

```python
### Moonshot API KEY
import os
os.environ['API_KEY'] = "sk-***"

import gptbioinsightor as gbi

background = "Cells are PBMCs from a Healthy Donor" 

res = gbi.get_celltype(adata, background=background, out="gbi.md", provider="moonshot", model="moonshot-v1-8k")

```

three models：
- moonshot-v1-8k
- moonshot-v1-32k
- moonshot-v1-128k

## Custom API URL

The custom API URL needs to comply with the OpenAI interface specification.
```
### set API KEY
import os
os.environ['API_KEY'] = "sk-***"

# https://api.flybirdsci.com
BASE_URL = "https://api.flybirdsci.com/v1"

background = "Cells are PBMCs from a Healthy Donor" 

res = gbi.get_celltype(adata, background=background, out="gbi.md", base_url=BASE_URL, model="gpt-4o")

```
