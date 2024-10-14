# 模型选择

当前GPTBioinsightor支持多个大语言模型API：
- openai的ChatGPT
- anthropic的Claude
- 阿里云的通义千问
- 智谱清言 
- 月之暗面
- 深度求索
- siliconflow
- groq

GPTBioinsightor 提供3个参数来设置模型的选择：
- `provider`, 设置模型提供者，例如: "openai", "aliyun", "anthropic"
- `model`, 具体模型
- `base_url` 自定义模型API接口


## ChatGPT

示例用法：
```python
### 设置ChatGPT的 API KEY
import os
os.environ['OPENAI_API_KEY'] = "sk-***"

import gptbioinsightor as gbi


# 设置数据的背景信息
background = "Cells are PBMCs from a Healthy Donor" 

res = gbi.gpt_celltype(adata, background=background, out="gbi.md", provider="openai", model="gpt-4o-mini")

```
API KEY 获取方式：https://platform.openai.com/api-keys
更多模型参考：https://platform.openai.com/docs/models/


## 通义千问

示例用法：
```python
### 设置通义千问的 API KEY
import os
os.environ['ALIYUN_API_KEY'] = "sk-***"

import gptbioinsightor as gbi

# 设置数据的背景信息
background = "Cells are PBMCs from a Healthy Donor" 

res = gbi.gpt_celltype(adata, background=background, out="gbi.md", provider="aliyun", model="qwen2-72b-instruct")

```
API KEY 获取方式：https://help.aliyun.com/zh/model-studio/developer-reference/get-api-key
更多模型参考：https://help.aliyun.com/zh/dashscope/developer-reference/tongyi-qianwen-7b-14b-72b-api-detailes


## Claude

示例用法：

```python
### 设置月之暗面的 API KEY
import os
os.environ['ANTHROPIC_API_KEY'] = "sk-***"

import gptbioinsightor as gbi

# 设置数据的背景信息
background = "Cells are PBMCs from a Healthy Donor" 

res = gbi.gpt_celltype(adata, background=background, out="gbi.md", provider="anthropic", model="claude-3-5-sonnet-20240620")

```
API KEY 获取方式：https://console.anthropic.com/settings/keys
更多模型参考：https://console.anthropic.com/settings/keys

## 自定义API URL

自定义API URL 需要满足 openai 的接口规范
```
### 设置API KEY
import os
os.environ['API_KEY'] = "sk-***"

# https://api.flybirdsci.com
BASE_URL = "https://api.flybirdsci.com/v1"

background = "Cells are PBMCs from a Healthy Donor" 

res = gbi.gpt_celltype(adata, background=background, out="gbi.md", base_url=BASE_URL, model="gpt-4o")

```
