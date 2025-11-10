from __future__ import annotations

import os
import sys
from collections import defaultdict
from html import escape
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

import scanpy as sc
import pandas as pd
from anndata import AnnData

from .exception import ApiKeyMissingError
from .constant import API_SOURCE
from .logging_utils import logger


def get_marker_from_seurat(path: str | Path) -> dict:
    """\
    generate a gene dict from Seurat FindAllMarkers output csv file

    Parameters
    ----------
    path : str | Path
        gene marker csv path

    Returns
    -------
    dict
        gene marker dict
    """
    logger.info("Loading marker data from '%s'.", path)
    df = pd.read_csv(path)
    marker_dict = df.groupby("cluster", observed=True)["gene"].agg(list).to_dict()
    logger.info("Loaded markers for %d clusters.", len(marker_dict))
    return marker_dict


def get_gene_dict(input, group, key, topnumber, rm_genes):
    if isinstance(input, AnnData):
        logger.info(
            "Preparing gene dictionary from AnnData using group '%s' and key '%s'.",
            group or "default",
            key,
        )
        deg_df = sc.get.rank_genes_groups_df(input, group=group, key=key)
        gene_dic = {}
        for gid, sdf in deg_df.groupby("group", observed=True):
            gene_dic[gid] = sdf["names"].tolist()
    elif isinstance(input, dict):
        gene_dic = input.copy()
        logger.info("Using existing gene dictionary with %d entries.", len(gene_dic))
    if rm_genes:
        logger.info("Filtering out mitochondrial and ribosomal genes.")
        for k in gene_dic.keys():
            gene_dic[k] = [
                g
                for g in gene_dic[k]
                if not g.startswith(("MT-", "RPL", "RPS", "ENSG"))
            ]
    for k in gene_dic.keys():
        gene_dic[k] = gene_dic[k][:topnumber]
    logger.info("Prepared gene dictionary with %d gene sets.", len(gene_dic))
    return gene_dic


def parse_api(provider, model, base_url):
    if provider == "ollama":
        OLLAMA_HOST = os.getenv("OLLAMA_HOST")
        if OLLAMA_HOST is not None:
            base_url = os.getenv("OLLAMA_HOST")
        else:
            base_url = API_SOURCE[provider]
    elif provider is not None:
        base_url = API_SOURCE[provider]
    else:
        items = model.split(":")
        provider = items[0]
        model = ":".join(items[1:])
        if provider in API_SOURCE:
            base_url = API_SOURCE[provider]
        else:
            base_url = base_url
    return provider, model, base_url


def get_api_key(provider=None):
    if provider is not None:
        if provider == "ollama":
            API_KEY = "ollama"
        API_KEY = os.getenv(f"{provider.upper()}_API_KEY")
        if API_KEY is None:
            API_KEY = os.getenv("API_KEY")
    else:
        API_KEY = os.getenv("API_KEY")
    if API_KEY is None:
        raise ApiKeyMissingError(
            f"Note: API key not found, please set {provider.upper()}_API_KEY or API_KEY"
        )
    return API_KEY


def set_api_key(api_key: str, provider: str | None = None):
    """\
    set api key for different providers

    Parameters
    ----------
    api_key : str
        api key of the LLM provider 
    provider : str | None, optional
        LLM provider, by default None
    """
    if provider is None:
        os.environ["API_KEY"] = api_key
    else:
        os.environ[f"{provider.upper()}_API_KEY"] = api_key


def get_celltype_name(text):
    for line in text.split("\n"):
        if line.startswith("####"):
            try:
                return line.split(":")[1]
            except IndexError:
                print("LLM doesn't output result accroding predefined format")


class Outputor:
    def __init__(self, path: str | Path | None = None) -> None:
        self.path = path
        if self.path is None:
            self.handle = sys.stdout
            logger.info("Streaming analysis output to standard output.")
        else:
            self.handle = open(self.path, "w", encoding="utf-8")
            logger.info("Writing analysis output to '%s'.", self.path)

    def write(self, text):
        print(text, file=self.handle)

    def close(self):
        if self.path is not None:
            self.handle.close()
            logger.info("Saved analysis output to '%s'.", self.path)


# @lru_cache(maxsize=500)
def list_celltype(num, background, provider, model, base_url, sys_prompt):
    from .core import Agent
    from .prompt import (
        PRE_CELLTYPE_PROMPT1,
        PRE_CELLTYPE_PROMPT2,
        PRE_CELLTYPE_MERGE_PROMPT,
    )

    # query_num = 3
    text = PRE_CELLTYPE_PROMPT1.format(num=num, background=background)
    agent = Agent(
        model=model, provider=provider, sys_prompt=sys_prompt, base_url=base_url
    )
    logger.info(
        "Requesting candidate cell types for %d clusters with provider '%s' and model '%s'.",
        num,
        provider or "default",
        model or "default",
    )
    agent.repeat_query(text, n=3, use_context=False)
    chat_msg = [
        {"role": "user", "content": PRE_CELLTYPE_PROMPT2.format(background=background)},
        {"role": "assistant", "content": agent.query(PRE_CELLTYPE_MERGE_PROMPT)},
    ]
    return chat_msg


def agent_pipe(agent, pct_txt, score_prompt: str | None = None, cluster_id=None):
    from .prompt import CELLTYPE_SCORE, CELLTYPE_REPORT

    cluster_label = (
        f"cluster {cluster_id}" if cluster_id is not None else "current gene set"
    )
    logger.info(f"Generating detailed reasoning for {cluster_label}.")
    agent.query(pct_txt, use_context=True, add_context=True, use_cache=True)
    if score_prompt is None:
        score_instruction = CELLTYPE_SCORE
    else:
        from .prompt import USER_DEFINED_CELLTYPE_SCORE

        score_instruction = USER_DEFINED_CELLTYPE_SCORE.replace(
            "Scoring_Criteria_PROMPT", score_prompt
        )
    logger.info(f"Requesting scoring details for {cluster_label}.")
    scores = agent.query(
        score_instruction, use_context=True, add_context=True, use_cache=False
    )
    report_prompt = CELLTYPE_REPORT.format(score=scores)
    logger.info(f"Compiling final cell type report segment for {cluster_label}.")
    agent.query(report_prompt, use_context=True, add_context=True, use_cache=False)
    return agent.get_history(role="assistant")


def ensemble_agent_pipe(
    agent,
    pct_txt,
    score_prompt: str | None = None,
    cluster_id=None,
    model_configs=None,
    return_metadata=False,
):
    from .prompt import CELLTYPE_SCORE, CELLTYPE_REPORT

    cluster_label = (
        f"cluster {cluster_id}" if cluster_id is not None else "current gene set"
    )
    logger.info(f"Generating detailed reasoning for {cluster_label}.")
    original_provider = agent.provider
    original_model = agent.model
    original_base_url = agent.base_url
    metadata = {
        "primary_label": None,
        "reasoning_by_model": [],
        "score_by_model": [],
        "aggregated_scores": None,
    }
    base_label_parts = []
    if original_provider:
        base_label_parts.append(str(original_provider))
    if original_model:
        base_label_parts.append(str(original_model))
    default_label = ":".join(base_label_parts) if base_label_parts else "primary_model"
    metadata["primary_label"] = default_label
    configured_models: list[dict] | None = None
    if model_configs is None:
        reasoning_response = agent.query(
            pct_txt, use_context=True, add_context=True, use_cache=True
        )
        metadata["reasoning_by_model"].append(
            {"label": default_label, "response": reasoning_response}
        )
    else:
        configured_models = []
        reasoning_blocks = []
        if not isinstance(model_configs, (list, tuple)):
            raise TypeError("model_configs must be an iterable of configuration dicts.")
        for idx, cfg in enumerate(model_configs, start=1):
            if not isinstance(cfg, dict):
                raise TypeError(
                    f"Model config at position {idx} must be a dict, got {type(cfg)}."
                )
            model_name = cfg.get("model")
            if model_name is None:
                raise ValueError(
                    f"Missing 'model' entry in model config at position {idx}."
                )
            provider_name = cfg.get("provider", original_provider)
            base_url = cfg.get("base_url", original_base_url)
            label = cfg.get("label")
            if label is None:
                label_parts = []
                if provider_name:
                    label_parts.append(str(provider_name))
                label_parts.append(str(model_name))
                label = ":".join(label_parts) if label_parts else f"model_{idx}"
            agent.provider = provider_name
            agent.model = model_name
            agent.base_url = base_url
            logger.info(
                "Querying model '%s' (provider '%s') for %s.",
                model_name,
                provider_name or "default",
                cluster_label,
            )
            response = agent.query(
                pct_txt, use_context=True, add_context=False, use_cache=False
            )
            reasoning_blocks.append((label, response))
            configured_models.append(
                {
                    "label": label,
                    "provider": provider_name,
                    "model": model_name,
                    "base_url": base_url,
                }
            )
        if not reasoning_blocks:
            raise ValueError("model_configs must contain at least one configuration.")
        agent.provider = original_provider
        agent.model = original_model
        agent.base_url = original_base_url
        combination_prompt = (
            f"{pct_txt}\n\n"
            "The above user task has already been answered independently by multiple models. "
            "Please synthesize a concise, reconciled reasoning by taking the strengths of each response, "
            "highlighting agreements, resolving conflicts, and pointing out any novel insights contributed by individual models.\n\n"
        )
        for label, response in reasoning_blocks:
            combination_prompt += (
                f"### Response from {label}\n{(response or '').strip()}\n\n"
            )
        combination_prompt += (
            "### Instructions\n"
            "1. Merge the complementary strengths of the above responses into a single coherent analysis.\n"
            "2. Explicitly mention any consensus points across models.\n"
            "3. Address conflicting statements, clarifying which interpretation is better supported.\n"
            "4. Note any unique but compelling insights contributed by individual models.\n"
            "5. Keep the reasoning focused and structured for downstream scoring.\n"
            "6. After the merged reasoning, add a `Candidate Roster` section that lists every plausible cell type referenced across the responses (one entry per cell type, no duplicates, descending confidence if possible).\n"
        )
        agent.provider = original_provider
        agent.model = original_model
        agent.base_url = original_base_url
        combined_response = agent.query(
            combination_prompt,
            use_context=False,
            add_context=True,
            use_cache=False,
        )
        agent.provider = original_provider
        agent.model = original_model
        agent.base_url = original_base_url
        agent.create_conversation(pct_txt, combined_response)
        metadata["reasoning_by_model"] = [
            {"label": label, "response": response}
            for label, response in reasoning_blocks
        ]
    score_instruction = CELLTYPE_SCORE if score_prompt is None else score_prompt
    logger.info(f"Requesting scoring details for {cluster_label}.")
    if model_configs is None:
        score_response = agent.query(
            score_instruction, use_context=True, add_context=True, use_cache=False
        )
        scores = score_response
        parsed_entries = []
        try:
            from .structure import extract_score  # local import to avoid cycles
        except Exception as exc:  # pragma: no cover - defensive
            logger.warning(
                "Failed to import extract_score for single-model scoring: %s", exc
            )
            extract_score = None  # type: ignore
        if extract_score is not None and score_response:
            try:
                parsed = extract_score(
                    score_response,
                    original_provider,
                    original_model,
                    original_base_url,
                )
            except Exception as exc:  # pragma: no cover - defensive
                logger.warning(
                    "Failed to parse scores for single-model response: %s", exc
                )
            else:
                parsed_entries = [
                    {"celltype": entry.celltype, "score": float(entry.score)}
                    for entry in parsed.score_ls
                ]
                metadata["aggregated_scores"] = {
                    entry["celltype"]: entry["score"] for entry in parsed_entries
                }
        metadata["score_by_model"].append(
            {
                "label": default_label,
                "response": score_response,
                "scores": parsed_entries,
            }
        )
    else:
        assert configured_models is not None
        score_blocks = []
        try:
            from .structure import extract_score  # local import to avoid cycles
        except Exception as exc:  # pragma: no cover - defensive
            logger.warning(
                "Failed to import extract_score for ensemble averaging: %s", exc
            )
            extract_score = None  # type: ignore
        for cfg in configured_models:
            agent.provider = cfg["provider"]
            agent.model = cfg["model"]
            agent.base_url = cfg["base_url"]
            logger.info(
                "Scoring with model '%s' (provider '%s') for %s.",
                cfg["model"],
                cfg["provider"] or "default",
                cluster_label,
            )
            response = agent.query(
                score_instruction, use_context=True, add_context=False, use_cache=False
            )
            score_blocks.append((cfg, response))
            parsed_entries = []
            if response and extract_score is not None:
                try:
                    parsed = extract_score(
                        response,
                        cfg["provider"] or original_provider,
                        cfg["model"],
                        cfg["base_url"],
                    )
                except Exception as exc:  # pragma: no cover - defensive
                    logger.warning(
                        "Failed to parse scores from model '%s': %s", cfg["label"], exc
                    )
                else:
                    parsed_entries = [
                        {"celltype": entry.celltype, "score": float(entry.score)}
                        for entry in parsed.score_ls
                    ]
            metadata["score_by_model"].append(
                {
                    "label": cfg["label"],
                    "response": response,
                    "scores": parsed_entries,
                }
            )
        agent.provider = original_provider
        agent.model = original_model
        agent.base_url = original_base_url
        aggregated = False
        raw_score_dict = {
            entry["label"]: {
                score["celltype"]: score["score"] for score in entry["scores"]
            }
            for entry in metadata["score_by_model"]
        }
        if raw_score_dict:
            base_cfg = configured_models[0]
            unify_model = original_model or base_cfg["model"]
            unify_provider = original_provider or base_cfg["provider"]
            unify_base_url = original_base_url or base_cfg["base_url"]
            try:
                unified_scores = unify_name(  # type: ignore[name-defined]
                    raw_score_dict,
                    unify_model,
                    provider=unify_provider,
                    base_url=unify_base_url,
                )
            except Exception as exc:  # pragma: no cover - defensive
                logger.warning("Failed to unify cell type names across models: %s", exc)
                unified_scores = raw_score_dict
            aggregated_scores_map: defaultdict[str, list[float]] = defaultdict(list)
            support_counts: defaultdict[str, int] = defaultdict(int)
            for entry in metadata["score_by_model"]:
                updated_scores = unified_scores.get(entry["label"], {})
                entry["scores"] = [
                    {"celltype": name, "score": float(value)}
                    for name, value in updated_scores.items()
                ]
                for name, value in updated_scores.items():
                    aggregated_scores_map[name].append(float(value))
                    support_counts[name] += 1
            if aggregated_scores_map:
                aggregated = True
                averaged_scores = {
                    name: sum(values) / len(values)
                    for name, values in aggregated_scores_map.items()
                }
                metadata["aggregated_scores"] = averaged_scores
                sorted_scores = sorted(
                    (
                        (name, averaged_scores[name], support_counts[name])
                        for name in aggregated_scores_map.keys()
                    ),
                    key=lambda item: item[1],
                    reverse=True,
                )
                summary_lines = [
                    "Ensemble scoring summary (averaged across contributing models):"
                ]
                for idx, (celltype, avg_score, count) in enumerate(
                    sorted_scores, start=1
                ):
                    summary_lines.append(
                        f"CELLTYPE{idx}: {celltype} (ensemble score: {avg_score:.2f}; models: {count})"
                    )
                scores = "\n".join(summary_lines)
            else:
                aggregated = False
        if not aggregated:
            logger.warning(
                "Falling back to concatenated per-model scores because averaging failed."
            )
            scores = "\n\n".join(
                f"### Score from {cfg['label']}\n{(resp or '').strip()}"
                for cfg, resp in score_blocks
            )
        per_model_sections = [
            f"### Score from {cfg['label']}\n{(resp or '').strip()}"
            for cfg, resp in score_blocks
        ]
        if aggregated and per_model_sections:
            scores = f"{scores}\n\n" + "\n\n".join(per_model_sections)
        agent.create_conversation(score_instruction, scores)
    report_prompt = CELLTYPE_REPORT.format(score=scores)
    logger.info(f"Compiling final cell type report segment for {cluster_label}.")
    agent.query(report_prompt, use_context=True, add_context=True, use_cache=False)
    agent.provider = original_provider
    agent.model = original_model
    agent.base_url = original_base_url
    history = agent.get_history(role="assistant")
    if return_metadata:
        return {"history": history, "metadata": metadata}
    return history


def get_score_prompt() -> str:
    """Return the default scoring prompt used for cell type evaluation."""
    from .prompt import CELLTYPE_SCORE

    return CELLTYPE_SCORE


def score_heatmap(score_dic, cutoff=0, figsize=(10, 6), cmap="viridis"):
    import seaborn as sns
    import matplotlib.pyplot as plt

    df = pd.DataFrame(score_dic).T.apply(pd.to_numeric).round(1)
    df = df[df >= cutoff].dropna(axis=1, how="all")
    plt.figure(figsize=figsize)
    base_size = min(figsize) * 2
    font_size = max(base_size / max(df.shape), 8)
    heatmap = sns.heatmap(
        df,
        annot=True,
        cmap=cmap,
        fmt="g",
        linewidths=0.5,
        annot_kws={"size": font_size},
    )
    heatmap.set_xticklabels(
        heatmap.get_xticklabels(), rotation=45, ha="right", fontsize=12
    )
    plt.title("CellType Score Heatmap")
    plt.xlabel("CellTypes")
    plt.ylabel("Cluster")
    return heatmap


def unify_name(dic, model, provider=None, base_url=None):
    from .core import Agent

    agent = Agent(model=model, provider=provider, sys_prompt=None, base_url=base_url)

    fmt_demo = """
    {
        '0': {'Plasma': '85', 'Memory B Cells': '70', 'Activated B Cells': '60'},
        '1': {'T Cells': '85', 'NK Cells': '25', 'Dendritic Cells': '5'}
    }
    """
    correct_txt = """eg. if you meet Dendritic Cells, DC, Dendritic Cell, you should correct them as one same name, such as Dendritic Cell;
    If you meet Platelets/Megakaryocytes, Megakaryocytes/Platelets, you should correct them as one same name, such as Megakaryocytes/Platelets."""
    text = f"""
    ```JSON
    {str(dic)}
    ```
    Unify the cell type names in this JSON data, using the same format and term or name to represent the same cell type.
    {correct_txt}
    Only return the corrected JSON format data,without any additional characters or text, such as "", ```or ', like:

    {fmt_demo}

    """
    logger.info("Harmonizing cell type names provided by the language model.")
    new_dic_str = agent.query(
        text, use_context=False, add_context=False, use_cache=True
    )
    try:
        new_dic = eval(new_dic_str)
    except Exception:
        print("Failed to unify the cell type names")
        new_dic = dic
    else:
        logger.info("Cell type names have been unified successfully.")
    return new_dic


def add_obs(adata, score_dic, add_key="gbi_celltype", cluster_key="leiden"):
    logger.info(
        "Adding predicted cell types to AnnData.obs using key '%s' (source cluster key '%s').",
        add_key,
        cluster_key,
    )
    new_dic = {}
    for key, cell_dict in score_dic.items():
        max_cell = max(cell_dict, key=lambda k: float(cell_dict[k]))
        new_dic[key] = max_cell
    adata.obs[add_key] = adata.obs[cluster_key].map(new_dic)
    return adata


def generate_html_report(report_data, output_path: str | Path):
    """Create a standalone HTML report summarizing cell type predictions."""

    def _render_text_block(text: str | None):
        if not text:
            return '<div class="empty">No additional details.</div>'
        clean = text.strip()
        if not clean:
            return '<div class="empty">No additional details.</div>'
        return f"<pre>{escape(clean)}</pre>"

    def _render_pathways(pathways: dict | None):
        if not pathways:
            return '<div class="empty">No pathway enrichment provided.</div>'
        sections = []
        for db_name, terms in pathways.items():
            term_tags = "".join(
                f'<span class="tag">{escape(str(term))}</span>' for term in terms
            )
            sections.append(
                f'<div class="pathway-block"><strong>{escape(str(db_name))}</strong><div class="tag-wrap">{term_tags}</div></div>'
            )
        return "".join(sections)

    def _render_scores(scores: list | None):
        if not scores:
            return '<div class="empty">No structured scores available.</div>'
        rows = []
        for entry in scores:
            celltype = escape(str(entry.get("celltype", "Unknown")))
            try:
                score_val = float(entry.get("score", 0.0))
            except (TypeError, ValueError):
                score_val = 0.0
            rows.append(f"<tr><td>{celltype}</td><td>{score_val:.1f}</td></tr>")
        return (
            '<table class="score-table"><thead><tr><th>Cell type</th><th>Score</th>'
            "</tr></thead><tbody>" + "".join(rows) + "</tbody></table>"
        )

    def _render_citations(citations):
        if not citations:
            return "<li>No citations provided.</li>"
        if isinstance(citations, str):
            citations = [citations]
        items = []
        for cite in citations:
            if isinstance(cite, dict):
                title = (
                    cite.get("title")
                    or cite.get("lock_title")
                    or cite.get("description")
                    or cite.get("url")
                    or cite.get("source")
                    or "Citation"
                )
                url = cite.get("url") or cite.get("source")
                if url:
                    items.append(
                        f"<li>{escape(str(title))}: "
                        f'<a href="{escape(str(url))}" target="_blank" rel="noopener noreferrer">'
                        f"{escape(str(url))}</a></li>"
                    )
                else:
                    items.append(f"<li>{escape(str(title))}</li>")
            else:
                items.append(f"<li>{escape(str(cite))}</li>")
        return "".join(items)

    def _confidence_label(score: float | None):
        if score is None:
            return "Unknown"
        if score >= 85:
            return "High"
        if score >= 60:
            return "Moderate"
        return "Low"

    clusters = report_data.get("clusters") or []
    nav_items = []
    cards = []
    for idx, cluster in enumerate(clusters):
        cid = cluster.get("cluster_id", idx)
        scores = cluster.get("scores") or []
        if scores:
            try:
                top_score = float(scores[0].get("score", 0.0))
            except (TypeError, ValueError):
                top_score = None
            top_label = scores[0].get("celltype") or "Not provided"
        else:
            top_score = None
            top_label = "Not provided"
        confidence = _confidence_label(top_score)
        score_display = "N/A" if top_score is None else f"{top_score:.0f}%"
        nav_items.append(
            (
                f'<li class="cluster-item{" active" if idx == 0 else ""}" '
                f'data-target="cluster-{idx}">'
                f'<div><div class="cluster-title">Cluster {escape(str(cid))}</div>'
                f'<div class="cluster-label">{escape(str(top_label))}</div></div>'
                f'<span class="cluster-score">{score_display} • {confidence}</span>'
                "</li>"
            )
        )
        genes = cluster.get("genes") or []
        gene_tags = "".join(
            f'<span class="tag">{escape(str(gene))}</span>' for gene in genes
        )
        gene_block = gene_tags or '<div class="empty">No genes provided.</div>'
        citations = cluster.get("citations") or []
        citation_html = _render_citations(citations)
        online_block = cluster.get("online_search") or ""
        cards.append(
            f'<article id="cluster-{idx}" class="cluster-card{" active" if idx == 0 else ""}">'
            f"<header><div><h3>Cluster {escape(str(cid))}</h3>"
            f'<p class="lead">{escape(str(top_label))}</p>'
            f'<span class="badge">Confidence: {confidence}</span></div></header>'
            "<section><h4>Top Genes</h4>"
            f'<div class="tag-wrap">{gene_block}</div></section>'
            "<section><h4>Enrichment Pathway</h4>"
            f"{_render_pathways(cluster.get('pathways'))}</section>"
            "<section><h4>Online Search</h4>"
            f"{_render_text_block(online_block)}"
            "<h5>Citations</h5>"
            f"<ul>{citation_html}</ul></section>"
            "<section><h4>Thinking</h4>"
            f"{_render_text_block(cluster.get('thinking'))}</section>"
            "<section><h4>Score Rationale</h4>"
            f"{_render_text_block(cluster.get('score_text'))}"
            f"{_render_scores(scores)}</section>"
            "<section><h4>Report</h4>"
            f"{_render_text_block(cluster.get('report_text'))}</section>"
            "</article>"
        )

    nav_html = "".join(nav_items) or '<li class="empty">No clusters available.</li>'
    cards_html = "".join(cards) or '<div class="empty">No cluster details.</div>'
    summary_block = _render_text_block(report_data.get("potential_celltype"))
    disclaimer = escape(report_data.get("disclaimer") or "")
    background = escape(report_data.get("background") or "Unspecified")
    online_summary = _render_text_block(report_data.get("online_search_summary"))
    online_section = ""
    if report_data.get("online_search_summary"):
        online_section = (
            '<details class="collapsible">'
            "<summary>Online Search Overview</summary>"
            f"{online_summary}"
            "</details>"
        )

    html_doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>{escape(report_data.get("title") or "GPTBioInsightor Report")}</title>
    <style>
        body {{
            font-family: "Segoe UI", Arial, sans-serif;
            margin: 0;
            background: #f5f6fb;
            color: #1f2937;
        }}
        header.hero {{
            background: linear-gradient(90deg, #0a84ff, #4c6ef5);
            color: #fff;
            padding: 32px;
        }}
        header.hero h1 {{
            margin: 0 0 8px;
        }}
        header.hero p {{
            margin: 4px 0;
        }}
        main {{
            padding: 24px;
        }}
        .summary-card {{
            background: #fff;
            border-radius: 16px;
            padding: 20px;
            box-shadow: 0 8px 30px rgba(15, 23, 42, 0.08);
            margin-bottom: 24px;
            word-break: break-word;
            overflow-wrap: anywhere;
        }}
        details.collapsible {{
            border: 1px solid #e2e8f0;
            border-radius: 12px;
            padding: 12px 16px;
            background: #f8fafc;
        }}
        .summary-card details.collapsible {{
            margin-top: 16px;
        }}
        .summary-card details.collapsible:first-of-type {{
            margin-top: 0;
        }}
        details.collapsible summary {{
            font-weight: 600;
            font-size: 1rem;
            color: #0f172a;
            cursor: pointer;
            outline: none;
            list-style: none;
        }}
        details.collapsible summary::marker,
        details.collapsible summary::-webkit-details-marker {{
            display: none;
        }}
        details.collapsible summary::after {{
            content: "▼";
            float: right;
            transition: transform 0.2s ease;
        }}
        details.collapsible[open] summary::after {{
            transform: rotate(-180deg);
        }}
        .grid {{
            display: grid;
            gap: 24px;
        }}
        .layout {{
            display: grid;
            grid-template-columns: 320px 1fr;
            gap: 24px;
        }}
        .cluster-nav {{
            background: #fff;
            border-radius: 16px;
            box-shadow: 0 8px 30px rgba(15, 23, 42, 0.08);
        }}
        .cluster-nav ul {{
            list-style: none;
            margin: 0;
            padding: 0;
        }}
        .cluster-item {{
            padding: 16px 20px;
            border-bottom: 1px solid #eef2ff;
            cursor: pointer;
            display: flex;
            justify-content: space-between;
            align-items: center;
            transition: background 0.2s ease;
        }}
        .cluster-item.active {{
            background: #eef2ff;
        }}
        .cluster-item:hover {{
            background: #f5f7ff;
        }}
        .cluster-title {{
            font-weight: 600;
            color: #1d4ed8;
        }}
        .cluster-label {{
            font-size: 0.9rem;
            color: #475569;
        }}
        .cluster-score {{
            font-weight: 600;
            color: #0f172a;
        }}
        .cluster-details {{
            display: flex;
            flex-direction: column;
            gap: 16px;
        }}
        .cluster-card {{
            display: none;
            background: #fff;
            border-radius: 16px;
            padding: 24px;
            box-shadow: 0 8px 30px rgba(15, 23, 42, 0.08);
            word-break: break-word;
            overflow-wrap: anywhere;
        }}
        .cluster-card.active {{
            display: block;
        }}
        .badge {{
            display: inline-block;
            background: #fee2e2;
            color: #991b1b;
            padding: 4px 10px;
            border-radius: 999px;
            font-size: 0.85rem;
            margin-top: 8px;
        }}
        section {{
            margin-top: 16px;
        }}
        h4 {{
            margin-bottom: 8px;
            color: #111827;
        }}
        pre {{
            background: #f8fafc;
            padding: 12px;
            border-radius: 8px;
            overflow: auto;
            font-family: "JetBrains Mono", "Fira Code", monospace;
            font-size: 0.9rem;
            white-space: pre-wrap;
            word-break: break-word;
            overflow-wrap: anywhere;
        }}
        .tag-wrap {{
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
        }}
        .tag {{
            display: inline-block;
            background: #eef2ff;
            color: #312e81;
            padding: 4px 10px;
            border-radius: 999px;
            font-size: 0.85rem;
        }}
        .pathway-block {{
            margin-bottom: 12px;
        }}
        .score-table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 12px;
        }}
        .score-table th, .score-table td {{
            border: 1px solid #e2e8f0;
            padding: 8px;
            text-align: left;
        }}
        .empty {{
            color: #94a3b8;
            font-style: italic;
        }}
        @media (max-width: 1024px) {{
            .layout {{
                grid-template-columns: 1fr;
            }}
        }}
    </style>
</head>
<body>
    <header class="hero">
        <h1>{escape(report_data.get("title") or "CellType Analysis")}</h1>
        <p>Background: {background}</p>
        {"<p>" + disclaimer + "</p>" if disclaimer else ""}
    </header>
    <main>
        <div class="summary-card">
            <details class="collapsible">
                <summary>Potential Cell Types</summary>
                {summary_block}
            </details>
            {online_section}
        </div>
        <div class="layout">
            <aside class="cluster-nav">
                <h2 style="padding: 20px 20px 0;">Clusters</h2>
                <ul>{nav_html}</ul>
            </aside>
            <section class="cluster-details">
                {cards_html}
            </section>
        </div>
    </main>
    <script>
        const navItems = document.querySelectorAll('.cluster-item');
        const cards = document.querySelectorAll('.cluster-card');
        navItems.forEach((item) => {{
            item.addEventListener('click', () => {{
                const target = item.getAttribute('data-target');
                navItems.forEach((nav) => nav.classList.remove('active'));
                item.classList.add('active');
                cards.forEach((card) => {{
                    card.classList.toggle('active', card.id === target);
                }});
            }});
        }});
    </script>
</body>
</html>"""

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html_doc, encoding="utf-8")
    logger.info("HTML report saved to '%s'.", output_path)
    return output_path


def _get_perplexity_key():
    api_key = (
        os.getenv("PERPLEXITY_API_KEY")
        or os.getenv("PPLX_API_KEY")
        or os.getenv("PERPLEXITYAI_API_KEY")
    )
    if api_key is None:
        raise ApiKeyMissingError(
            "Note: API key not found, please set PERPLEXITY_API_KEY or PPLX_API_KEY"
        )
    return api_key


def search_celltype(
    background,
    genes,
    *,
    search_model="sonar",
    system_prompt: str | None = None,
    timeout: int = 60,
):
    """Query Perplexity once for a single cluster gene set."""
    if not genes:
        logger.warning("Skipping Perplexity search because no genes were provided.")
        return {"content": "", "citations": [], "completion": None}

    try:
        from perplexity import Perplexity
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise ImportError(
            "The 'perplexity' package is required when using search_model. "
            "Please install it and try again."
        ) from exc

    if isinstance(genes, (list, tuple, set)):
        gene_text = ", ".join(str(g) for g in genes)
    else:
        gene_text = str(genes)

    query_txt = f"""
    <Input>
    <Biological_context>
        {background or "unknown"}
    </Biological_context>
    <Gene_markers>
        {gene_text}
    </Gene_markers>
    </Input>

    <Task>
    Identify the specific cell type(s) within the provided biological context
    that show high expression of these gene markers, focusing on the most likely cell type and state.
    </Task>
    """.strip()

    _get_perplexity_key()
    messages = []
    if system_prompt:
        messages.append({"role": "system", "content": system_prompt})
    messages.append({"role": "user", "content": query_txt})
    client = Perplexity()

    try:
        completion = client.chat.completions.create(
            messages=messages,
            model=search_model,
            timeout=timeout,
        )
    except TypeError as exc:
        if "timeout" in str(exc):
            completion = client.chat.completions.create(
                messages=messages,
                model=search_model,
            )
        else:
            logger.exception("Perplexity request failed: %s", exc)
            return {
                "content": "",
                "citations": [],
                "completion": {"error": str(exc)},
            }
    except Exception as exc:  # pragma: no cover - defensive
        logger.exception("Perplexity request failed: %s", exc)
        return {
            "content": "",
            "citations": [],
            "completion": {"error": str(exc)},
        }

    def _get_attr(obj, key, default=None):
        if obj is None:
            return default
        if isinstance(obj, dict):
            return obj.get(key, default)
        return getattr(obj, key, default)

    choices = _get_attr(completion, "choices", []) or []
    if not choices:
        logger.warning("Perplexity response did not contain any choices.")
        return {"content": "", "citations": [], "completion": completion}

    first_choice = choices[0]
    message = _get_attr(first_choice, "message", None)
    content = ""
    if message:
        raw_content = _get_attr(message, "content", "")
        if isinstance(raw_content, list):
            content = " ".join(str(part) for part in raw_content)
        else:
            content = str(raw_content)
    else:
        content = str(_get_attr(first_choice, "content", "") or "")
    content = content.strip()

    citations = (
        _get_attr(message, "citations", None)
        or _get_attr(first_choice, "citations", None)
        or _get_attr(completion, "citations", None)
        or []
    )
    if citations and not isinstance(citations, list):
        citations = [citations]
    logger.debug("Received Perplexity content with %d citations.", len(citations))
    return {"content": content, "citations": citations, "completion": completion}


def search_celltype_for_clusters(
    background,
    gene_dict,
    *,
    search_model,
    system_prompt: str | None = None,
    max_workers: int | None = None,
    timeout: int = 60,
):
    """\
    Run Perplexity search for every cluster in parallel and build a short summary.

    Returns a mapping containing per-cluster search outputs and a concise summary
    that can be reused in downstream prompts.
    """
    if not gene_dict:
        return {"details": {}, "summary": "", "messages": []}

    if max_workers is None:
        cpu_count = os.cpu_count() or 1
        max_workers = min(max(1, cpu_count // 2), len(gene_dict))
        max_workers = max(max_workers, 1)

    logger.info(
        "Running Perplexity (%s) search for %d clusters with %d workers.",
        search_model,
        len(gene_dict),
        max_workers,
    )
    details = {}
    futures = {}
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for cluster_id, genes in gene_dict.items():
            futures[
                executor.submit(
                    search_celltype,
                    background,
                    genes,
                    search_model=search_model,
                    system_prompt=system_prompt,
                    timeout=timeout,
                )
            ] = cluster_id
        for future in as_completed(futures):
            cluster_id = futures[future]
            try:
                details[cluster_id] = future.result()
            except Exception as exc:  # pragma: no cover - defensive
                logger.exception(
                    "Perplexity search failed for cluster %s: %s", cluster_id, exc
                )
                details[cluster_id] = {
                    "content": "",
                    "citations": [],
                    "completion": {"error": str(exc)},
                }

    summary_lines = []
    for cluster_id in sorted(details, key=str):
        content = details[cluster_id].get("content") or "No response."
        summary_lines.append(f"cluster {cluster_id}: {content}")
    summary = "\n".join(summary_lines)
    messages = [
        {
            "role": "user",
            "content": (
                f"Summarize potential cell types for {len(gene_dict)} clusters in "
                f"the context of '{background or 'unspecified'}'."
            ),
        },
        {"role": "assistant", "content": summary},
    ]
    return {"details": details, "summary": summary, "messages": messages}
