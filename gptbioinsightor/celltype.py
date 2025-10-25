from __future__ import annotations

import os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from collections.abc import Iterable, Mapping
from typing import Any

from anndata import AnnData

from .core import query_model, Agent
from . import utils as ul
from .prompt import (
    SYSTEM_PROMPT,
    CELLTYPE_PROMPT,
    SUBTYPE_PROMPT,
    CHECK_TYPE_PROMPT,
    CELLSTATE_PROMPT,
)
from .structure import extract_score
from .logging_utils import logger


def get_celltype(
    input: AnnData | dict,
    out: Path | str = None,
    background: str = None,
    pathway: dict | None = None,
    key: str = "rank_genes_groups",
    topnumber: int = 15,
    n_jobs: int | None = None,
    provider: str | None = None,
    model: str | None = None,
    search_model: str | None = None,
    group: str | Iterable[str] | None = None,
    base_url: str | None = None,
    rm_genes=True,
    score_prompt: str | None = None,
) -> dict:
    """\
    Annotating genesets using LLM, providing cell types, supporting gene markers, reasons, and potential cell state annotations.

    Parameters
    ----------
    input : AnnData | dict
        An AnnData object or geneset dict
    out : Path | str, optional
        output path, by default None
    background : str, optional
        background information of input data, by default None
    key : str, optional
        rank_genes_groups key, by default "rank_genes_groups"
    topnumber : int, optional
        select top gene for analysis, by default 15
    n_jobs : int | None, optional
        set multiple jobs for querying LLM, by default None
    provider : str| None, optional
        LLM provider, by default None
        "openai" for chatgpt
        "aliyun" for qwen
        "deepseek" for DeepSeek
        "anthropic" for claude
    model : str | None, optional
        set a model based on LLM provider, by default None
    search_model : str | None, optional
        If provided, run Perplexity-based cell type search using this model
        before the main annotation workflow, by default None.
    group : str | Iterable, optional
         Which group, by default None
    base_url : str | None, optional
        customized LLM API url, by default None
    rm_genes : bool, optional
        remove rb and mt genes, by default True
    score_prompt : str | None, optional
        custom scoring criteria prompt to override CELLTYPE_SCORE, by default None

    Returns
    -------
    dict
        a celltypes dict
    """
    sys_prompt = SYSTEM_PROMPT
    logger.info("Starting cell type identification workflow.")
    gene_dic = ul.get_gene_dict(input, group, key, topnumber, rm_genes)
    logger.info("Prepared %d cluster gene sets for analysis.", len(gene_dic))
    genes_ls = []
    for ck, genes in gene_dic.items():
        genes_ls.append(f"   - cluster {ck}: {','.join(genes[:topnumber])}")
    # all_gene_txt = "\n".join(genes_ls)
    search_results = None
    if search_model:
        search_results = ul.search_celltype_for_clusters(
            background,
            gene_dic,
            search_model=search_model,
            max_workers=n_jobs,
        )
        summary_agent = Agent(
            model=model, provider=provider, sys_prompt=sys_prompt, base_url=base_url
        )
        perplexity_summary = search_results.get("summary", "")
        detail_lines = []
        for cluster_id in sorted(search_results["details"], key=str):
            content = search_results["details"][cluster_id].get("content") or ""
            detail_lines.append(f"Cluster {cluster_id}: {content}")
        detail_block = "\n".join(detail_lines)
        summary_prompt = (
            "You are assisting with scRNA-Seq cluster annotation.\n"
            f"Background: {background or 'unspecified'}.\n"
            "External evidence suggests the following potential identities for each cluster:\n"
            f"{perplexity_summary}\n\n"
            "Provide a concise, well-structured overview of the most likely cell types per cluster, "
            "highlighting the leading hypothesis and key rationale for each cluster in 1-2 sentences."
        )
        try:
            candidate_content = summary_agent.query(
                summary_prompt,
                use_context=False,
                add_context=False,
                use_cache=False,
            ).strip()
        except Exception as exc:  # pragma: no cover - defensive
            logger.exception(
                "Failed to summarize Perplexity search outputs, falling back to raw summary: %s",
                exc,
            )
            candidate_content = perplexity_summary or detail_block
        chat_msg = [
            {"role": "user", "content": summary_prompt},
            {"role": "assistant", "content": candidate_content},
        ]
        logger.info(
            "Leveraging Perplexity search results (%s) for downstream prompts.",
            search_model,
        )
    else:
        chat_msg = ul.list_celltype(
            len(gene_dic), background, provider, model, base_url, sys_prompt
        )
        candidate_content = chat_msg[-1]["content"]
    logger.info(
        "Received initial candidate cell types for %d clusters from provider '%s' using model '%s'.",
        len(gene_dic),
        provider or "default",
        model or "default",
    )
    ot = ul.Outputor(out)
    ot.write("# CellType Analysis")
    ot.write(
        "GPTBioInsightor is powered by AI, so mistakes are possible. Review output carefully before use"
    )
    ot.write("## Potential CellType")
    ot.write(
        f"In scRNA-Seq data background of '{background}', the following Potential CellType to be identified:"
    )
    if not search_model:
        candidate_content = chat_msg[-1]["content"] if chat_msg else ""
    candidate_content = candidate_content or ""
    reminder_text = (
        "Reminder: If the current cluster serves as the broad parent lineage for other clusters "
        "(e.g., multiple clusters roll up to this umbrella cell type), zoom in on this cluster with "
        "more granular reasoning and spell out the evidence that sets it apart."
    )
    if reminder_text not in candidate_content:
        if candidate_content:
            candidate_content = f"{candidate_content.rstrip()}\n{reminder_text}"
        else:
            candidate_content = reminder_text
    ot.write(candidate_content)
    if search_results:
        ot.write("## Online Search Overview")
        ot.write(search_results.get("summary", ""))

    worker_count = n_jobs
    if worker_count is None:
        cpu_count = os.cpu_count() or 1
        worker_count = min(max(1, cpu_count // 2), max(1, len(gene_dic)))

    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        futures = {}
        pathway = {} if pathway is None else pathway
        pathway_txt_dic = {}
        for k, genes in gene_dic.items():
            agent = Agent(
                model=model, provider=provider, sys_prompt=sys_prompt, base_url=base_url
            )
            gene_txt = f"   - cluster {k}: {','.join(genes[:topnumber])}"
            cluster_pathway = pathway.get(k, {})
            pw_txt = ""
            for db, pw in cluster_pathway.items():
                pw_txt += f"    - {db}: {','.join(pw)}\n"
            pathway_txt_dic[k] = pw_txt
            pct_txt = CELLTYPE_PROMPT.format(
                candidate=candidate_content,
                setid=k,
                gene=gene_txt,
                setnum=len(gene_dic),
                background=background,
                pathway=pw_txt,
            )
            logger.info("Infer celltype of cluster %s with %d genes.", k, len(genes))
            future = executor.submit(
                ul.agent_pipe,
                agent,
                pct_txt,
                score_prompt=score_prompt,
                cluster_id=k,
            )
            futures[k] = future
        score_dic = {}
        for k, future in futures.items():
            reps = future.result()
            ot.write(f"## cluster geneset {k}\n")
            ot.write("### Gene List\n")
            ot.write(f"Top genes\n```\n{','.join(gene_dic[k])}\n```\n")
            ot.write(f"enrichment pathway\n```\n{pathway_txt_dic[k]}```\n")
            if search_results:
                detail = search_results["details"].get(k, {})
                ot.write("### Online Search Output\n")
                per_content = (
                    detail.get("content") or "No response returned by Perplexity."
                )
                ot.write(per_content)
                ot.write("### Citations\n")
                citations = detail.get("citations") or []
                if citations:
                    for cite in citations:
                        if isinstance(cite, dict):
                            title = (
                                cite.get("title")
                                or cite.get("lock_title")
                                or cite.get("description")
                                or cite.get("url")
                                or "Citation"
                            )
                            url = cite.get("url") or cite.get("source")
                            if url:
                                ot.write(f"- {title}: {url}")
                            else:
                                ot.write(f"- {title}")
                        else:
                            ot.write(f"- {cite}")
                else:
                    ot.write("No citations provided.")
            ot.write("### celltype thinking\n")
            ot.write(reps[0])
            ot.write("### Score\n")
            ot.write(reps[1])
            ot.write("### Report\n")
            ot.write(reps[2])
            score_ls = extract_score(reps[1], provider, model, base_url).score_ls
            score_dic[k] = {entry.celltype: entry.score for entry in score_ls}
            logger.info("Finished scoring cluster %s.", k)
    score_dic = ul.unify_name(score_dic, model, provider, base_url)
    logger.info("Cell type identification completed successfully.")
    return score_dic


def get_celltype_ensemble(
    input: AnnData | dict,
    models: Iterable[Mapping[str, Any]],
    *,
    out: Path | str | None = None,
    background: str | None = None,
    pathway: dict | None = None,
    key: str = "rank_genes_groups",
    topnumber: int = 15,
    n_jobs: int | None = None,
    group: str | Iterable[str] | None = None,
    rm_genes: bool = True,
    score_prompt: str | None = None,
    leader_model: Mapping[str, Any] | None = None,
    max_workers: int | None = None,
    search_model: str | None = None,
    return_details: bool = False,
) -> dict:
    """\
    Run multi-model cell type annotation using :func:`ensemble_agent_pipe` and aggregate the results.

    Parameters
    ----------
    input : AnnData | dict
        An AnnData object or geneset dict that will be analysed cluster by cluster.
    models : Iterable[Mapping[str, Any]]
        Iterable of configuration mappings. Each mapping must provide at least a ``model`` key and may
        optionally include ``provider``, ``base_url``, ``label``, or ``score_prompt``.
    background : str | None, optional
        Background information shared by all runs, by default None.
    pathway : dict | None, optional
        Pathway enrichment information per cluster, by default None.
    key : str, optional
        Rank genes groups key, by default "rank_genes_groups".
    topnumber : int, optional
        Select top gene number for analysis, by default 15.
    n_jobs : int | None, optional
        Number of worker threads for per-cluster queries (legacy argument, retained for API compatibility).
    group : str | Iterable[str] | None, optional
        Which group, by default None.
    rm_genes : bool, optional
        Remove ribosomal and mitochondrial genes before querying, by default True.
    score_prompt : str | None, optional
        Override scoring prompt for every model. If not provided, a shared ``score_prompt`` from the model
        configs may be used, provided they are consistent.
    leader_model : Mapping[str, Any] | None, optional
        Optional configuration for a leader model that performs consensus summarisation and report generation.
    max_workers : int | None, optional
        Deprecated compatibility argument; retained for API stability and ignored.
    search_model : str | None, optional
        If provided, perform an external Perplexity search to seed the ensemble workflow with a summarised candidate roster.
    return_details : bool, optional
        When True, include per-model raw results in ``raw_results``.

    Returns
    -------
    dict
        A dictionary containing ``combined_scores`` (cluster -> cell type -> mean score),
        ``vote_counts`` (cluster -> cell type -> votes), ``total_models`` (number of contributing models),
        and optionally ``raw_results`` when ``return_details`` is True.
    """
    configs = list(models)
    if not configs:
        raise ValueError("models must contain at least one configuration.")

    allowed_keys = {"provider", "model", "base_url", "label", "out", "score_prompt"}
    normalized_configs: list[dict[str, Any]] = []
    seen_labels: set[str] = set()
    score_prompts_found: set[str] = set()
    for idx, raw_cfg in enumerate(configs):
        if not isinstance(raw_cfg, Mapping):
            raise TypeError(
                "Each model configuration must be a mapping containing at least a 'model' key."
            )
        unknown = set(raw_cfg.keys()) - allowed_keys
        if unknown:
            raise ValueError(
                f"Unsupported keys {unknown} in model configuration at position {idx}."
            )
        model_name = raw_cfg.get("model")
        if model_name is None:
            raise ValueError(
                f"Missing 'model' entry in model configuration at position {idx}."
            )
        provider_name = raw_cfg.get("provider")
        label = raw_cfg.get("label") or f"{provider_name or 'default'}:{model_name}"
        if label in seen_labels:
            label = f"{label}_{idx}"
        seen_labels.add(label)
        normalized_configs.append(
            {
                "provider": provider_name,
                "model": model_name,
                "base_url": raw_cfg.get("base_url"),
                "label": label,
            }
        )
        cfg_score_prompt = raw_cfg.get("score_prompt")
        if cfg_score_prompt is not None:
            score_prompts_found.add(cfg_score_prompt)

    if score_prompt is not None:
        score_instruction = score_prompt
    else:
        if len(score_prompts_found) > 1:
            raise ValueError(
                "Model-specific score_prompt values differ. Please provide a single score_prompt argument."
            )
        score_instruction = (
            next(iter(score_prompts_found)) if score_prompts_found else None
        )

    if leader_model is not None and not isinstance(leader_model, Mapping):
        raise TypeError(
            "leader_model must be a mapping defining at least a 'model' entry."
        )

    if leader_model is not None:
        leader_model_name = leader_model.get("model")
        if leader_model_name is None:
            raise ValueError("leader_model configuration must include a 'model' entry.")
        default_provider = leader_model.get("provider")
        default_model = leader_model_name
        default_base_url = leader_model.get("base_url")
    else:
        default_provider = normalized_configs[0]["provider"]
        default_model = normalized_configs[0]["model"]
        default_base_url = normalized_configs[0]["base_url"]

    sys_prompt = SYSTEM_PROMPT
    logger.info(
        "Running ensemble cell type annotation with %d contributing model(s).",
        len(normalized_configs),
    )
    gene_dic = ul.get_gene_dict(input, group, key, topnumber, rm_genes)
    logger.info("Prepared %d cluster gene sets for ensemble analysis.", len(gene_dic))

    search_results = None
    if search_model:
        search_results = ul.search_celltype_for_clusters(
            background,
            gene_dic,
            search_model=search_model,
            max_workers=max_workers,
        )
        summary_agent = Agent(
            model=default_model,
            provider=default_provider,
            sys_prompt=sys_prompt,
            base_url=default_base_url,
        )
        perplexity_summary = search_results.get("summary", "")
        detail_lines = []
        details = search_results.get("details", {}) or {}
        for cluster_id in sorted(gene_dic, key=str):
            detail = details.get(cluster_id) or details.get(str(cluster_id)) or {}
            content = detail.get("content") or ""
            detail_lines.append(f"Cluster {cluster_id}: {content}")
        detail_block = "\n".join(detail_lines)
        summary_prompt = (
            "You are assisting with scRNA-Seq cluster annotation.\n"
            f"Background: {background or 'unspecified'}.\n"
            "External evidence suggests the following potential identities for each cluster:\n"
            f"{perplexity_summary}\n\n"
            "Provide a concise, well-structured overview of the most likely cell types per cluster, "
            "highlighting the leading hypothesis and key rationale for each cluster in 1-2 sentences."
        )
        try:
            candidate_content = summary_agent.query(
                summary_prompt,
                use_context=False,
                add_context=False,
                use_cache=False,
            ).strip()
        except Exception as exc:  # pragma: no cover - defensive
            logger.exception(
                "Failed to summarize Perplexity search outputs for ensemble run, falling back: %s",
                exc,
            )
            candidate_content = perplexity_summary or detail_block
        chat_msg = [
            {"role": "user", "content": summary_prompt},
            {"role": "assistant", "content": candidate_content},
        ]
    else:
        chat_msg = ul.list_celltype(
            len(gene_dic),
            background,
            default_provider,
            default_model,
            default_base_url,
            sys_prompt,
        )
        candidate_content = chat_msg[-1]["content"]

    candidate_content = candidate_content or ""
    reminder_text = (
        "Reminder: If the current cluster serves as the broad parent lineage for other clusters "
        "(e.g., multiple clusters roll up to this umbrella cell type), zoom in on this cluster with "
        "more granular reasoning and spell out the evidence that sets it apart."
    )
    if reminder_text not in candidate_content:
        if candidate_content:
            candidate_content = f"{candidate_content.rstrip()}\n{reminder_text}"
        else:
            candidate_content = reminder_text

    ot = ul.Outputor(out)
    ot.write("# CellType Analysis (Ensemble)")
    ot.write(
        "GPTBioInsightor is powered by AI, so mistakes are possible. Review output carefully before use"
    )
    ot.write("## Potential CellType")
    ot.write(
        f"In scRNA-Seq data background of '{background}', the following Potential CellType to be identified:"
    )
    ot.write(candidate_content)
    if search_results:
        ot.write("## Online Search Overview")
        ot.write(search_results.get("summary", ""))

    ensemble_configs = [
        {
            "label": cfg["label"],
            "provider": cfg["provider"],
            "model": cfg["model"],
            "base_url": cfg["base_url"],
        }
        for cfg in normalized_configs
    ]

    def _normalize_name(name: str) -> str:
        return "".join(ch.lower() for ch in name if ch.isalnum())

    total_models = len(normalized_configs)
    combined_scores: dict[str, dict[str, float]] = {}
    vote_counts: dict[str, dict[str, int]] = {}
    raw_results: dict[str, dict[str, dict[str, float]]] = {} if return_details else {}

    pathway = {} if pathway is None else pathway
    cluster_order = list(gene_dic.keys())
    if not cluster_order:
        summary = {
            "combined_scores": {},
            "vote_counts": {},
            "total_models": len(normalized_configs),
        }
        if return_details:
            summary["raw_results"] = {}
        return summary

    def _prepare_pathway_text(cluster_pathway: dict) -> str:
        if not cluster_pathway:
            return ""
        lines = []
        for db, pw in cluster_pathway.items():
            lines.append(f"    - {db}: {','.join(pw)}")
        return "\n".join(lines) + ("\n" if lines else "")

    def _process_cluster(cluster_id) -> tuple[str, dict]:
        genes = gene_dic[cluster_id]
        agent = Agent(
            model=default_model,
            provider=default_provider,
            sys_prompt=sys_prompt,
            base_url=default_base_url,
        )
        gene_txt = f"   - cluster {cluster_id}: {','.join(genes[:topnumber])}"
        cluster_pathway = pathway.get(cluster_id, {})
        pw_txt = _prepare_pathway_text(cluster_pathway)
        pct_txt = CELLTYPE_PROMPT.format(
            candidate=candidate_content,
            setid=cluster_id,
            gene=gene_txt,
            setnum=len(gene_dic),
            background=background,
            pathway=pw_txt,
        )
        logger.info(
            "Running ensemble inference for cluster %s with %d genes.",
            cluster_id,
            len(genes),
        )
        result = ul.ensemble_agent_pipe(
            agent,
            pct_txt,
            score_prompt=score_instruction,
            cluster_id=cluster_id,
            model_configs=ensemble_configs,
            return_metadata=True,
        )
        history = result["history"]
        metadata = result["metadata"]
        return cluster_id, {
            "history": history,
            "metadata": metadata,
            "genes": genes,
            "pathway_txt": pw_txt,
        }

    if max_workers is None:
        max_workers = min(8, len(cluster_order))
    max_workers = max(1, max_workers)

    cluster_results: dict[str, dict] = {}
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_map = {
            executor.submit(_process_cluster, cluster_id): cluster_id
            for cluster_id in cluster_order
        }
        for future in future_map:
            cluster_id, payload = future.result()
            cluster_results[str(cluster_id)] = payload

    for cluster_id in cluster_order:
        cluster_key = str(cluster_id)
        payload = cluster_results[cluster_key]
        history = payload["history"]
        metadata = payload["metadata"]
        genes = payload["genes"]
        pw_txt = payload["pathway_txt"]

        ot.write(f"## cluster geneset {cluster_id}\n")
        ot.write("### Gene List\n")
        ot.write(f"Top genes\n```\n{','.join(genes)}\n```\n")
        ot.write(f"enrichment pathway\n```\n{pw_txt}```\n")
        if search_results:
            detail = search_results.get("details", {}).get(
                cluster_id
            ) or search_results.get("details", {}).get(str(cluster_id), {})
            ot.write("### Online Search Output\n")
            per_content = detail.get("content") if detail else None
            ot.write(per_content or "No response returned by Perplexity.")
            ot.write("### Online Search Citations\n")
            citations = detail.get("citations") if detail else None
            if citations:
                for cite in citations:
                    if isinstance(cite, dict):
                        title = (
                            cite.get("title")
                            or cite.get("lock_title")
                            or cite.get("description")
                            or cite.get("url")
                            or "Citation"
                        )
                        url = cite.get("url") or cite.get("source")
                        if url:
                            ot.write(f"- {title}: {url}")
                        else:
                            ot.write(f"- {title}")
                    else:
                        ot.write(f"- {cite}")
            else:
                ot.write("No citations provided.")
        ot.write("### celltype thinking\n")
        if history:
            ot.write(history[0])
        if len(history) > 1:
            ot.write("### Score\n")
            ot.write(history[1])
        if len(history) > 2:
            ot.write("### Report\n")
            ot.write(history[2])

        score_bucket: dict[str, list[float]] = {}
        vote_bucket: dict[str, int] = {}

        per_model_scores: dict[str, dict[str, float]] = {}
        for model_score in metadata.get("score_by_model", []):
            label = model_score.get("label") or "model"
            scores = model_score.get("scores") or []
            per_model_scores[label] = {
                entry["celltype"]: float(entry["score"]) for entry in scores
            }

        if per_model_scores:
            try:
                unified_per_model = ul.unify_name(
                    per_model_scores,
                    default_model,
                    provider=default_provider,
                    base_url=default_base_url,
                )
            except Exception as exc:  # pragma: no cover - defensive
                logger.warning(
                    "Failed to unify names across models for cluster %s: %s",
                    cluster_id,
                    exc,
                )
                unified_per_model = per_model_scores
        else:
            unified_per_model = {}

        for label, scores_dict in unified_per_model.items():
            if return_details:
                raw_results.setdefault(label, {})
                raw_results[label][cluster_key] = {
                    celltype: float(score) for celltype, score in scores_dict.items()
                }
            for entry in metadata.get("score_by_model", []):
                if (entry.get("label") or "model") == label:
                    entry["scores"] = [
                        {"celltype": celltype, "score": float(score)}
                        for celltype, score in scores_dict.items()
                    ]
            for celltype, score in scores_dict.items():
                score_bucket.setdefault(celltype, []).append(float(score))
                vote_bucket[celltype] = vote_bucket.get(celltype, 0) + 1

        if metadata.get("aggregated_scores") and not score_bucket:
            for celltype, value in metadata["aggregated_scores"].items():
                score_bucket.setdefault(celltype, []).append(float(value))
                vote_bucket[celltype] = vote_bucket.get(celltype, 0) + 1

        if score_bucket:
            combined_scores[cluster_key] = dict(
                sorted(
                    (
                        (
                            celltype,
                            sum(values) / total_models,
                        )
                        for celltype, values in score_bucket.items()
                    ),
                    key=lambda item: item[1],
                    reverse=True,
                )
            )
            vote_counts[cluster_key] = dict(
                sorted(
                    (
                        (celltype, vote_bucket.get(celltype, len(values)))
                        for celltype, values in score_bucket.items()
                    ),
                    key=lambda item: (-item[1], item[0]),
                )
            )
        else:
            combined_scores[cluster_key] = {}
            vote_counts[cluster_key] = {}

    summary_scores = ul.unify_name(
        combined_scores,
        default_model,
        provider=default_provider,
        base_url=default_base_url,
    )
    remapped_vote_counts: dict[str, dict[str, int]] = {}
    for cluster_key, scores in summary_scores.items():
        original_votes = vote_counts.get(cluster_key, {})
        cluster_vote: dict[str, int] = {}
        for celltype in scores.keys():
            norm_key = _normalize_name(celltype)
            total_votes = sum(
                count
                for name, count in original_votes.items()
                if _normalize_name(name) == norm_key
            )
            if total_votes:
                cluster_vote[celltype] = total_votes
        remapped_vote_counts[cluster_key] = dict(
            sorted(cluster_vote.items(), key=lambda item: (-item[1], item[0]))
        )
    for cluster_key in vote_counts.keys():
        remapped_vote_counts.setdefault(cluster_key, vote_counts[cluster_key])
    vote_counts = remapped_vote_counts
    summary = {
        "combined_scores": summary_scores,
        "vote_counts": vote_counts,
        "total_models": len(normalized_configs),
    }
    if return_details:
        summary["raw_results"] = raw_results
    logger.info("Ensemble cell type annotation finished.")
    return summary


def get_subtype(
    input,
    out: Path | str | None = None,
    celltype: str = None,
    background: str = None,
    group: Iterable[str] | None = None,
    key: str = "rank_genes_groups",
    topnumber: int = 15,
    n_jobs: int | None = None,
    provider: str | None = None,
    model: str | None = None,
    base_url: str | None = None,
    rm_genes=True,
) -> dict:
    """\
    Annotating cell subtypes using LLM, providing cell types, supporting gene markers, reasons, and potential cell state annotations.
    
    Parameters
    ----------
    input : _type_
        An AnnData object or geneset dict
    out : Path | str | None, optional
        output path, by default None
    celltype : str, optional
        major cell type, by default None
    background : str, optional
        background information of input data, by default None
    group : Iterable[str] | None, optional
        which group, by default None
    key : str, optional
        deg group key, by default "rank_genes_groups"
    topnumber : int, optional
        select top gene for analysis, by default 15
    n_jobs : int | None, optional
        set multiple jobs for querying LLM, by default None
    provider : str | None, optional
        LLM provider, by default None
        "openai" for chatgpt
        "aliyun" for qwen
        "moonshot" for kimi
    model : str | None, optional
        set a model based on LLM provider, by default None
    base_url : str | None, optional
        customized LLM API url, by default None
    rm_genes : bool, optional
        remove rb and mt genes, by default True

    Returns
    -------
    dict
        a cell subtypes dict
    """
    logger.info("Starting cell subtype annotation workflow.")
    gene_dic = ul.get_gene_dict(input, group, key, topnumber, rm_genes)
    logger.info("Prepared %d gene sets for subtype analysis.", len(gene_dic))
    ot = ul.Outputor(out)
    ot.write("# Cell Subtype")

    genesets = []
    for k in gene_dic.keys():
        genestr = ",".join(gene_dic[k])
        genesetstr = f"geneset {k}: {genestr}"
        genesets.append(genesetstr)
    genesets_txt = "\n".join(genesets)
    msgs = [
        {
            "role": "user",
            "content": SUBTYPE_PROMPT.format(
                celltype=celltype, genesets=genesets_txt, background=background
            ),
        }
    ]

    logger.info(
        "Requesting subtype suggestions for '%s' from provider '%s' using model '%s'.",
        celltype or "unspecified cell type",
        provider or "default",
        model or "default",
    )
    response = query_model(
        msgs,
        provider=provider,
        model=model,
        base_url=base_url,
        sys_prompt=SYSTEM_PROMPT,
    )
    res_content = response.strip("```").strip("'''")
    ot.write(res_content)
    ot.close()
    subtype_lines = [
        line for line in res_content.split("\n") if line.startswith("###")
    ][-len(gene_dic) :]
    subtype_ls = [line.split(":")[1].strip() for line in subtype_lines]

    if len(gene_dic.keys()) == len(subtype_ls):
        subtype_dic = {k: subtype_ls[idx] for idx, k in enumerate(gene_dic.keys())}
    else:  # low-capability model may not get correct output
        subtype_dic = {}
        print(
            "The model may not be producing correct outputs; please try using a better model"
        )
    logger.info("Cell subtype annotation finished.")
    return subtype_dic


def check_celltype(
    input: AnnData | dict,
    out: Path | str = None,
    background: str = None,
    key: str = "rank_genes_groups",
    topnumber: int = 15,
    n_jobs: int | None = None,
    provider: str | None = None,
    model: str | None = None,
    group: str | Iterable[str] | None = None,
    base_url: str | None = None,
    rm_genes=True,
):
    """\
    Check the reason why genesets are annotated as these celltypes.

    Parameters
    ----------
    input : AnnData | dict
        An AnnData object or geneset dict
    out : Path | str, optional
        output path, by default None
    background : str, optional
        background information of input data, by default None
    key : str, optional
        deg group key, by default "rank_genes_groups"
    topnumber : int, optional
        select top gene for analysis, by default 15
    n_jobs : int | None, optional
        set multiple jobs for querying LLM, by default None
    provider : str | None, optional
        LLM provider, by default None
        "openai" for chatgpt
        "aliyun" for qwen
        "moonshot" for kimi
    model : str | None, optional
        set a model based on LLM provider, by default None
    group : str | Iterable[str] | None, optional
        _description_, by default None
    base_url : str | None, optional
        customized LLM API url, by default None
    rm_genes : bool, optional
        remove rb and mt genes, by default True

    Returns
    -------
    None
    """
    sys_prompt = SYSTEM_PROMPT
    logger.info("Starting cell type checking workflow.")
    gene_dic = ul.get_gene_dict(input, group, key, topnumber, rm_genes)
    logger.info("Prepared %d gene sets for checking.", len(gene_dic))

    ot = ul.Outputor(out)
    ot.write("# CellType checking")

    genesets = []
    for k in gene_dic.keys():
        genestr = ",".join(gene_dic[k])
        genesetstr = f"Celltype: {k}; geneset: {genestr}"
        genesets.append(genesetstr)
    genesets_txt = "\n".join(genesets)
    msgs = [
        {
            "role": "user",
            "content": CHECK_TYPE_PROMPT.format(
                genesets=genesets_txt, background=background
            ),
        }
    ]

    logger.info(
        "Submitting cell type review request to provider '%s' using model '%s'.",
        provider or "default",
        model or "default",
    )
    response = query_model(
        msgs, provider=provider, model=model, base_url=base_url, sys_prompt=sys_prompt
    )
    res_content = response.strip("```").strip("'''")
    ot.write(res_content)
    ot.close()
    logger.info("Cell type checking completed.")


def get_cellstate(
    input: AnnData | dict,
    out: Path | str = None,
    background: str = None,
    pathway: dict | None = None,
    deg_key: str = "rank_genes_groups",
    topnumber: int = 15,
    n_jobs: int | None = None,
    provider: str | None = None,
    model: str | None = None,
    group: str | Iterable[str] | None = None,
    base_url: str | None = None,
    rm_genes=True,
) -> dict:
    """\
    Annotating cell type state using LLM.

    Parameters
    ----------
    input : AnnData | dict
        An AnnData object or geneset dict
    out : Path | str, optional
        output path, by default None
    background : str, optional
        background information of input data, by default None
    deg_key : str, optional
        deg key, by default "rank_genes_groups"
    topnumber : int, optional
        select top gene for analysis, by default 15
    n_jobs : int | None, optional
        set multiple jobs for querying LLM, by default None
    provider : str| None, optional
        LLM provider, by default None
        "openai" for chatgpt
        "aliyun" for qwen
        "deepseek" for DeepSeek
        "anthropic" for claude
    model : str | None, optional
        set a model based on LLM provider, by default None
    group : str | Iterable, optional
         Which group, by default None
    base_url : str | None, optional
        customized LLM API url, by default None
    rm_genes : bool, optional
        remove rb and mt genes, by default True

    Returns
    -------
    None
    """
    sys_prompt = SYSTEM_PROMPT
    logger.info("Starting cell state analysis workflow.")
    gene_dic = ul.get_gene_dict(input, group, deg_key, topnumber, rm_genes)
    logger.info("Prepared %d gene sets for state analysis.", len(gene_dic))
    ot = ul.Outputor(out)
    ot.write("# Cell State Analysis")
    ot.write(
        "GPTBioInsightor is powered by AI, so mistakes are possible. Review output carefully before use"
    )

    with ThreadPoolExecutor(max_workers=n_jobs) as executor:
        futures = {}
        pathway = {} if pathway is None else pathway
        pathway_txt_dic = {}
        for k, genes in gene_dic.items():
            agent = Agent(
                model=model, provider=provider, sys_prompt=sys_prompt, base_url=base_url
            )
            gene_txt = f"   - cluster {k}: {','.join(genes[:topnumber])}"
            cluster_pathway = pathway.get(k, {})
            pw_txt = ""
            for db, pw in cluster_pathway.items():
                pw_txt += f"    - {db}: {','.join(pw)}\n"
            pathway_txt_dic[k] = pw_txt
            pct_txt = CELLSTATE_PROMPT.format(
                celltype=k, gene=gene_txt, background=background, pathway=pw_txt
            )
            logger.info("Generating cell state narrative for cluster %s.", k)
            future = executor.submit(agent.query, pct_txt)
            futures[k] = future
        # score_dic = {}
        for k, future in futures.items():
            reps = future.result()
            ot.write(f"## cluster geneset {k}\n")
            ot.write("### Gene List\n")
            ot.write(f"Top genes\n```\n{','.join(gene_dic[k])}\n```\n")
            ot.write(f"enrichment pathway\n```\n{pathway_txt_dic[k]}```\n")
            ot.write("### cell state thinking\n")
            ot.write(reps)
            logger.info("Finished cell state report for cluster %s.", k)
    logger.info("Cell state analysis completed.")
