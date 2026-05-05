# src/network_builder.py
# ─────────────────────────────────────────────────────────────
# Builds a patient-specific protein–protein interaction (PPI)
# network from STRING database.
#
# Nodes  = genes present in patient expression profile
# Edges  = STRING interactions with confidence >= threshold
# Weight = STRING combined score (0–1)
#
# Usage:
#   from src.network_builder import build_ppi_network, attach_expression
#   G = build_ppi_network(gene_list)
#   G = attach_expression(G, patient_series)
# ─────────────────────────────────────────────────────────────

import os
import logging
import requests
import pandas as pd
import networkx as nx

from src.config import STRING_CONFIDENCE_MIN

logger = logging.getLogger(__name__)

STRING_API = "https://string-db.org/api/tsv/network"
STRING_SPECIES = 9606        # Homo sapiens
NETWORK_CACHE_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "networks")


def _cache_path(gene_set_hash: str) -> str:
    os.makedirs(NETWORK_CACHE_DIR, exist_ok=True)
    return os.path.join(NETWORK_CACHE_DIR, f"ppi_{gene_set_hash}.csv")


def fetch_string_interactions(
    genes: list[str],
    confidence_min: float = STRING_CONFIDENCE_MIN,
) -> pd.DataFrame:
    """
    Query the STRING REST API for interactions among a list of genes.
    Returns a DataFrame with columns: gene_a, gene_b, score.
    Results are cached to avoid redundant API calls.
    """
    gene_hash = str(abs(hash(frozenset(genes))))[:12]
    cache = _cache_path(gene_hash)

    if os.path.exists(cache):
        logger.info("Loading cached STRING interactions from %s", cache)
        return pd.read_csv(cache)

    logger.info("Querying STRING API for %d genes ...", len(genes))
    params = {
        "identifiers": "%0d".join(genes),
        "species": STRING_SPECIES,
        "required_score": int(confidence_min * 1000),
        "caller_identity": "oncopredict_research",
    }
    try:
        resp = requests.post(STRING_API, data=params, timeout=30)
        resp.raise_for_status()
        lines = [l.split("\t") for l in resp.text.strip().split("\n")]
        if len(lines) < 2:
            logger.warning("STRING returned no interactions for given gene list.")
            return pd.DataFrame(columns=["gene_a", "gene_b", "score"])
        header = lines[0]
        rows = lines[1:]
        df = pd.DataFrame(rows, columns=header)
        # STRING returns preferredName_A / preferredName_B / score columns
        df = df.rename(columns={
            "preferredName_A": "gene_a",
            "preferredName_B": "gene_b",
            "score": "score",
        })
        df["score"] = pd.to_numeric(df["score"], errors="coerce") / 1000.0
        df = df[["gene_a", "gene_b", "score"]].dropna()
        df = df[df["score"] >= confidence_min]
        df.to_csv(cache, index=False)
        logger.info("Fetched %d interactions, saved to cache.", len(df))
        return df
    except requests.RequestException as e:
        logger.error("STRING API request failed: %s", e)
        return pd.DataFrame(columns=["gene_a", "gene_b", "score"])


def build_ppi_network(
    genes: list[str],
    confidence_min: float = STRING_CONFIDENCE_MIN,
) -> nx.Graph:
    """
    Build an undirected weighted PPI network for a list of genes.

    Returns:
        nx.Graph with nodes=genes and edge attr 'weight'=confidence score.
    """
    interactions = fetch_string_interactions(genes, confidence_min)
    G = nx.Graph()
    G.add_nodes_from(genes)

    for _, row in interactions.iterrows():
        a, b, w = row["gene_a"], row["gene_b"], float(row["score"])
        if a in G and b in G:
            G.add_edge(a, b, weight=w)

    logger.info(
        "PPI network: %d nodes, %d edges (confidence >= %.2f)",
        G.number_of_nodes(), G.number_of_edges(), confidence_min,
    )
    return G


def attach_expression(G: nx.Graph, expression: pd.Series) -> nx.Graph:
    """
    Attach patient-specific expression values to each node of the graph.
    Nodes without expression data receive the median value.
    """
    median_val = float(expression.median())
    for node in G.nodes():
        G.nodes[node]["expression"] = float(expression.get(node, median_val))
    return G


def build_fallback_network(genes: list[str], expression: pd.Series) -> nx.Graph:
    """
    Lightweight fallback network when STRING API is unavailable.
    Connects biologically related genes based on known pathway groupings
    hardcoded from KEGG breast cancer pathway (hsa05224).
    """
    KNOWN_EDGES = [
        ("ERBB2", "AKT1", 0.9), ("ERBB2", "MAPK3", 0.85), ("ERBB2", "PIK3CA", 0.88),
        ("EGFR", "AKT1", 0.82), ("EGFR", "MAPK1", 0.83), ("EGFR", "ERBB2", 0.91),
        ("PIK3CA", "AKT1", 0.95), ("AKT1", "MTOR", 0.90), ("AKT1", "GSK3B", 0.85),
        ("AKT1", "FOXO1", 0.82), ("PTEN", "PIK3CA", 0.88), ("PTEN", "AKT1", 0.86),
        ("MAPK1", "MAPK3", 0.93), ("MAPK3", "MYC", 0.78), ("MYC", "CCND1", 0.80),
        ("CCND1", "CDK4", 0.92), ("CCND1", "CDK6", 0.90), ("CDK4", "RB1", 0.87),
        ("CDK6", "RB1", 0.85), ("RB1", "E2F1", 0.88), ("E2F1", "PCNA", 0.80),
        ("TP53", "CDKN2A", 0.85), ("TP53", "BCL2", 0.78), ("TP53", "MYC", 0.72),
        ("ESR1", "GATA3", 0.82), ("ESR1", "PGR", 0.90), ("ESR1", "MYC", 0.75),
        ("BRCA1", "TP53", 0.80), ("STAT3", "MYC", 0.76), ("VEGFA", "MTOR", 0.71),
    ]
    G = nx.Graph()
    G.add_nodes_from(genes)
    for a, b, w in KNOWN_EDGES:
        if a in G and b in G:
            G.add_edge(a, b, weight=w)
    return attach_expression(G, expression)
