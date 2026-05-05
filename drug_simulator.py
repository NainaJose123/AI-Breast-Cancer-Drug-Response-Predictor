# src/drug_simulator.py
# ─────────────────────────────────────────────────────────────
# Simulates drug perturbation on a patient PPI network.
#
# Step 1 — Apply inhibition:
#   new_expr = orig_expr × (1 − inhibition_factor)
#
# Step 2 — Cascade propagation (iterative):
#   gene_new = gene_old + damping × Σ(neighbor_expr × edge_weight)
#
# Usage:
#   from src.drug_simulator import simulate_drug
#   pre_state, post_state = simulate_drug(G, drug_config, iterations=5)
# ─────────────────────────────────────────────────────────────

import copy
import logging
import numpy as np
import networkx as nx
import pandas as pd

from src.config import CASCADE_ITERATIONS_DEFAULT, CASCADE_DAMPING, GENE_PATHWAY_MAP

logger = logging.getLogger(__name__)


def simulate_drug(
    G: nx.Graph,
    drug_config: dict,
    iterations: int = CASCADE_ITERATIONS_DEFAULT,
    damping: float = CASCADE_DAMPING,
) -> tuple[dict, dict]:
    """
    Simulate drug action on a patient PPI network.

    Args:
        G:           networkx Graph with node attr 'expression'
        drug_config: dict from config.DRUGS (contains target_genes, inhibition_factor)
        iterations:  number of cascade propagation steps
        damping:     step size for signal propagation (0 < damping <= 1)

    Returns:
        pre_state:  {gene: expression_value} before drug
        post_state: {gene: expression_value} after drug + cascade
    """
    # ── Capture pre-drug state ─────────────────────────────
    pre_state = {
        node: G.nodes[node].get("expression", 0.0)
        for node in G.nodes()
    }

    # ── Step 1: Apply direct inhibition to target genes ────
    G_post = copy.deepcopy(G)
    inhibition = drug_config["inhibition_factor"]

    for target in drug_config["target_genes"]:
        if target in G_post.nodes:
            orig = G_post.nodes[target]["expression"]
            G_post.nodes[target]["expression"] = orig * (1.0 - inhibition)
            logger.debug(
                "Inhibited %s: %.3f → %.3f (factor=%.2f)",
                target, orig, G_post.nodes[target]["expression"], inhibition,
            )

    # ── Step 2: Iterative cascade propagation ─────────────
    for iteration in range(iterations):
        new_expressions = {}
        for node in G_post.nodes():
            current = G_post.nodes[node]["expression"]
            neighbors = list(G_post.neighbors(node))
            if not neighbors:
                new_expressions[node] = current
                continue
            # Weighted sum of neighbor expressions
            neighbor_signal = sum(
                G_post.nodes[nb]["expression"] * G_post[node][nb].get("weight", 1.0)
                for nb in neighbors
            )
            new_val = current + damping * neighbor_signal
            # Clamp to biologically reasonable range [0, 15]
            new_expressions[node] = float(np.clip(new_val, 0.0, 15.0))

        for node, val in new_expressions.items():
            G_post.nodes[node]["expression"] = val

        logger.debug("Cascade iteration %d complete.", iteration + 1)

    # ── Capture post-drug state ────────────────────────────
    post_state = {
        node: G_post.nodes[node].get("expression", 0.0)
        for node in G_post.nodes()
    }

    return pre_state, post_state


def compute_gene_changes(
    pre_state: dict,
    post_state: dict,
    key_genes: list[str],
) -> list[dict]:
    """
    Compute per-gene expression changes and return structured list
    matching the frontend's expected format:

    [
        {
            "gene": "ERBB2",
            "pathway": "RTK Signaling",
            "pre": 8.42,
            "post": 2.18,
            "delta": -6.24,
            "delta_pct": -74.1
        },
        ...
    ]
    """
    results = []
    for gene in key_genes:
        if gene not in pre_state:
            continue
        pre_val = round(pre_state[gene], 4)
        post_val = round(post_state.get(gene, pre_val), 4)
        delta = round(post_val - pre_val, 4)
        delta_pct = round((delta / pre_val * 100) if pre_val != 0 else 0.0, 2)
        results.append({
            "gene": gene,
            "pathway": GENE_PATHWAY_MAP.get(gene, "Other"),
            "pre": pre_val,
            "post": post_val,
            "delta": delta,
            "delta_pct": delta_pct,
        })

    # Sort by absolute delta descending (most affected genes first)
    results.sort(key=lambda x: abs(x["delta"]), reverse=True)
    return results
