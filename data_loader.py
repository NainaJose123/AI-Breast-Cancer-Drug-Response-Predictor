# src/data_loader.py
# ─────────────────────────────────────────────────────────────
# Downloads and preprocesses breast cancer gene expression data
# from NCBI GEO dataset GSE81538 (405 tumor samples, ~20K genes).
#
# Usage:
#   from src.data_loader import load_geo_dataset, preprocess_expression
#   expr_df = load_geo_dataset()
#   expr_clean = preprocess_expression(expr_df)
# ─────────────────────────────────────────────────────────────

import os
import io
import logging
import numpy as np
import pandas as pd
import GEOparse

logger = logging.getLogger(__name__)

GEO_ACCESSION = "GSE81538"
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
CACHED_EXPR_PATH = os.path.join(DATA_DIR, "gse81538_expression.csv")
CACHED_META_PATH = os.path.join(DATA_DIR, "gse81538_metadata.csv")


def load_geo_dataset(force_download: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Download GSE81538 from GEO and return:
      - expression DataFrame:  genes × samples  (log2 values)
      - metadata DataFrame:    samples × attributes
    Results are cached to disk after first download.
    """
    os.makedirs(DATA_DIR, exist_ok=True)

    if not force_download and os.path.exists(CACHED_EXPR_PATH):
        logger.info("Loading cached expression data from %s", CACHED_EXPR_PATH)
        expr = pd.read_csv(CACHED_EXPR_PATH, index_col=0)
        meta = pd.read_csv(CACHED_META_PATH, index_col=0)
        return expr, meta

    logger.info("Downloading %s from NCBI GEO ...", GEO_ACCESSION)
    gse = GEOparse.get_GEO(geo=GEO_ACCESSION, destdir=DATA_DIR, silent=True)

    # ── Expression matrix ──────────────────────────────────
    pivot_df = gse.pivot_samples("VALUE")   # genes × samples
    expr = pivot_df.copy()

    # ── Sample metadata ────────────────────────────────────
    meta_rows = []
    for gsm_name, gsm in gse.gsms.items():
        row = {"sample_id": gsm_name}
        row.update({k: v[0] if isinstance(v, list) else v
                    for k, v in gsm.metadata.items()})
        meta_rows.append(row)
    meta = pd.DataFrame(meta_rows).set_index("sample_id")

    # Keep only tumor samples (filter out "normal")
    if "source_name_ch1" in meta.columns:
        tumor_mask = ~meta["source_name_ch1"].str.lower().str.contains("normal", na=False)
        tumor_samples = meta.index[tumor_mask].tolist()
        expr = expr[[c for c in tumor_samples if c in expr.columns]]
        meta = meta.loc[tumor_samples]
        logger.info("Retained %d tumor samples after filtering", len(tumor_samples))

    expr.to_csv(CACHED_EXPR_PATH)
    meta.to_csv(CACHED_META_PATH)
    logger.info("Saved expression data → %s", CACHED_EXPR_PATH)
    return expr, meta


def preprocess_expression(
    expr: pd.DataFrame,
    method: str = "log2",
    min_threshold: float = 1.0,
) -> pd.DataFrame:
    """
    Clean and normalize a gene expression DataFrame.

    Args:
        expr:          genes × samples DataFrame (raw counts or log values)
        method:        'log2' | 'zscore' | 'minmax' | 'none'
        min_threshold: drop genes where ALL samples are below this value

    Returns:
        Cleaned, normalized DataFrame with same shape convention.
    """
    # 1. Convert to numeric, drop rows with all-NaN
    expr = expr.apply(pd.to_numeric, errors="coerce")
    expr = expr.dropna(how="all")

    # 2. Drop low-expression genes
    expr = expr[expr.max(axis=1) >= min_threshold]

    # 3. Fill remaining NaN with per-gene median
    expr = expr.T.fillna(expr.median(axis=1)).T

    # 4. Normalize
    if method == "log2":
        # Shift to ensure all values > 0 before log
        shift = max(0, -expr.values.min()) + 1
        expr = np.log2(expr + shift)
    elif method == "zscore":
        from scipy.stats import zscore as sp_zscore
        expr = expr.apply(lambda row: sp_zscore(row), axis=1, result_type="expand")
        expr.columns = expr.columns  # preserve column names
    elif method == "minmax":
        mn = expr.min(axis=1)
        mx = expr.max(axis=1)
        expr = expr.subtract(mn, axis=0).divide((mx - mn).replace(0, 1), axis=0)

    logger.info("Preprocessed expression matrix: %d genes × %d samples", *expr.shape)
    return expr


def load_patient_csv(file_content: bytes | str, method: str = "log2") -> pd.Series:
    """
    Parse an uploaded patient CSV file.

    Expected format (no header):
        ERBB2,8.4231
        TP53,6.1902
        ...

    Returns:
        pd.Series indexed by gene symbol, values = expression (normalized).
    """
    if isinstance(file_content, bytes):
        file_content = file_content.decode("utf-8")

    df = pd.read_csv(io.StringIO(file_content), header=None, names=["gene", "value"])
    df = df.dropna()
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df = df.dropna(subset=["value"])
    df = df.drop_duplicates(subset="gene")
    series = df.set_index("gene")["value"]

    # Normalize if not already log2
    if method == "log2" and series.max() > 30:
        shift = max(0, -series.min()) + 1
        series = np.log2(series + shift)

    return series
