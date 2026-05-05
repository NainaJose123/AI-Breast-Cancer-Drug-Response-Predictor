# src/predictor.py
# ─────────────────────────────────────────────────────────────
# Core inference engine — orchestrates all pipeline stages:
#
#   1. Load patient expression (from CSV or pre-loaded Series)
#   2. Build PPI network (STRING or fallback)
#   3. Simulate drug cascade
#   4. Build ML feature vector
#   5. Run XGBoost classifier → class + confidence
#   6. Return structured PredictionResult
#
# Usage:
#   from src.predictor import Predictor
#   predictor = Predictor()
#   result = predictor.predict(csv_bytes, drug_key="trastuzumab")
# ─────────────────────────────────────────────────────────────

import logging
import numpy as np
import pandas as pd
from dataclasses import dataclass, field

from src.config import DRUGS, KEY_GENES, CASCADE_ITERATIONS_DEFAULT
from src.data_loader import load_patient_csv
from src.network_builder import build_ppi_network, attach_expression, build_fallback_network
from src.drug_simulator import simulate_drug, compute_gene_changes
from src.model_trainer import build_feature_vector, load_model

logger = logging.getLogger(__name__)


@dataclass
class PredictionResult:
    """
    Structured prediction result — mirrors the frontend's data schema.
    """
    drug_key: str
    drug_name: str
    drug_target: str
    is_sensitive: bool
    classification: str          # "Sensitive" | "Resistant"
    confidence: float            # 0.0 – 1.0
    effectiveness_score: float   # 0.0 – 2.0 (log-odds-style score)
    sensitivity_pct: float       # 0 – 100

    gene_changes: list[dict] = field(default_factory=list)
    # Each entry: {gene, pathway, pre, post, delta, delta_pct}

    cascade_iterations: int = CASCADE_ITERATIONS_DEFAULT
    n_genes_input: int = 0
    n_genes_matched: int = 0
    error: str | None = None


class Predictor:
    """
    Loads the trained model once and exposes a predict() method.
    Thread-safe for use inside FastAPI.
    """

    def __init__(self, use_fallback_network: bool = True):
        """
        Args:
            use_fallback_network: If True, use hardcoded PPI network
                instead of querying STRING (faster, no internet needed).
                Set to False in production when STRING API is available.
        """
        self.use_fallback = use_fallback_network
        self._model = None
        self._drug_keys = None

    def _ensure_model_loaded(self):
        if self._model is None:
            logger.info("Loading ML model from disk ...")
            self._model, self._drug_keys = load_model()

    def predict(
        self,
        csv_bytes: bytes,
        drug_key: str,
        cascade_iterations: int = CASCADE_ITERATIONS_DEFAULT,
        normalization: str = "log2",
    ) -> PredictionResult:
        """
        Full pipeline: CSV → network → simulation → ML → result.

        Args:
            csv_bytes:          Raw bytes of uploaded CSV file
            drug_key:           One of the keys in config.DRUGS
            cascade_iterations: Number of signal propagation steps
            normalization:      'log2' | 'zscore' | 'minmax' | 'none'

        Returns:
            PredictionResult dataclass
        """
        if drug_key not in DRUGS:
            return PredictionResult(
                drug_key=drug_key,
                drug_name=drug_key,
                drug_target="",
                is_sensitive=False,
                classification="Unknown",
                confidence=0.0,
                effectiveness_score=0.0,
                sensitivity_pct=0.0,
                error=f"Unknown drug key: {drug_key}",
            )

        drug_conf = DRUGS[drug_key]

        # ── 1. Parse patient CSV ──────────────────────────
        try:
            patient_series = load_patient_csv(csv_bytes, method=normalization)
        except Exception as e:
            logger.error("Failed to parse CSV: %s", e)
            return self._error_result(drug_key, drug_conf, str(e))

        n_genes_input = len(patient_series)

        # ── 2. Match genes to KEY_GENES ───────────────────
        # Use patient data where available, fill with median otherwise
        median_val = float(patient_series.median())
        full_expression = pd.Series(
            {g: patient_series.get(g, median_val) for g in KEY_GENES}
        )
        n_genes_matched = int(patient_series.index.isin(KEY_GENES).sum())

        # ── 3. Build PPI network ──────────────────────────
        try:
            if self.use_fallback:
                G = build_fallback_network(KEY_GENES, full_expression)
            else:
                G = build_ppi_network(KEY_GENES)
                G = attach_expression(G, full_expression)
        except Exception as e:
            logger.warning("Network build failed, using fallback: %s", e)
            G = build_fallback_network(KEY_GENES, full_expression)

        # ── 4. Simulate drug cascade ──────────────────────
        pre_state, post_state = simulate_drug(
            G, drug_conf, iterations=cascade_iterations
        )

        # ── 5. Build feature vector ───────────────────────
        self._ensure_model_loaded()
        feature_vec = build_feature_vector(
            pre_state, post_state, drug_key, KEY_GENES, self._drug_keys
        )

        # ── 6. ML Inference ───────────────────────────────
        X = feature_vec.reshape(1, -1)
        proba = self._model.predict_proba(X)[0]  # [P(resistant), P(sensitive)]
        pred_class = int(self._model.predict(X)[0])
        confidence = float(proba[pred_class])
        sensitive_proba = float(proba[1])

        is_sensitive = pred_class == 1
        classification = "Sensitive" if is_sensitive else "Resistant"

        # Effectiveness score: log-odds style, range approx 0–2.0
        # Maps confidence 0.5→1.0 to score 0.0→2.0
        effectiveness_score = round(
            np.log1p(sensitive_proba * 2.718) / np.log1p(2.718) * 2.0, 2
        )
        sensitivity_pct = round(sensitive_proba * 100, 1)

        # ── 7. Gene change table ──────────────────────────
        gene_changes = compute_gene_changes(pre_state, post_state, KEY_GENES)

        logger.info(
            "Prediction complete: drug=%s class=%s confidence=%.3f",
            drug_key, classification, confidence,
        )

        return PredictionResult(
            drug_key=drug_key,
            drug_name=drug_conf["name"],
            drug_target=drug_conf["target_label"],
            is_sensitive=is_sensitive,
            classification=classification,
            confidence=round(confidence, 4),
            effectiveness_score=effectiveness_score,
            sensitivity_pct=sensitivity_pct,
            gene_changes=gene_changes,
            cascade_iterations=cascade_iterations,
            n_genes_input=n_genes_input,
            n_genes_matched=n_genes_matched,
        )

    @staticmethod
    def _error_result(drug_key: str, drug_conf: dict, msg: str) -> PredictionResult:
        return PredictionResult(
            drug_key=drug_key,
            drug_name=drug_conf["name"],
            drug_target=drug_conf.get("target_label", ""),
            is_sensitive=False,
            classification="Error",
            confidence=0.0,
            effectiveness_score=0.0,
            sensitivity_pct=0.0,
            error=msg,
        )
