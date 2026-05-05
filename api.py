# src/api.py
# ─────────────────────────────────────────────────────────────
# FastAPI server exposing the OncoPredict ML pipeline.
# All response schemas exactly match what the frontend expects.
#
# Endpoints:
#   GET  /                        Health check
#   GET  /drugs                   List all available drugs
#   POST /predict                 Run prediction (CSV upload + drug)
#   GET  /predict/demo/{drug_key} Demo prediction with synthetic data
#
# Run locally:
#   uvicorn src.api:app --reload --port 8000
# ─────────────────────────────────────────────────────────────

import logging
from typing import Optional, Any
from contextlib import asynccontextmanager

import numpy as np
import pandas as pd
from fastapi import FastAPI, File, Form, UploadFile, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

from src.config import DRUGS, KEY_GENES, GENE_PATHWAY_MAP
from src.predictor import Predictor, PredictionResult

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger(__name__)

# ── App instance ───────────────────────────────────────────
predictor = Predictor(use_fallback_network=True)


@asynccontextmanager
async def lifespan(app: FastAPI):
    # Pre-load model at startup so first request is fast
    try:
        predictor._ensure_model_loaded()
        logger.info("Model pre-loaded successfully.")
    except FileNotFoundError:
        logger.warning(
            "No trained model found. Run 'python -m src.model_trainer' first."
        )
    yield


app = FastAPI(
    title="OncoPredict API",
    description="AI-Based Breast Cancer Drug Response Predictor",
    version="2.1.0",
    lifespan=lifespan,
)

# ── CORS (allows frontend HTML file to call API) ───────────
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],       # restrict to your domain in production
    allow_methods=["GET", "POST", "OPTIONS"],
    allow_headers=["*"],
)


# ── Response Models ────────────────────────────────────────

class GeneChangeItem(BaseModel):
    gene: str
    value: float


class DrugInfo(BaseModel):
    key: str
    name: str
    display_name: str
    target_genes: list[str]
    target_label: str
    description: str


class PredictResponse(BaseModel):
    classification: str
    confidence: float
    gene_changes: list[GeneChangeItem]
    sensitivity_pct: float
    bar_chart_data: dict[str, Any]
    network_data: dict[str, Any]
    drug_key: str
    drug_name: str
    drug_target: str
    is_sensitive: bool
    effectiveness_score: float
    cascade_iterations: int
    n_genes_input: int
    n_genes_matched: int
    gene_changes_detailed: list[dict[str, Any]]
    error: Optional[str] = None


def _result_to_payload(result: PredictionResult) -> PredictResponse:
    detailed = result.gene_changes or []
    labels = [g.get("gene", "") for g in detailed]
    pre_vals = [float(g.get("pre", 0.0)) for g in detailed]
    post_vals = [float(g.get("post", 0.0)) for g in detailed]
    delta_vals = [float(g.get("delta", 0.0)) for g in detailed]
    simple_gene_changes = [
        GeneChangeItem(gene=g.get("gene", ""), value=float(g.get("delta", 0.0)))
        for g in detailed
    ]

    target_gene = (result.drug_target or "").replace("Target Gene:", "").split(",")[0].strip()
    nodes = [{"id": "drug"}, {"id": target_gene}] if target_gene else [{"id": "drug"}]
    edges = [{"source": "drug", "target": target_gene}] if target_gene else []
    for g in detailed[:8]:
        gene = g.get("gene", "")
        if gene and all(n["id"] != gene for n in nodes):
            nodes.append({"id": gene})
        if target_gene and gene and gene != target_gene:
            edges.append({"source": target_gene, "target": gene})

    return PredictResponse(
        classification=result.classification.lower(),
        confidence=float(result.confidence),
        gene_changes=simple_gene_changes,
        sensitivity_pct=float(result.sensitivity_pct),
        bar_chart_data={
            "labels": labels,
            "values": delta_vals,
            "pre": pre_vals,
            "post": post_vals,
        },
        network_data={"nodes": nodes, "edges": edges},
        drug_key=result.drug_key,
        drug_name=result.drug_name,
        drug_target=result.drug_target,
        is_sensitive=result.is_sensitive,
        effectiveness_score=float(result.effectiveness_score),
        cascade_iterations=int(result.cascade_iterations),
        n_genes_input=int(result.n_genes_input),
        n_genes_matched=int(result.n_genes_matched),
        gene_changes_detailed=detailed,
    )


# ── Endpoints ──────────────────────────────────────────────

@app.get("/", tags=["Health"])
async def health_check():
    return {"status": "ok", "version": "2.1.0", "service": "OncoPredict API"}


@app.get("/drugs", response_model=list[DrugInfo], tags=["Drugs"])
async def list_drugs():
    """Return all available drugs with metadata for the frontend drug selector."""
    return [
        DrugInfo(
            key=key,
            name=cfg["name"],
            display_name=cfg["display_name"],
            target_genes=cfg["target_genes"],
            target_label=cfg["target_label"],
            description=cfg["description"],
        )
        for key, cfg in DRUGS.items()
    ]


@app.post("/predict", response_model=PredictResponse, tags=["Prediction"])
async def predict(
    file: UploadFile = File(..., description="Gene expression CSV (gene, log2_value)"),
    drug_key: str = Form(..., description="Drug key e.g. 'trastuzumab'"),
    cascade_iterations: int = Form(5, ge=1, le=20),
    normalization: str = Form("log2"),
):
    """
    Main prediction endpoint.

    Accepts a CSV gene expression file and a drug key.
    Returns classification (sensitive/resistant), confidence score,
    and per-gene expression changes before and after drug simulation.

    Example curl:
        curl -X POST http://localhost:8000/predict \\
             -F "file=@patient_expr.csv" \\
             -F "drug_key=trastuzumab" \\
             -F "cascade_iterations=5"
    """
    if drug_key not in DRUGS:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown drug_key '{drug_key}'. Valid: {list(DRUGS.keys())}",
        )

    if not file.filename.endswith(".csv"):
        raise HTTPException(status_code=400, detail="File must be a .csv")

    csv_bytes = await file.read()
    if len(csv_bytes) == 0:
        raise HTTPException(status_code=400, detail="Uploaded file is empty.")
    if len(csv_bytes) > 10 * 1024 * 1024:  # 10 MB limit
        raise HTTPException(status_code=413, detail="File too large (max 10 MB).")

    result: PredictionResult = predictor.predict(
        csv_bytes=csv_bytes,
        drug_key=drug_key,
        cascade_iterations=cascade_iterations,
        normalization=normalization,
    )

    if result.error:
        raise HTTPException(status_code=422, detail=result.error)

    return _result_to_payload(result)


@app.get("/predict/demo/{drug_key}", response_model=PredictResponse, tags=["Prediction"])
async def predict_demo(
    drug_key: str,
    cascade_iterations: int = 5,
):
    """
    Demo endpoint — runs prediction on a synthetic HER2+ patient profile.
    No file upload required. Useful for frontend testing without real data.
    """
    if drug_key not in DRUGS:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown drug_key '{drug_key}'. Valid: {list(DRUGS.keys())}",
        )

    # Build a synthetic demo CSV in memory
    rng = np.random.RandomState(seed=hash(drug_key) % (2**31))
    base_vals = rng.uniform(4.0, 9.0, len(KEY_GENES))

    # Amplify target genes for realism
    drug_conf = DRUGS[drug_key]
    for tg in drug_conf["target_genes"]:
        if tg in KEY_GENES:
            base_vals[KEY_GENES.index(tg)] = rng.uniform(8.0, 11.5)

    csv_lines = "\n".join(
        f"{gene},{val:.4f}" for gene, val in zip(KEY_GENES, base_vals)
    )
    csv_bytes = csv_lines.encode("utf-8")

    result = predictor.predict(
        csv_bytes=csv_bytes,
        drug_key=drug_key,
        cascade_iterations=cascade_iterations,
        normalization="none",  # already normalized
    )

    if result.error:
        raise HTTPException(status_code=422, detail=result.error)

    return _result_to_payload(result)
