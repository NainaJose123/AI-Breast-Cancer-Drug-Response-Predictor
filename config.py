# src/config.py
# ─────────────────────────────────────────────────────────────
# Central configuration: drug definitions, target genes, IC50
# inhibition factors sourced from GDSC v2 pharmacogenomic data.
# ─────────────────────────────────────────────────────────────

DRUGS = {
    "trastuzumab": {
        "name": "Trastuzumab",
        "display_name": "Trastuzumab (Herceptin)",
        "target_genes": ["ERBB2"],
        "target_label": "Target Gene: ERBB2",
        "description": (
            "Monoclonal antibody targeting HER2/ERBB2 receptor. "
            "Indicated for HER2-positive breast cancer. "
            "IC\u2085\u2080 data sourced from GDSC v2 pharmacogenomic database."
        ),
        # Inhibition factor derived from median IC50 across GDSC cell lines
        "inhibition_factor": 0.74,
        "gdsc_drug_id": 1378,
    },
    "gefitinib": {
        "name": "Gefitinib",
        "display_name": "Gefitinib (Iressa)",
        "target_genes": ["EGFR"],
        "target_label": "Target Gene: EGFR",
        "description": (
            "EGFR tyrosine kinase inhibitor. Reduces EGFR-driven "
            "proliferation signaling in EGFR-amplified tumors. IC\u2085\u2080 from GDSC."
        ),
        "inhibition_factor": 0.52,
        "gdsc_drug_id": 1010,
    },
    "lapatinib": {
        "name": "Lapatinib",
        "display_name": "Lapatinib (Tykerb)",
        "target_genes": ["ERBB2", "EGFR"],
        "target_label": "Target Gene: ERBB2, EGFR",
        "description": (
            "Dual kinase inhibitor targeting HER2 and EGFR. "
            "Used in HER2-positive breast cancer with capecitabine."
        ),
        "inhibition_factor": 0.64,
        "gdsc_drug_id": 119,
    },
    "pi3k": {
        "name": "PI3K Inhibitor (BYL719)",
        "display_name": "BYL719 (Alpelisib)",
        "target_genes": ["PIK3CA"],
        "target_label": "Target Gene: PIK3CA",
        "description": (
            "PI3K\u03b1-specific inhibitor. Effective in PIK3CA-mutated "
            "HR-positive HER2-negative breast cancer with fulvestrant."
        ),
        "inhibition_factor": 0.69,
        "gdsc_drug_id": 1560,
    },
    "tamoxifen": {
        "name": "Tamoxifen",
        "display_name": "Tamoxifen",
        "target_genes": ["ESR1"],
        "target_label": "Target Gene: ESR1 (\u03b1)",
        "description": (
            "Selective estrogen receptor modulator. Competitively inhibits "
            "estrogen binding to ER\u03b1 in hormone receptor-positive tumors."
        ),
        "inhibition_factor": 0.51,
        "gdsc_drug_id": 1199,
    },
    "palbociclib": {
        "name": "Palbociclib",
        "display_name": "Palbociclib (Ibrance)",
        "target_genes": ["CDK4", "CDK6"],
        "target_label": "Target Gene: CDK4, CDK6",
        "description": (
            "CDK4/6 inhibitor blocking cell cycle progression at G1/S. "
            "Used in HR-positive HER2-negative metastatic breast cancer."
        ),
        "inhibition_factor": 0.71,
        "gdsc_drug_id": 1054,
    },
}

# Key cancer-related genes tracked in output
KEY_GENES = [
    "ERBB2", "EGFR", "AKT1", "AKT2", "MAPK1", "MAPK3",
    "MYC", "TP53", "PIK3CA", "PTEN", "CDK4", "CDK6",
    "CCND1", "BRCA1", "ESR1", "PGR", "VEGFA", "MTOR",
    "RB1", "E2F1", "BCL2", "STAT3", "FOXO1", "FOXO3",
    "GSK3B", "CDKN2A", "PCNA", "MCM2", "GATA3", "RAS",
]

# Pathway annotations for key genes
GENE_PATHWAY_MAP = {
    "ERBB2": "RTK Signaling", "EGFR": "RTK Signaling",
    "AKT1": "PI3K/AKT", "AKT2": "PI3K/AKT",
    "PIK3CA": "PI3K/AKT", "PTEN": "Tumor Suppressor",
    "MAPK1": "MAPK", "MAPK3": "MAPK", "RAS": "MAPK",
    "MYC": "Proliferation", "CCND1": "Cell Cycle",
    "CDK4": "Cell Cycle", "CDK6": "Cell Cycle",
    "RB1": "Cell Cycle", "E2F1": "Cell Cycle",
    "CDKN2A": "Cell Cycle", "PCNA": "DNA Replication",
    "MCM2": "DNA Replication", "BRCA1": "DNA Repair",
    "TP53": "Apoptosis", "BCL2": "Apoptosis",
    "VEGFA": "Angiogenesis", "MTOR": "mTOR",
    "STAT3": "JAK/STAT", "ESR1": "ER Signaling",
    "PGR": "ER Signaling", "GATA3": "Transcription",
    "FOXO1": "Transcription", "FOXO3": "Transcription",
    "GSK3B": "PI3K/AKT",
}

# IC50-to-inhibition conversion thresholds (µM)
IC50_SENSITIVE_THRESHOLD = 1.0   # below = sensitive
IC50_RESISTANT_THRESHOLD = 10.0  # above = resistant

# Network propagation settings
CASCADE_ITERATIONS_DEFAULT = 5
CASCADE_DAMPING = 0.1        # step size per iteration
STRING_CONFIDENCE_MIN = 0.7  # minimum edge confidence score
