# src/model_trainer.py
import os
import logging
import joblib
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import classification_report
from xgboost import XGBClassifier

from src.config import DRUGS, KEY_GENES

logger = logging.getLogger(__name__)

MODEL_DIR = os.path.join(os.path.dirname(__file__), "..", "models")
MODEL_PATH = os.path.join(MODEL_DIR, "xgb_drug_sensitivity.joblib")
DRUG_ENCODER_PATH = os.path.join(MODEL_DIR, "drug_label_encoder.joblib")
FEATURE_NAMES_PATH = os.path.join(MODEL_DIR, "feature_names.joblib")


def build_feature_vector(pre_state: dict, post_state: dict, drug_key: str, gene_list: list[str], drug_keys: list[str]) -> np.ndarray:
    pre_vec = np.array([pre_state.get(g, 0.0) for g in gene_list], dtype=np.float32)
    post_vec = np.array([post_state.get(g, 0.0) for g in gene_list], dtype=np.float32)
    delta_vec = post_vec - pre_vec

    drug_onehot = np.zeros(len(drug_keys), dtype=np.float32)
    if drug_key in drug_keys:
        drug_onehot[drug_keys.index(drug_key)] = 1.0

    return np.concatenate([pre_vec, post_vec, delta_vec, drug_onehot])


def generate_synthetic_training_data(n_samples_per_drug: int = 200, seed: int = 42) -> tuple[np.ndarray, np.ndarray]:
    from src.network_builder import build_fallback_network
    from src.drug_simulator import simulate_drug

    rng = np.random.RandomState(seed)
    drug_keys = list(DRUGS.keys())
    X_rows, y_rows = [], []

    for drug_key, drug_conf in DRUGS.items():
        for _ in range(n_samples_per_drug):
            base_expr = rng.uniform(3.0, 10.0, len(KEY_GENES))
            if "ERBB2" in drug_conf["target_genes"] and rng.random() < 0.5:
                if "ERBB2" in KEY_GENES:
                    base_expr[KEY_GENES.index("ERBB2")] = rng.uniform(8.5, 12.0)
            patient_series = pd.Series(base_expr, index=KEY_GENES)

            G = build_fallback_network(KEY_GENES, patient_series)
            pre_state, post_state = simulate_drug(G, drug_conf)

            target_suppression = np.mean([
                max(0, pre_state.get(t, 5.0) - post_state.get(t, 5.0))
                for t in drug_conf["target_genes"] if t in pre_state
            ])
            label = 1 if target_suppression > 1.5 else 0

            feat = build_feature_vector(pre_state, post_state, drug_key, KEY_GENES, drug_keys)
            X_rows.append(feat)
            y_rows.append(label)

    return np.array(X_rows), np.array(y_rows)


def train_model(X: np.ndarray, y: np.ndarray, n_estimators: int = 200, max_depth: int = 5, learning_rate: float = 0.05, cv_folds: int = 5) -> XGBClassifier:
    os.makedirs(MODEL_DIR, exist_ok=True)

    clf = XGBClassifier(
        n_estimators=n_estimators,
        max_depth=max_depth,
        learning_rate=learning_rate,
        subsample=0.8,
        colsample_bytree=0.8,
        use_label_encoder=False,
        eval_metric="logloss",
        random_state=42,
        n_jobs=-1,
    )

    # Cross-validation
    cv_scores = cross_val_score(clf, X, y, cv=StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=42), scoring="roc_auc")
    logger.info("Cross-validation AUC: %.3f ± %.3f", cv_scores.mean(), cv_scores.std())

    clf.fit(X, y)
    logger.info("Model trained on %d samples.", len(X))

    # Handle single-class situation safely
    unique_classes = np.unique(y)
    if len(unique_classes) == 1:
        logger.warning("Only one class present in y. Skipping detailed classification report.")
        print(f"All samples belong to class: {unique_classes[0]}")
    else:
        print(classification_report(y, clf.predict(X), target_names=["Resistant", "Sensitive"]))

    # Save model
    drug_keys = list(DRUGS.keys())
    joblib.dump(clf, MODEL_PATH)
    joblib.dump(drug_keys, DRUG_ENCODER_PATH)

    feature_names = [f"{prefix}_{g}" for prefix in ["pre", "post", "delta"] for g in KEY_GENES] + [f"drug_{k}" for k in drug_keys]
    joblib.dump(feature_names, FEATURE_NAMES_PATH)
    logger.info("Model saved to %s", MODEL_PATH)

    return clf


def load_model() -> tuple[XGBClassifier, list[str]]:
    if not os.path.exists(MODEL_PATH):
        raise FileNotFoundError(f"No trained model found at {MODEL_PATH}. Run: python -m src.model_trainer")
    clf = joblib.load(MODEL_PATH)
    drug_keys = joblib.load(DRUG_ENCODER_PATH)
    return clf, drug_keys


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    logger.info("Generating synthetic training data ...")
    X, y = generate_synthetic_training_data(n_samples_per_drug=300)
    logger.info("Training XGBoost model (%d samples, %d features) ...", *X.shape)
    train_model(X, y)
    logger.info("Training complete.")