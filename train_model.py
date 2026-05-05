import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, mean_squared_error
import joblib

print("Loading merged dataset...")
df = pd.read_csv("data/merged_dataset.csv")

print("Dataset shape:", df.shape)

print("Keeping only numeric columns...")
df = df.select_dtypes(include=["number"])

if "LN_IC50" not in df.columns:
    raise Exception("Target column LN_IC50 not found")

print("Preparing features and target...")
X = df.drop(columns=["LN_IC50"])
y = df["LN_IC50"]

print("Splitting data...")
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)

print("Training Random Forest model...")

model = RandomForestRegressor(
    n_estimators=30,
    max_depth=12,
    random_state=42,
    n_jobs=-1,
    verbose=1
)

model.fit(X_train, y_train)

print("Evaluating model...")
preds = model.predict(X_test)

r2 = r2_score(y_test, preds)
rmse = mean_squared_error(y_test, preds) ** 0.5

print("R2 Score:", r2)
print("RMSE:", rmse)

print("Saving trained model...")
joblib.dump(model, "model.pkl")

print("Training complete. Model saved as model.pkl")