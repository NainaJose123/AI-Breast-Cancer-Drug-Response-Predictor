import pandas as pd

def load_expression(path):
    print("Loading CCLE expression...")
    df = pd.read_csv(path)

    if "Unnamed: 0" in df.columns:
        df = df.rename(columns={"Unnamed: 0": "DepMap_ID"})

    # Keep only first 1000 genes to avoid memory crash
    cols = ["DepMap_ID"] + list(df.columns[1:1001])
    df = df[cols]

    print("Expression shape after reduction:", df.shape)
    return df


def load_ic50(path):
    print("Loading GDSC IC50...")
    df = pd.read_excel(path)
    return df


def load_metadata(path):
    print("Loading sample metadata...")
    df = pd.read_csv(path)
    return df


def merge_datasets(expr, ic50, meta):
    print("Merging datasets...")

    meta_small = meta[["DepMap_ID", "stripped_cell_line_name"]]

    merged = expr.merge(meta_small, on="DepMap_ID", how="inner")

    if "CELL_LINE_NAME" in ic50.columns:
        ic50 = ic50.rename(columns={"CELL_LINE_NAME": "stripped_cell_line_name"})

    merged = merged.merge(ic50, on="stripped_cell_line_name", how="inner")

    return merged


def main():
    expr = load_expression("data/CCLE_expression.csv")
    ic50 = load_ic50("data/GDSC2_fitted_dose_response_27Oct23.xlsx")
    meta = load_metadata("data/sample_info.csv")

    merged = merge_datasets(expr, ic50, meta)

    print("Final dataset shape:", merged.shape)

    merged.to_csv("data/merged_dataset.csv", index=False)
    print("Saved merged dataset.")


if __name__ == "__main__":
    main()