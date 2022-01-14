import numpy as np
import pandas as pd


def read_data(path):
    return pd.read_table(path, dtype={"swept_frequencies": str})


def save_data(df, path):
    df.to_csv(path, index=False, sep="\t")


def add_regime_kind_column(df):
    kinds = {"hard": "hard", "rnm (true)": "soft", "sgv (true)": "soft"}
    results = df.assign(regime_kind=[
        kinds[x]
        if x in kinds
        else None
        for x in df.sweep_mode])
    return results


def balance_hard_vs_soft(df, seed=None):
    results = add_regime_kind_column(df)
    results = results.groupby("regime_kind").sample(
        n=results.regime_kind.value_counts().min(), random_state=seed
    )
    return results


def balance_rnm_vs_sgv(df, seed=None):
    results = df.loc[df.sweep_mode != "hard"]
    results = results.groupby("sweep_mode").sample(
        n=df.sweep_mode.value_counts().min(), random_state=seed
    )
    return results


balancing_functions = {
    "log-sel-strength": None,
    "sweep-mode": None,
    "hard-vs-soft": balance_hard_vs_soft,
    "rnm-vs-sgv": balance_rnm_vs_sgv,
}

# The actual columns in the DataFrames corresponding to the true labels for each inference target.
target_columns = {
    "log-sel-strength": "log_selection_coefficient",
    "sweep-mode": "sweep_mode",
    "hard-vs-soft": "regime_kind",
    "rnm-vs-sgv": "sweep_mode",
}
