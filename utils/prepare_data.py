import numpy as np
import pandas as pd

def read_data(path):
    return pd.read_table(path, dtype={"swept_frequencies": str})

def save_data(df, path):
    df.to_csv(path, index=False, sep='\t')

def balance_log_sel_strength(df, seed=None):
    results = (
        df
        .assign(log_selection_coefficient=np.log10(df.selection_coefficient))
    )
    return results

def balance_sweep_mode(df, seed=None):
    results = (
        df
        .groupby("regime")
        .sample(n=df.regime.value_counts().min(), random_state=seed)
    )
    return results

def balance_hard_vs_soft(df, seed=None):
    kinds = {"hard": "hard", "rnm": "soft", "sgv": "soft"}
    results = df.assign(regime_kind=[kinds[x] for x in df.regime])
    results = (
        results
        .groupby("regime_kind")
        .sample(n=results.regime_kind.value_counts().min(), random_state=seed)
    )
    return results

def balance_rnm_vs_sgv(df, seed=None):
    results = df.loc[df.regime != "hard"]
    results = (
        results
        .groupby("regime")
        .sample(n=df.regime.value_counts().min(), random_state=seed)
    )
    return results

balancing_functions = {
    "log-sel-strength": balance_log_sel_strength,
    "sweep-mode": balance_sweep_mode,
    "hard-vs-soft": balance_hard_vs_soft,
    "rnm-vs-sgv": balance_rnm_vs_sgv,
}

# The actual columns in the DataFrames corresponding to the true labels for each inference target.
# These final columns may only appear after the balancing functions have been applied, I
# don't know yet.
target_columns = {
    "log-sel-strength": "log_selection_coefficient",
    "sweep-mode": "regime",
    "hard-vs-soft": "regime_kind",
    "rnm-vs-sgv": "regime",
}
