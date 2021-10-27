import pandas as pd

def read_data(path):
    return pd.read_table(path, dtype={"swept_frequencies": str})

def save_data(df, path):
    df.to_csv(path, index=False, sep='\t')

def balance_log_sel_strength(df):
    return df

def balance_sweep_mode(df):
    return df

def balance_hard_vs_soft(df):
    return df

def balance_rnm_vs_sgv(df):
    return df

balancing_functions = {
    "log-sel-strength": balance_log_sel_strength,
    "sweep-mode": balance_sweep_mode,
    "hard-vs-soft": balance_hard_vs_soft,
    "rnm-vs-sgv": balance_rnm_vs_sgv,
}
