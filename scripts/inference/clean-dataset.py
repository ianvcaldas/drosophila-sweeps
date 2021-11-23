import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype

sys.path.append(snakemake.scriptdir + "/../..")
from utils.prepare_data import read_data, save_data


def get_sweep_mode(df):
    x = df.regime.copy()
    one_adaptive_allele = df.sample_num_adaptive_copies == 1
    failed = df.sample_num_adaptive_copies.isna()
    rnm = df.regime == "rnm"
    sgv = df.regime == "sgv"
    x.loc[one_adaptive_allele & sgv] = "sgv (hard-like)"
    x.loc[one_adaptive_allele & rnm] = "rnm (hard-like)"
    x.loc[~one_adaptive_allele & sgv & ~failed] = "sgv (true)"
    x.loc[~one_adaptive_allele & rnm & ~failed] = "rnm (true)"
    x.loc[~one_adaptive_allele & sgv & failed] = "sgv (failed)"
    x.loc[~one_adaptive_allele & rnm & failed] = "rnm (failed)"
    x.loc[(df.regime == "hard") & failed] = "hard (failed)"
    return x


def add_key_columns(df):
    return df.assign(
        log_selection_coefficient=np.log10(df.selection_coefficient),
        sweep_mode=get_sweep_mode(df),
    )


def get_distribution(s):
    always_unique = ["demography", "seed", "swept_frequencies"]
    if s.name in always_unique:
        x = f"{s.nunique()} unique values."
    elif is_numeric_dtype(s):
        if s.nunique() == 0:
            x = f"{len(s)}x NaN"
        elif s.nunique() == 1:
            x = f"{s.notna().sum()}x {s.loc[s.notna()].unique()[0]}"
        else:
            x = s.describe().to_string()
    elif s.loc[s.notna()].is_unique:
        x = f"{s.nunique()} unique values."
    else:
        if s.nunique() == 1:
            x = f"{s.notna().sum()}x {s.loc[s.notna()].unique()[0]}"
        else:
            x = s.value_counts().to_string()
    return x


def get_columns_report(df):
    lines = []
    for column in sorted(df.columns):
        lines.append(column)
        lines.append(get_distribution(df[column]))
        lines.append("")
    return "\n".join(lines)


def get_df_report(df, name):
    lines = []
    title = f"{name} simulations"
    lines.append(title.upper())
    lines.append("-" * len(title))
    lines.append("")
    lines.append(get_columns_report(df))
    lines.append("")
    return lines


def make_data_report(df):
    lines = []
    constant_cols = df.columns[df.nunique(dropna=False) == 1].tolist()
    title = "CONSTANT PARAMETERS"
    lines.append(title)
    lines.append("-" * len(title))
    lines.append("")
    for col in constant_cols:
        lines.append(col)
        s = df[col]
        lines.append(f"{len(s)}x {s.unique()[0]}")
        lines.append("")
        if col != "sweep_mode":
            df = df.drop(col, axis="columns")
    lines.append("")
    modes = {sm: df.loc[df.sweep_mode == sm] for sm in sorted(df.sweep_mode.unique())}
    lines.extend(get_df_report(df, "All"))
    for mode, this_df in modes.items():
        lines.extend(get_df_report(this_df, mode))
    return "\n".join(lines)


raw = add_key_columns(read_data(snakemake.input["sim_params"]))
successful = raw.loc[raw.simulation_status == "ok"]
failed = raw.loc[raw.simulation_status == "failed"]

save_data(successful, snakemake.output["cleaned_parameters"])
successful_report = make_data_report(successful)
with open(snakemake.output["successful_report"], "w") as f:
    f.write(f"Successful simulations of {snakemake.wildcards['sim_id']}\n\n\n")
    f.write(successful_report)

failed_report = make_data_report(failed)
with open(snakemake.output["failed_report"], "w") as f:
    f.write(f"Failed simulations of {snakemake.wildcards['sim_id']}\n\n\n")
    f.write(failed_report)

with open(snakemake.output["sweep_mode_report"], "w") as f:
    sweep_modes = successful.sweep_mode.value_counts().sort_index().to_string()
    failed_modes = failed.sweep_mode.value_counts().sort_index().to_string()
    f.write(f"Sweep mode report for {snakemake.wildcards['sim_id']}" + "\n\n")
    f.write(sweep_modes)
    f.write("\n\n")
    f.write(failed_modes)
