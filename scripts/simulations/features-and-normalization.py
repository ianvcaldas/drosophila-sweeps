from pathlib import Path
import json

import numpy as np

sys.path.append(snakemake.scriptdir + "/../..")
from utils.project_parameters import default_summary_statistics, summary_statistic_order
from utils.feature_calculation import (
    get_windows,
    calculate_features,
    get_normalization_stats,
    normalize_features,
)


NAME = snakemake.wildcards["sim_id"]
NORM_STATS = snakemake.input["normalization_stats"]

ms_file = snakemake.input["ms"]
npy_file = snakemake.output["npy"]
feature_file = snakemake.output["features"]
stats_file = snakemake.output["stats"]
npy_log_file = snakemake.output["npy_log"]

with open(snakemake.input["params_file"], "r") as f:
    params = json.load(f)

window_sizes, center_pos_dict = get_windows(
    int(params["locus-size"]),
    int(params["data-dimension"]),
    start_pos=1,
    smallest_window=int(params["smallest-window"]),
)

calculate_features(
    ms_file=ms_file,
    summary_statistics=default_summary_statistics,
    center_pos=center_pos_dict,
    window_sizes=window_sizes,
    last_position=params["locus-size"],
    output_file=feature_file,
    output_stats=stats_file,
)


norm_stats_dict = get_normalization_stats(NORM_STATS, summary_statistic_order)
normalized = normalize_features(
    feature_file,
    summary_statistic_order,
    norm_stats_dict,
    reshape=True,
    convert_to_uint8=True,
    log_scale=False,
)

for data_id, normalized_array in normalized.items():
    np.save(npy_file, normalized_array)

# Also do log-scale normalization
normalized_log = normalize_features(
    feature_file,
    summary_statistic_order,
    norm_stats_dict,
    reshape=True,
    convert_to_uint8=True,
    log_scale=True,
)
for data_id, normalized_array in normalized_log.items():
    np.save(npy_log_file, normalized_array)
