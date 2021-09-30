from pathlib import Path
import sys

import numpy as np
import pandas as pd

from .conversion import ms_to_numpy


def get_windows(locus_size, dimension, start_pos=1, smallest_window=1000):
    """Get logarithmically distributed window sizes and their locations, for the given parameters.
    start_pos offsets the center of every window.
    """
    # Window sizes are logarithmically distributed.
    window_sizes = np.exp(
        np.linspace(
            np.log(smallest_window),
            np.log(locus_size / ((dimension + 1) / 2)),
            dimension,
        )
    )
    window_sizes = [np.round(i) for i in window_sizes]
    # Check that the parameters produce meaningful window sizes.
    if (
        len(np.unique(window_sizes)) != len(window_sizes)
        or window_sizes[0] > window_sizes[-1]
        or smallest_window > locus_size
    ):
        lines = [
            f"Your combination of parameters does not produce meaningful window sizes.",
            f"Locus size: {locus_size}, dimension: {dimension}, starting position: {start_pos}, smallest window: {smallest_window}.",
            "Window sizes:",
        ]
        lines.extend([str(ws) for ws in window_sizes])
        raise Exception("\n".join(lines))
    # The positions for each window size depend on that window size.
    center_pos_dict = {
        ws: (
            np.linspace(
                locus_size / 2 - (dimension - 1) * ws / 4,
                locus_size / 2 + (dimension - 1) * ws / 4,
                dimension,
            )
            - 1
            + int(start_pos)
        ).tolist()
        for ws in window_sizes
    }
    return window_sizes, center_pos_dict


def calculate_features(
    ms_file,
    summary_statistics,
    center_pos=None,
    window_sizes=None,
    snp_window_size=None,
    snp_window_step=None,
    last_position=None,
    output_file=None,
    output_stats=None,
):
    """Calculates the desired population genetic summary statistics for the provided ms
    file, using the given window sizes and dimensions.
    The output is a long table, written either to stdout or to a file.
    """
    ms_file = Path(ms_file)
    sim_stats = dict()
    if output_file is None:
        output = sys.stdout
    else:
        output = open(output_file, "w")
    output.write(f"data_id\tfeature\tcenter_pos_ix\tcenter_pos\tws\tvalue\n")
    for stat in summary_statistics:
        sim_stats[stat] = dict()
        sim_stats[stat]["min"] = np.inf
        sim_stats[stat]["max"] = -np.inf
    snp_coordinates, genotypes = ms_to_numpy(ms_file)
    if last_position is None:
        # This allows for 0-1 or raw base-pair coordinates.
        try:
            last_position = max(snp_coordinates[-1], 1)
        except IndexError:
            s = "Can't determine locus size; you probably have 0 segregating sites and no last_position was passed to calculate_features."
            raise Exception(s)
    for stat in summary_statistics:
        # these if statements are a poor way to check if the user wants
        # a bp window or a snp window.
        if center_pos is not None:
            results = stat.run_bp(
                snp_coordinates, genotypes, center_pos, window_sizes, last_position
            )
        elif snp_window_size is not None:
            results = stat.run_snp(
                snp_coordinates, genotypes, snp_window_size, snp_window_step
            )
        else:
            raise ValueError("Bad window specifications.")
        for ws, ws_res in results.items():
            for cp_ix, (cp, res) in enumerate(ws_res.items()):
                if res > sim_stats[stat]["max"]:
                    sim_stats[stat]["max"] = res
                if res < sim_stats[stat]["min"]:
                    sim_stats[stat]["min"] = res
                output.write(
                    f"{ms_file.stem}\t{stat.name}\t{cp_ix}\t{cp}\t{ws}\t{res}\n"
                )
    if output_file is not None:
        output.close()
    if output_stats is not None:
        outs = open(output_stats, "w")
    else:
        outs = sys.stdout
    for stat in sim_stats:
        for key, value in sim_stats[stat].items():
            outs.write(f"{stat}\t{key}\t{value}\n")
    if output_stats is not None:
        outs.close()


def get_normalization_stats(stats_file, feat_list):
    """Reads normalization stats file, modifies it with appropriate biological priors, and
    returns a dictionary. feat_list is the list of population genetics summary statistics
    in the right order they were calculated in.
    """
    stats = pd.read_table(stats_file, index_col=0)
    result = dict()
    for row in stats.iterrows():
        name = row[0]
        # If there's a statistic in the stats file that you don't want to keep.
        if name not in feat_list:
            continue
        this_min = row[1]["feat_min"]
        this_max = row[1]["feat_max"]

        result[name] = dict()
        result[name]["ix"] = feat_list.index(name)
        result[name]["min"] = this_min
        result[name]["max"] = this_max

    # Hard-coded minimums and maximums
    result["H1"]["min"] = 0
    result["H1"]["max"] = 1
    result["H12"]["min"] = 0
    result["H12"]["max"] = 1
    result["H2/H1"]["min"] = 0
    result["H2/H1"]["max"] = 1

    # For features with NaN, what value do we want to use instead?
    result["H2/H1"]["nan_fill"] = 0
    result["taj_D"]["nan_fill"] = result["taj_D"]["min"]

    return result


def features_long_to_matrix(data_id, feature_df, feature_order, return_df=False):
    """Converts a list of feature values in long form into a stacked 3d matrix.
    Only does this for the requested data_id analysis, if there's more than one in the
    features file. feature_df is the long-form dataframe, feature_order is the desired
    order of features in the resulting stack.
    
    If return_df is True, get a dictionary of individual 2d feature matrices as
    dataframes. Otherwise, return as a 3d numpy array (without labels for position and
    window size).
    """
    features_of_data = feature_df.loc[feature_df["data_id"] == data_id]
    feature_2d_matrices = []
    for feature in feature_order:
        feature_values = features_of_data.loc[features_of_data["feature"] == feature]
        feature_2d_matrices.append(
            (
                feature_values.pivot(
                    index="ws", columns="center_pos_ix", values="value"
                ).sort_index(ascending=False)
            )
        )
    if return_df:
        result = {
            feature: feature_2d_matrices[ix] for ix, feature in enumerate(feature_order)
        }
        return result
    else:
        return np.stack([df.values for df in feature_2d_matrices])


def normalize_feature_array(
    features_np,
    feature_order,
    normalization_stats,
    log_scale=False,
    convert_to_255=True,
):
    """Normalize values in a 3d numpy array of features.
    If log_scale, normalize on a log instead of a linear scale.
    If convert_to_255, transpose values to range 0-255.
    """
    result = features_np[:]  # Copy of original array
    for ix, feature in enumerate(feature_order):
        x = result[ix]
        # Fill NaNs
        if "nan_fill" in normalization_stats[feature]:
            x[np.where(np.isnan(x))] = normalization_stats[feature]["nan_fill"]
        a = normalization_stats[feature]["min"]
        b = normalization_stats[feature]["max"]
        # Truncate values outside the normalization distribution
        x[x > b] = b
        x[x < a] = a
        if log_scale:
            x = np.log1p(abs(x)) * np.sign(x)
            a = np.log1p(abs(a)) * np.sign(a)
            b = np.log1p(abs(b)) * np.sign(b)
        if convert_to_255:
            result[ix] = (x - a) * (255 / (b - a))
        else:
            result[ix] = x
    return result


def normalize_features(
    feature_file,
    feature_order,
    normalization_stats,
    reshape=True,
    convert_to_uint8=True,
    log_scale=False,
):
    """Normalize features for each data_id in feature_df.
    If reshape is True, convert from a ZxHxW array to a HxWxZ array.
    """
    feature_df = pd.read_table(feature_file)
    normalized = dict()
    for data_id in feature_df.data_id.unique():
        unnormalized_matrix = features_long_to_matrix(
            data_id, feature_df, feature_order, return_df=False
        )
        normalized_matrix = normalize_feature_array(
            unnormalized_matrix,
            feature_order,
            normalization_stats,
            convert_to_255=True,
            log_scale=log_scale,
        )
        if reshape:
            normalized_matrix = np.transpose(normalized_matrix, axes=[1, 2, 0])
        if convert_to_uint8:
            normalized_matrix = normalized_matrix.round().astype(np.uint8)
        normalized[data_id] = normalized_matrix
    return normalized
