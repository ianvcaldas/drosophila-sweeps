from collections.abc import Mapping
from numbers import Number

import sys
import json

import numpy as np
import scipy.stats as stats

sys.path.append(snakemake.scriptdir + "/../..")
from utils.project_parameters import default_simulation_parameters


def draw_from_distribution(distr_dict):
    if "lower" in distr_dict:
        lower = float(distr_dict["lower"])
    if "upper" in distr_dict:
        upper = float(distr_dict["upper"])
    distr = distr_dict["distribution"]
    if distr == "choice":
        value = np.random.choice(distr_dict["options"])
    elif distr == "uniform":
        value = stats.uniform(lower, upper - lower).rvs()
    elif distr == "integer-range":
        value = stats.randint(lower, upper + 1).rvs()
    elif distr == "log-uniform":
        log_value = stats.uniform(
            np.log10(lower), np.log10(upper) - np.log10(lower),
        ).rvs()
        value = np.power(10, log_value)
    else:
        raise Exception(f"Distribution not understood:\n{distr_dict}")
    return value


def number_or_string(x):
    try:
        result = float(x)
    except ValueError as err:
        if isinstance(x, str):
            result = x
        else:
            raise err
    return result


def instantiate_parameter(key, value):
    if isinstance(value, Mapping):
        instantiated = draw_from_distribution(value)
    else:
        instantiated = value
    try:
        json_friendly_value = number_or_string(instantiated)
    except Exception as err:
        raise ValueError(
            f"Could not instantiate {key} from:\n{type(value)}\n{value}. Got this error:\n{err}"
        )
    return json_friendly_value


def instantiate_demography(events):
    resolved_events = []
    for event in events:
        resolved = {
            key: instantiate_parameter(key, value) for key, value in event.items()
        }
        resolved_events.append(resolved)
    return resolved_events


params = snakemake.config["simulation-parameters"]

results = dict()
for key, value in params.items():
    if key == "demography":
        instantiated = instantiate_demography(value)
    else:
        instantiated = instantiate_parameter(key, value)
    results[key] = instantiated

# We always want to start every parameter set with a seed.
# If not given explicitly in the desired parameters, we make one right here.
if "seed" not in results:
    results["seed"] = np.random.randint(1, 2 ** 32)

# Project-wide parameters are added now if not overwritten by the specific config file.
for key, value in default_simulation_parameters.items():
    if key not in results:
        results[key] = number_or_string(value)

# Determine sweep regime based on what parameters were specified.
if "selection-coefficient" in results and "adaptive-mutation-rate" in results:
    regime = "rnm"
elif "selection-coefficient" in results and "frequency-at-selection" in results:
    regime = "sgv"
elif "selection-coefficient" in results:
    regime = "hard"
else:
    regime = "neutral"

results["regime"] = regime

with open(snakemake.output["params_file"], "w") as f:
    json.dump(results, f, sort_keys=True, indent=2)
