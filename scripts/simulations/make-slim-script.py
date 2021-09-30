from pathlib import Path
from collections import defaultdict

import json

import numpy as np


with open(snakemake.input["params_file"], "r") as f:
    specified_params = json.load(f)

# Missing parameter values get a placeholder -1.
params = defaultdict(lambda: -1)
params.update(specified_params)

replacements = {
    # SLiM configuration
    "INITIAL_SEED": params["seed"],
    "LAST_GENERATION": int(params["last-generation"]),
    "SIMPLIFICATION_INTERVAL": int(params["simplification-interval"]),
    "RESTARTS_LIMIT": int(params["restarts-limit"]),
    # Biology
    "LOCUS_SIZE": int(params["locus-size"] - 1),
    "RECOMBINATION_RATE": params["recombination-rate"],
    "DOMINANCE_COEFFICIENT": params["dominance-coefficient"],
    "FREQUENCY_AT_SAMPLING": params["frequency-at-sampling"],
    "SELECTION_COEFFICIENT": params["selection-coefficient"],
    "SELECTION_COORDINATE": int(params["selection-coordinate"]),
    "SELECTION_GENERATION": int(params["selection-generation"]),
    "ADAPTIVE_MUTATION_RATE": params["adaptive-mutation-rate"],
    "FREQUENCY_AT_SELECTION": params["frequency-at-selection"],
    # Files
    "CURRENT_WORKDIR": Path(".").resolve(),
    "BURNIN_TREES_FILE": snakemake.input["burnin"],
    "OUTPUT_TREES_FILE": str(
        Path(snakemake.params["outdir"])
        / (snakemake.wildcards["sim_id"] + "_slim.trees")
    ),
    "SAVE_TREES_FILE": str(
        Path(snakemake.params["outdir"])
        / (snakemake.wildcards["sim_id"] + "_savefile.trees")
    ),
    "METRICS_OUTPUT_FILE": str(
        Path(snakemake.params["outdir"])
        / (snakemake.wildcards["sim_id"] + "_slim-metrics.txt")
    ),
}
# Special adjustments for RNM sweeps.
if params["regime"] == "rnm":
    s = params["selection-coefficient"]
    position = params["selection-coordinate"]
    window = params["selection-region-size"]
    # We have to set the dominance coefficient to match the selection coefficient, in case
    # an individual gets two separate mutations in their two genomes. Check section 10.5.1
    # of the SLiM Manual for details.
    replacements["DOMINANCE_COEFFICIENT"] = (np.sqrt(1 + s) - 1) / s
    # Coordinates of adaptive region for SLiM breakpoints of adaptive mutation rate.
    replacements["END_NEUTRAL_REGION"] = int(position - max(window // 2, 1))
    replacements["END_ADAPTIVE_REGION"] = int(position + window // 2)


def find_template(name):
    template = None
    for template_option in snakemake.input["slim_templates"]:
        if Path(template_option).stem == name:
            template = template_option
            break
    if template is None:
        raise FileNotFoundError(
            f"Could not find template {name}. Template options were:\n{snakemake.input['slim_templates']}"
        )
    return template


def load_template(template_file, replacements):
    with open(template_file, "r") as f:
        text = "".join(line for line in f)
    for key, replacer in replacements.items():
        text = text.replace(key, str(replacer))
    return text


def get_demography_text(event):
    if event["kind"] == "size-change":
        demog_replacements = dict()
        demog_replacements["GENERATION"] = int(event["generation"])
        demog_replacements["DIPLOID_SIZE"] = int(event["diploid-size"])
        template = find_template(event["kind"])
        demog_text = load_template(template, demog_replacements)
        return demog_text
    elif event["kind"] == "bottleneck":
        # A special event that resolves into two size-change events.
        population_crash = {
            "kind": "size-change",
            "generation": int(event["generation"]),
            "diploid-size": int(event["crash-size"]),
        }
        population_recovery = {
            "kind": "size-change",
            "generation": int(event["generation"] + event["duration"]),
            "diploid-size": int(event["recovery-size"]),
        }
        crash_text = get_demography_text(population_crash)
        recovery_text = get_demography_text(population_recovery)
        return crash_text + "\n" + recovery_text
    else:
        raise ValueError(f"Demography event not recognized:\n{event}")


template = find_template(params["regime"])
script_text = load_template(template, replacements)
if "demography" in params:
    for event in params["demography"]:
        demography_text = get_demography_text(event)
        script_text += "\n" + demography_text
with open(snakemake.output["slim_script"], "w") as f:
    f.write(script_text)
