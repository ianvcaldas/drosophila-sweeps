import pandas as pd

sys.path.append(snakemake.scriptdir + "/../..")
from utils.prepare_data import balancing_functions, read_data, save_data

try:
    seed = snakemake.params["random_seed"]
except AttributeError:
    seed = None

balance = balancing_functions[snakemake.wildcards["target"]]
for case in ["training", "validation"]:
    dataset = read_data(snakemake.input[case])
    if balance is not None:
        result = balance(dataset, seed)
    else:
        result = dataset
    save_data(result, snakemake.output[f"balanced_{case}"])
