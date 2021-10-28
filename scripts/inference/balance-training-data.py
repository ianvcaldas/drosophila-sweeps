
import pandas as pd

sys.path.append(snakemake.scriptdir + "/../..")
from utils.prepare_data import balancing_functions, read_data, save_data

balance = balancing_functions[snakemake.wildcards["target"]]
for case in ["training", "validation"]:
    dataset = read_data(snakemake.input[case])
    result = balance(dataset)
    save_data(result, snakemake.output[f"balanced_{case}"])
