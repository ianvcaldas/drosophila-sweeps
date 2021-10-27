
import pandas as pd

sys.path.append(snakemake.scriptdir + "/../..")
from utils.prepare_data import balancing_functions, read_data, save_data

balance = balancing_functions[snakemake.wildcards["target"]]
simulation_parameters = read_data(snakemake.input["sim_params"])
result = balance(simulation_parameters)
save_data(result, snakemake.output["balanced_params"])
