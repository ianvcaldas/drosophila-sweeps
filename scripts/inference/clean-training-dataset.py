
sys.path.append(snakemake.scriptdir + "/../..")
from utils.prepare_data import read_data, save_data

raw = read_data(snakemake.input["sim_params"])

# Remove failed and neutral simulations
df = raw.loc[(raw.simulation_status == 'ok') & (raw.regime != 'neutral')]

NUM_SIMS = snakemake.config["sims_per_regime"]

if not all(df.regime.value_counts() > NUM_SIMS):
    raise Exception(f"One or more selection regimes did not have {NUM_SIMS} simulations:\n{df.regime.value_counts()}")
    
df = df.groupby('regime').sample(n=NUM_SIMS, random_state=snakemake.params["random_seed"])

save_data(df, snakemake.output["cleaned"])
