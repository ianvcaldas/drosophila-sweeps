from sklearn.model_selection import train_test_split

sys.path.append(snakemake.scriptdir + "/../..")
from utils.prepare_data import read_data, save_data, target_columns

try:
    seed = snakemake.params["random_seed"]
except AttributeError:  # No random seed given
    seed = None

raw = read_data(snakemake.input["sim_params"])

# Remove neutral simulations; remove SGV and RNM sweeps that ended up hard-like.
df = raw.loc[raw.sweep_mode.isin(("hard", "rnm (true)", "sgv (true)"))]
# Downsample to the desired amount of simulations per sweep mode.
# If there are not enough simulations, just use however many there are.
sims_per_regime = min(
    df.sweep_mode.value_counts().min(), snakemake.config["sims_per_regime"]
)
downsampled = df.groupby("sweep_mode").sample(n=sims_per_regime, random_state=seed)

train, valid = train_test_split(
    downsampled,
    test_size=snakemake.config["validation_percentage"],
    random_state=seed,
    shuffle=True,
    stratify=downsampled.sweep_mode,
)

save_data(train, snakemake.output["training"])
save_data(valid, snakemake.output["validation"])
