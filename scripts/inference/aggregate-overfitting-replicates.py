
from pathlib import Path
import pandas as pd

sys.path.append(snakemake.scriptdir + "/../..")
from utils.prepare_data import read_data, save_data

dfs = [
    read_data(f).assign(replicate=Path(f).stem)
    for f in snakemake.input
]

df = pd.concat(dfs)
save_data(df, snakemake.output["aggregated"])
