from pathlib import Path
import tarfile
import shutil

import pandas as pd

uuids = (
    pd.read_table(snakemake.input["params_file"], dtype={"data_id": str}).set_index(
        "data_id"
    )
)["uuid"]

original_npys = [Path(f) for f in snakemake.input["npy_files"]]
unique_folder = Path(str(snakemake.params["npy_folder"]) + "-unique")
unique_folder.mkdir(exist_ok=True)
for original in original_npys:
    unique = unique_folder / (uuids.loc[original.stem] + ".npy")
    shutil.copy(original, unique)

with tarfile.open(snakemake.output["tar_file"], "w") as f:
    f.add(unique_folder)
