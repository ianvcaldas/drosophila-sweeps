import json
import uuid
from pathlib import Path

import pandas as pd


def parse_demography(demog):
    """Converts list of dicts in demography specification into a string that can go into a DataFrame."""
    components = []
    for event in demog:
        components.append(",".join(f"{key}={value}" for key, value in event.items()))
    parsed = ";".join(components)
    return parsed


metrics = {
    "slim-metrics": [],
    "msprime-metrics": [],
}


for file in [file for file in snakemake.input if file.endswith("metrics.txt")]:
    with open(file) as f:
        lines = [line.strip().split(" = ") for line in f]
    params_dict = dict()
    for line in lines:
        try:
            value = float(line[1])
        except ValueError:
            value = line[1]
        params_dict[line[0]] = value
    params_dict["data-id"], metric_kind = Path(file).stem.split("_")
    metrics[metric_kind].append(params_dict)

dfs = []
for list_of_dicts in metrics.values():
    frame = pd.DataFrame.from_records(list_of_dicts)
    if "data-id" in frame:
        frame = frame.set_index("data-id")
    dfs.append(frame)

json_params = []
for file in [file for file in snakemake.input if file.endswith("json")]:
    with open(file) as f:
        params_dict = json.load(f)
    if "demography" in params_dict:
        params_dict["demography"] = parse_demography(params_dict["demography"])
    params_dict["data-id"] = Path(file).stem.split("_")[0]
    json_params.append(params_dict)
json_df = pd.DataFrame.from_records(json_params).set_index("data-id")
dfs.append(json_df)


benchmark_files = list(Path(snakemake.params["benchmarks_folder"]).glob("*.tsv"))
benchmark_timings = []
for file in benchmark_files:
    seconds = pd.read_table(file).loc[0, "s"]
    data_id, benchmark_kind = file.stem.split("_")
    benchmark_timings.append((data_id, benchmark_kind, seconds))
benchmark_df = pd.DataFrame.from_records(
    benchmark_timings, columns=["data-id", "benchmark", "seconds"]
).pivot(index="data-id", columns="benchmark", values="seconds")
benchmark_df = benchmark_df.rename(
    {col: "benchmark_" + col for col in benchmark_df}, axis="columns"
)
dfs.append(benchmark_df)

failed_sims = [
    file.stem.split("_")[0]
    for file in Path(snakemake.params["simulations_folder"]).glob("*failed")
]
failed_series = pd.Series("failed", index=failed_sims, name="simulation_status")
failed_series.index.name = "data-id"
dfs.append(failed_series)

df = pd.concat(dfs, axis="columns")
df = df.fillna({"simulation_status": "ok"})
df = df.rename({col: col.replace("-", "_") for col in df.columns}, axis="columns")
df.index.name = "data_id"

# Finally, give a unique ID to each simulation
uuids = [uuid.uuid4() for _ in range(len(df))]
df = df.assign(uuid=uuids)

df.to_csv(snakemake.output["params_output"], index=True, sep="\t")
