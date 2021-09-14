import json

import msprime
import pyslim

with open(snakemake.input["params_file"], "r") as f:
    params = json.load(f)

treeseq = msprime.simulate(
    sample_size=int(params["diploid-population-size"] * 2),
    Ne=int(params["diploid-population-size"]),
    length=int(params["locus-size"]),
    mutation_rate=0,
    recombination_rate=params["recombination-rate"],
    random_seed=params["seed"],
)

slim_treeseq = pyslim.annotate_defaults(treeseq, model_type="WF", slim_generation=1)

slim_treeseq.dump(snakemake.output["trees"])
