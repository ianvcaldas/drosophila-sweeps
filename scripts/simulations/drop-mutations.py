import json

import numpy as np
import msprime
import pyslim

sys.path.append(snakemake.scriptdir + "/../..")
from utils.conversion import treeseq_to_ms


def sample_genomes(ts, num_genomes):
    # If there are less samples than we want, then take all the available ones.
    if num_genomes > ts.get_sample_size():
        return ts
    samples = np.random.choice(ts.get_samples(), size=num_genomes, replace=False)
    return ts.simplify(samples)


def get_slim_mut_freq(ts, slim_id):
    # Start by finding any tskit mutation record with the right slim_id, with stacked slim mutations or not
    mut = None
    for m in ts.mutations():
        if slim_id in m.derived_state.split(","):
            mut = m
            break
    if mut is None:
        return -1
    # Find the tskit record with the least stacked (more ancestral) version of the slim mutation
    while (mut.parent != -1) and (
        slim_id in ts.mutation(mut.parent).derived_state.split(",")
    ):
        mut = ts.mutation(mut.parent)
    # Find the node that tskit record is at and count the number of samples that are its children,
    # that should be the count of derived alleles of the slim mutation including all stacked versions of it
    position = ts.site(mut.site).position
    tree = ts.at(position)
    return tree.num_samples(mut.node)


def report_slim_mutations(ts):
    muts = []
    seen_slim_ids = []
    for mut in ts.mutations():
        slim_ids = mut.derived_state.split(",")
        for stacked_id, slim_id in enumerate(slim_ids):
            if slim_id in seen_slim_ids:
                continue
            counts = get_slim_mut_freq(ts, slim_id)
            result = {
                "slim_id": slim_id,
                "derived_count": counts,
                "frequency": counts / ts.num_samples,
            }
            result.update(mut.metadata["mutation_list"][stacked_id])
            muts.append(result)
            seen_slim_ids.append(slim_id)
    return muts


def msprime_metrics_premutation(ts, regime):
    """Calculate some statistics about a TreeSequence without neutral mutations and with
    the adaptive mutations still in there."""
    alleles_info = report_slim_mutations(ts)
    if regime == "sgv":
        mut_type = 2
    else:
        mut_type = 1
    derived_counts = [
        info["derived_count"]
        for info in alleles_info
        if info["mutation_type"] == mut_type
    ]
    return {
        "sample_num_adaptive_copies": len(derived_counts),
        "sample_adaptive_frequency": sum(derived_counts) / ts.sample_size,
    }


def clear_msprime_mutations(ts):
    """
    Empty out a msprime tree sequence from mutations and their sites. SLiM generates TSs
    with multi-character alleles and msprime doesn't support iterating haplotypes in that
    case.
    """
    tables = ts.dump_tables()
    tables.mutations.clear()
    tables.sites.clear()
    return tables.tree_sequence()


def drop_mutations(ts, mutation_rate, seed):
    mutated = msprime.mutate(ts, mutation_rate, random_seed=seed, keep=True)
    return mutated


def msprime_metrics_postmutation(ts):
    """Calculate some statistics about a TreeSequence with neutral mutations and without
    the adaptive mutations."""
    return {
        "pairwise_div": ts.diversity().item(),
        "num_snps": ts.num_sites,
        "actual_sample_size": ts.sample_size,
    }


# Load data

with open(snakemake.input["params_file"], "r") as f:
    params = json.load(f)
trees = pyslim.load(snakemake.input["slim_output"])

# Sample, get statistics, and drop mutations

trees = sample_genomes(trees, int(params["sample-size"]))
msprime_metrics = dict()
msprime_metrics.update(msprime_metrics_premutation(trees, params["regime"]))
trees = clear_msprime_mutations(trees)
trees = drop_mutations(trees, params["mutation-rate"], params["seed"])
msprime_metrics.update(msprime_metrics_postmutation(trees))

# Create output

trees.dump(snakemake.output["trees_file"])

with open(snakemake.output["metrics_file"], "w") as f:
    string = "\n".join([f"{key} = {value}" for key, value in msprime_metrics.items()])
    f.write(string)

treeseq_to_ms(
    trees,
    ms_file=snakemake.output["ms_file"],
    normalize_positions=False,
    header_line=f"SLiM with TreeSeq, simulation ID {snakemake.wildcards['sim_id']}",
    seed=params["seed"],
)
