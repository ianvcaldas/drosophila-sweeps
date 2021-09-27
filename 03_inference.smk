
configfile: '03_config.yaml'
report: 'captions/inference.rst'

# Think about how this interacts with the cluster, but do it later.
# workdir: config["workdir"]
# Cluster rules might have shadow: copy-minimal.

rule all:
    input:
        data = expand(
            "output/simulation-data/{sim_id}/{result}",
            sim_id=config["sim_folders"],
            result=["data.tar", "parameters.tsv", "info.txt"]
        )

rule combine_simulation_tasks:
    output:
        data = "output/simulation-data/{sim}/data.tar",
        parameters = "output/simulation-data/{sim}/parameters.tsv",
        info = "output/simulation-data/{sim}/info.txt"
    params:
        combine_features = True
    conda: "envs/simulate.yaml"
    notebook: "notebooks/combine-simulations.py.ipynb"
