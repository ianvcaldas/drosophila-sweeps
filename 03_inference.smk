
configfile: '03_config.yaml'
report: 'captions/inference.rst'

# Think about how this interacts with the cluster, but do it later.
# workdir: config["workdir"]

# Every rule should have shadow: copy-minimal.

rule all:
    input:
        data = expand(
            "output/simulations/{sim_id}/{result}",
            sim_id=config["sim_folders"],
            result=["data.tar", "parameters.tsv"]
        )

rule combine_simulation_tasks:
    output:
        data = "output/simulations/{sim}/data.tar",
        parameters = "output/simulations/{sim}/parameters.tsv"
    conda: "envs/simulate.yaml"
    notebook: "notebooks/combine-simulations.py.ipynb"
