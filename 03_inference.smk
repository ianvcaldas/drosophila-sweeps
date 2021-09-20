
configfile: '03_config.yaml'
report: 'captions/inference.rst'
workdir: config["workdir"]

# Every rule should have shadow: copy-minimal.

rule all:
    input:
        data = "output/simulations/simtest-2936549/data.tar",
        parameters = "output/simulations/simtest-2936549/parameters.tsv"

rule combine_simulation_tasks:
    output:
        data = "output/simulations/{sim}/data.tar",
        parameters = "output/simulations/{sim}/parameters.tsv"
    conda: "envs/simulate.yaml"
    script: "scripts/"
