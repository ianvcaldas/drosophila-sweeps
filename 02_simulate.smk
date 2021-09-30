
# UTILITY #
# ======= #

from random import choice

if 'slim' not in config:
    raise Exception("Need to provide a path to SLiM executable. Please invoke snakemake with --config slim=PATH_TO_SLIM.")

if 'normalization_stats' not in config:
    raise Exception("Need to provide a path to normalization stats file. Please invoke snakemake with --config normalization_stats=PATH.")

if 'use_subdirectory' not in config:
    config['use_subdirectory'] = False

def get_outdir():
    outdir = Path("simulations")
    if config['use_subdirectory']:
        vowels = 'aeiou'
        consonants = 'bcdfghjklmnprstvz'
        subdir = choice(consonants) + choice(vowels) + choice(consonants) + choice(vowels)
        return outdir / subdir
    else:
        return outdir

OUTDIR_ID = get_outdir()
OUTDIR = Path('output')/OUTDIR_ID
SIM_IDS = [
    str(i + 1).zfill(len(str(config['simulations'])))
    for i in range(config['simulations'])
]

def wait_for_slim_checkpoint(wildcards):
    # This line halts workflow execution until the slim rule runs:
    checkpoints.slim_has_run.get(**wildcards)
    # That rule creates some files; we can extract their filenames now.
    successful_sim_ids = glob_wildcards(OUTDIR/'{sim_id}_slim.trees').sim_id
    return successful_sim_ids


def successful_npy(wildcards):
    successful_sim_ids = wait_for_slim_checkpoint(wildcards)
    return expand(OUTDIR/'npy/{sim_id}.npy', sim_id=successful_sim_ids)

def successful_lognpy(wildcards):
    successful_sim_ids = wait_for_slim_checkpoint(wildcards)
    return expand(OUTDIR/'npy-log-scale/{sim_id}.npy', sim_id=successful_sim_ids)

def successful_metrics(wildcards):
    successful_sim_ids = wait_for_slim_checkpoint(wildcards)
    desired_inputs = expand(
        OUTDIR/"{sim_id}_{param_file}",
        sim_id=successful_sim_ids,
        param_file=["slim-metrics.txt", "msprime-metrics.txt"]
    ) + successful_npy(wildcards)
    return desired_inputs

def successful_features(wildcards):
    successful_sim_ids = wait_for_slim_checkpoint(wildcards)
    return expand(OUTDIR/'features/{sim_id}_features.tsv', sim_id=successful_sim_ids)
    


# RULES #
# ===== #


rule all:
    input:
        OUTDIR/"parameters.tsv",
        OUTDIR/"features.tar.gz",
        OUTDIR/"data.tar",
        OUTDIR/"logdata.tar"

rule logtar_npy:
    input:
        npy_files = successful_lognpy,
        params_file = OUTDIR/"parameters.tsv"
    output:
        tar_file = OUTDIR/"logdata.tar"
    params:
        npy_folder = OUTDIR/"npy-log-scale",
    conda: "envs/simulate.yaml"
    script: "scripts/simulations/tar_npy.py"

rule tar_npy:
    input:
        npy_files = successful_npy,
        params_file = OUTDIR/"parameters.tsv"
    output:
        tar_file = OUTDIR/"data.tar"
    params:
        npy_folder = OUTDIR/"npy",
    conda: "envs/simulate.yaml"
    script: "scripts/simulations/tar_npy.py"


rule compress_features:
    input: successful_features
    output: OUTDIR/"features.tar.gz"
    params:
        features_folder = OUTDIR/"features"
    shell:
        "mkdir -p {params.features_folder} ;"
        "tar -czf {output} {params.features_folder} ;"


rule aggregate_parameters:
    input:
        simulation_params = expand(
            OUTDIR/"{sim_id}_simulation-params.json",
            sim_id=SIM_IDS
        ),
        successful = successful_metrics
    output:
        params_output = OUTDIR/"parameters.tsv"
    params:
        benchmarks_folder = "benchmarks/" + str(OUTDIR_ID),
        simulations_folder = OUTDIR
    conda: "envs/simulate.yaml"
    script: "scripts/simulations/aggregate-parameters.py"


rule calculate_normalize_features:
    input:
        ms = OUTDIR/"{sim_id}_genotypes.ms",
        normalization_stats = config['normalization_stats'],
        params_file = OUTDIR/"{sim_id}_simulation-params.json"
    output:
        npy = OUTDIR/"npy/{sim_id}.npy",
        npy_log = OUTDIR/"npy-log-scale/{sim_id}.npy",
        features = OUTDIR/"features/{sim_id}_features.tsv",
        stats = OUTDIR/"features/{sim_id}_feature-stats.tsv"
    conda: 'envs/simulate.yaml'
    benchmark: 'benchmarks/' + str(OUTDIR_ID) + '/{sim_id}_features-and-normalization.tsv'
    script: "scripts/simulations/features-and-normalization.py"


rule drop_mutations:
    input:
        slim_output = OUTDIR/"{sim_id}_slim.trees",
        params_file = OUTDIR/"{sim_id}_simulation-params.json"
    output:
        trees_file = OUTDIR/"{sim_id}_mutation-dropped.trees",
        metrics_file = OUTDIR/"{sim_id}_msprime-metrics.txt",
        ms_file = OUTDIR/"{sim_id}_genotypes.ms"
    benchmark: 'benchmarks/' + str(OUTDIR_ID) + '/{sim_id}_drop-mutations.tsv'
    conda: "envs/simulate.yaml"
    script: "scripts/simulations/drop-mutations.py"


checkpoint slim_has_run:
    input:
        expand(OUTDIR/".{sim_id}_slim.done", sim_id=SIM_IDS)
    output:
        touch(OUTDIR/".slim_all.done")


rule slim:
    input:
        OUTDIR/"{sim_id}_script.slim"
    output:
        touch(OUTDIR/".{sim_id}_slim.done")
    params:
        slim = config['slim']
    log: 'logs/' + str(OUTDIR_ID) + '/{sim_id}_slim.log'
    benchmark: 'benchmarks/' + str(OUTDIR_ID) + '/{sim_id}_slim.tsv'
    shell:
        "{params.slim} {input} &> {log}"


rule slim_script:
    input:
        params_file = OUTDIR/"{sim_id}_simulation-params.json",
        burnin = OUTDIR/"{sim_id}_burnin.trees",
        slim_templates = expand(
            'resources/slim-templates/{regime}.slim',
            regime=['neutral', 'hard', 'rnm', 'sgv', 'size-change']
        )
    output:
        slim_script = OUTDIR/"{sim_id}_script.slim"
    params:
        outdir = OUTDIR
    conda: "envs/simulate.yaml"
    script: "scripts/simulations/make-slim-script.py"


rule burnin:
    input:
        params_file = OUTDIR/"{sim_id}_simulation-params.json"
    output:
        trees = OUTDIR/"{sim_id}_burnin.trees"
    conda: "envs/simulate.yaml"
    benchmark: 'benchmarks/' + str(OUTDIR_ID) + '/{sim_id}_burnin.tsv'
    script: "scripts/simulations/simulation-burnin.py"


rule instantiate_parameters:
    output:
        params_file = OUTDIR/"{sim_id}_simulation-params.json"
    conda: "envs/simulate.yaml"
    script: "scripts/simulations/instantiate-simulation-parameters.py"
