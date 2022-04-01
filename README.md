
Drosophila
==========

The project is split between 4 Snakemake workflows, for each of the main steps of the analysis:

1. Filter and prepare the empirical data
2. Simulate selective sweeps
3. Train models on simulated data and perform inferences on empirical data
4. Generate plots of results


### Step 1: Filter and prepare empirical data

    snakemake --use-conda --snakefile 01_prepare-drosophila.smk -c4
	
	
### Step 2: Simulate selective sweeps

To simulate locally:

    CONFIGFILES=$(ls resources/simulation-parameters/training/*.yaml)
	NUM_SIMS=5
    
    for file in ${CONFIGFILES}
    do
        snakemake --use-conda --snakefile 02_simulate.smk --configfile ${file} --config slim=bin/slim3.7 normalization_stats=resources/normalization-stats.tsv simulations=${NUM_SIMS} use_subdirectory=true
    done

`slurm_simulate.sh` was the script used to perform simulations on a SLURM cluster. `slurm_submit.sh` was the script used to submit the jobs to the cluster.


### Step 3: Machine learning inference

    snakemake --use-conda --snakefile 03_inference.smk -c8


### Step 4: Generate plots

    snakemake --use-conda --snakefile 04_plot.smk -c4
