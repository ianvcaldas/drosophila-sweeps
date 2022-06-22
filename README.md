
Sel. sweep parameter inference
==============================

This repository contains the entirety of the code used to infer sweep parameters in *Drosophila melanogaster* using machine learning. This code is part of an upcoming preprint.

In the paper, we propose a method for inferring selective sweep parameters based on simulations. The method's main steps are:

1. Simulate selective sweeps datasets
2. Calculate population genetic summary statistics for the simulated datasets and for empirical data
3. Train machine learning models on the simulated datasets; apply them to the empirical dataset

To reproduce this method and analysis, any implementation of those 3 steps is reasonable and our specific choices are only one of many possibilities. In this read me, first I'll describe how we approached those 3 important steps of the method and where to find the code for them so that our approach can serve as a basis or inspiration for further uses of the method. Then, I'll describe how we put the whole analysis pipeline together for the paper.


Reproducing the method
----------------------

### Simulation

To simulate selective sweeps, we use a combination of coalescence and forward-time simulation in SLiM. We pick simulation parameters from random distributions specified in YAML files; see `resources/simulation-parameters` for examples. Default parameters are in `utils/project-parameters.py`. Parameter distributions are then instantiated into concrete values and injected into SLiM templates. The templates are located in `resources/slim-templates`.

The simulation logic resides in the Snakefile `02_simulate.smk`. (See below for examples on how to run it.) The Snakefile takes as input a YAML configuration with simulation parameters and a number of simulations to run. It will then perform a series of steps for each simulation (each step is implemented as a script in `scripts/simulations`).

1. Instantiate a set of parameters from the specified random distributions.
2. Perform a coalescent burn-in.
3. Create a SLiM script from a template.
4. Run SLiM with the coalescent burn-in as a starting state and output a tree sequence.
5. Drop mutations on the tree sequence, then output the resulting genotypes in ms format.
6. Calculate population genetic summary statistics (see below).

At the end, information from all simulations is concatenated into a table of parameters, while the calculated features are compressed together in a `.tar` archive.


### Feature calculation

Our feature calculation method, described in the paper, is implemented in `utils/feature-calculation.py`. The calculation takes in a ms-format file as input and produces a table of feature values. We use our own implementation of population genetic summary statistics in `utils/popgen_summary_statistics.py`.

We then normalize the table of features and output the result as a three-dimensional matrix in binary `npy` format. The entire process is done for simulations in `scripts/simulations/features-and-normalization.py` and for the empirical data in `notebooks/prepare-drosophila/empirical-window-features.py.ipynb`. Parameter normalization values used in the paper are in `resources/normalization-stats.tsv`.

For feature calculation, we need to specify the data dimension, smallest subwindow size, and total locus size. We specify them in `utils/project-parameters.py`.


### Machine learning inference

Our implementation of deep learing uses the fastai v2 framework. All necessary tooling to run the framework is in `utils/deeplearning.py`. That module contains the neural network architecture, which is a CNN adapted to take in a number of input channels larger than 3, and the `SweepsDataset` class that loads data from the compressed `.tar` archives produced by earlier steps. The machine learning hyperparameters, such as the number of epochs to train for and the batch size for training, are in `config.yaml`.

We train CNNs to infer many different parameters, or "targets of inference", from the data. For each one of those, data was balanced such that classes or sweep modes were equally represented. The scripts are in `scripts/inference`.

After those steps are satisfied, the code that invokes fastai to perform model fitting and inference is in `notebooks/inference/fit-neural-network.py.ipynb`.


Workflow details
----------------

The paper is implemented as 4 separate Snakemake workflows:

1. Prepare the empirical data
2. Simulate selective sweep datasets
3. Train models on simulated data and perform inferences using them
4. Generate tables and plots of results

The workflows are set up to load the conda environments in `envs` automatically.

The workflows were developed with Snakemake 6.8, so that's the recommended version to use to run the code in the paper.


### Repository structure

- `envs`: Conda environment specifications.
- `notebooks`: Jupyter notebooks with the bulk of the analysis and plotting code.
- `resources`: Auxiliary files, including the SLiM script templates, parameters used for simulations, and bounds used for summary statistic normalization.
- `scripts`: Python scripts with code for preparing training and validation datasets for machine learning, as well as the main logic for the simulation workflow.
- `utils`: Analysis and plotting code that's reused across notebooks and scripts.


### Step 1: Prepare empirical data

The first workflow prepares the empirical data for the analyses downstream. It does filtering and imputation, then calculates metrics of the empirical data: the number of SNPs, pairwise heterozygosity, and SFS. These empirical statistics are used to determine parameters of our sweep simulations and to guarantee that the simulations have similar statistic distributions to the empirical data. This workflow will also perform the feature calculation step on the empirical data, producing input files that can be input to the machine learning models downstream.

The workflow can be run with:

    snakemake --use-conda --snakefile 01_prepare-drosophila.smk -c<X>

Where `<X>` is the number of cores.

	
### Step 2: Simulate

The second workflow performs simulations. As described above, one run of it will perform simulations with one set of simulation parameters specified in a YAML file. For practical reasons, simulations for the paper were split into many sets of parameter files, datasets, computers, and computational processes. We have one YAML file per sweep mode (hard, RNM, or SGV) per simulated dataset. There are many simulated datasets; some are used for training machine learning models, others aren't used for training but have the models applied to them for robustness analyses. Others never interact with machine learning but are used for plotting. All simulation parameter YAMLs can be found in `resources/simulation-parameters`. All simulations share a common set of default parameters, specified in `utils/project-parameters.py`.

In practice, we found that running the workflow in multiple replicate copies for each YAML file in a remote cluster worked best for generating simulated datasets. The scripts used for that (using Cornell's SLURM cluster architecture) are `slurm_submit.sh` and `slurm_simulate.sh`. Under this scenario, each time the workflow runs, it only needs to run a handful of simulations. Thousands replicates of results are then aggregated later, as the first step in the next workflow.

As an example of running the simulation workflow, the following snippet will run locally 5 simulations for every training YAML:

```
CONFIGFILES=$(ls resources/simulation-parameters/training/*.yaml)
NUM_SIMS=5

for file in ${CONFIGFILES}
do
	snakemake --use-conda --snakefile 02_simulate.smk --configfile ${file} --config slim=bin/slim3.7 normalization_stats=resources/normalization-stats.tsv simulations=${NUM_SIMS} use_subdirectory=true
done
```


### Step 3: Machine learning inference

The third workflow will first bring together all the output from the many simulation tasks run in Step 2 into coherent folders for each dataset. This workflow will then perform a split of successful simulations between a training a validation dataset, then balance the dataset for each target of inference, as described above.

For each target of inference and each simulated dataset tagged for training, this workflow trains a series of replicate machine learning models. These trained models are then applied to every simulated dataset tagged for testing, as well as the empirical data, to generate lists of inferences.

Run the whole workflow with:

    snakemake --use-conda --snakefile 03_inference.smk -c<X>

Where `<X>`is the number of cores.


### Step 4: Get results

The fourth workflow takes the series of inferences performed by the models in the previous steps and converts them into results tables and plots. Machine learning performance metrics, such as accuracy, RMSE, and ROC curves, are calculated at this step. Other auxiliary plots in the paper are also produced.

Run the whole workflow with:

    snakemake --use-conda --snakefile 04_results.smk -c<X>

Where `<X>`is the number of cores.
