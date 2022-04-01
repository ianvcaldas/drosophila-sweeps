
# SLURM Submissions
###################


# Neutral simulations, for finding the simulated neutral SFS

CONFIG=resources/simulation-parameters/main-fixedsweeps-neutral.yaml
NAME="neutral"
NUM_JOBS=100

sbatch --array=1-${NUM_JOBS} -J ${NAME} slurm_simulate.sh --configfile ${CONFIG}


# Training data and data for checking robustness

NUM_TRAINING_JOBS=1100
NUM_ROBUSTNESS_JOBS=200

for file in $(ls resources/simulation-parameters/training/*.yaml)
do
    NAME=$(basename ${file} .yaml)
    sbatch --array=1-${NUM_TRAINING_JOBS} -J ${NAME} slurm_simulate.sh --configfile ${file}
done

for file in $(ls resources/simulation-parameters/robustness/*.yaml)
do
    NAME=$(basename ${file} .yaml)
    sbatch --array=1-${NUM_ROBUSTNESS_JOBS} -J ${NAME} slurm_simulate.sh --configfile ${file}
done


# Extra training and robustness simulations

NUM_EXTRA_TRAINING_RNM=1100
NUM_EXTRA_TRAINING_SGV=550

NUM_EXTRA_ROBUSTNESS_RNM=200
NUM_EXTRA_ROBUSTNESS_SGV=100

for file in $(ls resources/simulation-parameters/training/*sgv.yaml)
do
    NAME=$(basename ${file} .yaml)
    sbatch --array=1-${NUM_EXTRA_TRAINING_SGV} -J ${NAME} slurm_simulate.sh --configfile ${file}
done

for file in $(ls resources/simulation-parameters/training/*rnm.yaml)
do
    NAME=$(basename ${file} .yaml)
    sbatch --array=1-${NUM_EXTRA_TRAINING_RNM} -J ${NAME} slurm_simulate.sh --configfile ${file}
done


for file in $(ls resources/simulation-parameters/robustness/*sgv.yaml)
do
    NAME=$(basename ${file} .yaml)
    sbatch --array=1-${NUM_EXTRA_ROBUSTNESS_SGV} -J ${NAME} slurm_simulate.sh --configfile ${file}
done

for file in $(ls resources/simulation-parameters/robustness/*rnm.yaml)
do
    NAME=$(basename ${file} .yaml)
    sbatch --array=1-${NUM_EXTRA_ROBUSTNESS_RNM} -J ${NAME} slurm_simulate.sh --configfile ${file}
done


# Hard sweeps with known s, for sanity checks

NUM_HARD_JOBS=20

for file in $(ls resources/simulation-parameters/hard-sweeps-with-known-s/*.yaml)
do
    NAME=$(basename ${file} .yaml)
    sbatch --array=1-${NUM_HARD_JOBS} -J ${NAME} slurm_simulate.sh --configfile ${file}
done
