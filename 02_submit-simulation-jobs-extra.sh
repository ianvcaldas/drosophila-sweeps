NUM_EXTRA_TRAINING_RNM=1100
NUM_EXTRA_TRAINING_SGV=550

NUM_EXTRA_ROBUSTNESS_RNM=200
NUM_EXTRA_ROBUSTNESS_SGV=100

for file in $(ls 02_configs/training/*sgv.yaml)
do
    NAME=$(basename ${file} .yaml)
    sbatch --array=1-${NUM_EXTRA_TRAINING_SGV} -J ${NAME} 02_slurm-simulate.sh --configfile ${file}
done

for file in $(ls 02_configs/training/*rnm.yaml)
do
    NAME=$(basename ${file} .yaml)
    sbatch --array=1-${NUM_EXTRA_TRAINING_RNM} -J ${NAME} 02_slurm-simulate.sh --configfile ${file}
done


for file in $(ls 02_configs/robustness/*sgv.yaml)
do
    NAME=$(basename ${file} .yaml)
    sbatch --array=1-${NUM_EXTRA_ROBUSTNESS_SGV} -J ${NAME} 02_slurm-simulate.sh --configfile ${file}
done

for file in $(ls 02_configs/robustness/*rnm.yaml)
do
    NAME=$(basename ${file} .yaml)
    sbatch --array=1-${NUM_EXTRA_ROBUSTNESS_RNM} -J ${NAME} 02_slurm-simulate.sh --configfile ${file}
done

