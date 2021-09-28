NUM_TRAINING_JOBS=1100
NUM_ROBUSTNESS_JOBS=200

for file in $(ls 02_configs/training/*.yaml)
do
    NAME=$(basename ${file} .yaml)
    sbatch --array=1-${NUM_TRAINING_JOBS} -J ${NAME} 02_slurm-simulate.sh --configfile ${file}
done

for file in $(ls 02_configs/robustness/*.yaml)
do
    NAME=$(basename ${file} .yaml)
    sbatch --array=1-${NUM_ROBUSTNESS_JOBS} -J ${NAME} 02_slurm-simulate.sh --configfile ${file}
done


