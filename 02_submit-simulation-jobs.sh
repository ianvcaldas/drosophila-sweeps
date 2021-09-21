
NUM_SLURM_JOBS=2

for file in $(ls 02_configs/*.yaml)
do
    NAME=$(basename ${file} .yaml)
    COMMAND="sbatch --array=1-${NUM_SLURM_JOBS} -J ${NAME} 02_slurm-simulate.sh --configfile ${file}"
    echo ${COMMAND}
done


