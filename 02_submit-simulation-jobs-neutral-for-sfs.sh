CONFIG=02_configs/main-fixedsweeps-neutral.yaml
NAME="neutral"
NUM_JOBS=100

sbatch --array=1-${NUM_JOBS} -J ${NAME} 02_slurm-simulate.sh --configfile ${CONFIG}
