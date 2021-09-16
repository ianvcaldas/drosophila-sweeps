
NUM_SLURM_JOBS=2

sbatch --array=1-${NUM_SLURM_JOBS} -J simtest 02_slurm-simulate.sh --configfile 02_configs/testing.yaml
