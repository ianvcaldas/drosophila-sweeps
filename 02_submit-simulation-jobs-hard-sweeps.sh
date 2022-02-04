NUM_HARD_JOBS=20

for file in $(ls 02_configs/hard-sweeps-with-known-s/*.yaml)
do
    NAME=$(basename ${file} .yaml)
    sbatch --array=1-${NUM_HARD_JOBS} -J ${NAME} 02_slurm-simulate.sh --configfile ${file}
done
