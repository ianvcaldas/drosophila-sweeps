
CONFIGFILES=$(ls 02_configs/test/*.yaml)

for file in ${CONFIGFILES}
do
    snakemake -c4 --use-conda --snakefile 02_simulate.smk --configfile ${file} "$@" --config slim=bin/slim3.6 normalization_stats=resources/normalization-stats.tsv simulations=3 use_subdirectory=true
done


# To re-run all the rules that depend on a script or notebook after it's changed, run
# snakemake --forcerun [rule with script or notebook].
