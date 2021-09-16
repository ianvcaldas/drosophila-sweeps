
CONFIGFILE=02_configs/testing-sgv.yaml

snakemake -c4 --use-conda --snakefile 02_simulate.smk --configfile ${CONFIGFILE} "$@" --config slim=bin/slim3.6 normalization_stats=resources/normalization-stats.tsv

# To re-run all the rules that depend on a script or notebook after it's changed, run
# snakemake --forcerun [rule with script or notebook].