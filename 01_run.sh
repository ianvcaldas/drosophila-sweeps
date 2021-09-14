
snakemake -c4 --use-conda --snakefile 01_prepare-drosophila.smk "$@"

# To re-run all the rules that depend on a script or notebook after it's changed, run
# snakemake --forcerun [rule with script or notebook].
