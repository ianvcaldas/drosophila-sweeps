
snakemake -c4 --use-conda --snakefile 03_inference.smk "output/inferences-empirical/sweep-mode_main-fixedsweeps_empirical.tsv" "$@"
