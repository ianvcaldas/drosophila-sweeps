
configfile: 'config.yaml'


rule all:
    input: "fig/clr.pdf"


rule plot_clr:
    input:
        expand(
            "output/clr/result/clr_{chrom}.tsv",
            chrom=["2L", "2R", "3L", "3R"]
        )
    output: "fig/clr.pdf"
    conda: "envs/plotting.yaml"
    notebook: "notebooks/plotting/clr.r.ipynb"


rule clr:
    input:
        data = "output/clr/input/clr_{chrom}.tsv",
        sfs = "output/clr/input/clr_genomewide-sfs.tsv"
    output: "output/clr/result/clr_{chrom}.tsv"
    log: "output/clr/result/clr_{chrom}.log"
    params:
        grid = 10000
    shell:
        "bin/sweed/SweeD -name {wildcards.chrom} -grid {params.grid} -folded -ploidy 1 -isfs {input.sfs} -input {input.data} &> {log} ;"
        "rm SweeD_Info.{wildcards.chrom} ;"
        "mv SweeD_Report.{wildcards.chrom} {output} ;"


rule clr_input:
    input:
        vcf = "output/dgrp2/imputed.vcf.gz"
    output:
        genomewide_sfs = "output/clr/input/clr_genomewide-sfs.tsv",
        input_per_chrom = expand(
            "output/clr/input/clr_{chrom}.tsv",
            chrom=["2L", "2R", "3L", "3R"]
        )
    conda: "envs/simulate.yaml"
    notebook: "notebooks/prepare-drosophila/clr-input.py.ipynb"
