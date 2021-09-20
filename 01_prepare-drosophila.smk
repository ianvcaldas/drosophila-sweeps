
configfile: '01_config.yaml'
report: 'captions/prepare-drosophila.rst'


def wait_for_empirical_npy(wildcards):
    # This line halts workflow execution until the empirical_windows_012 rule runs:
    checkpoints.empirical_windows_commands.get(**wildcards)
    # That rule creates some files; we can extract their filenames now.
    window_names = glob_wildcards('output/empirical-windows/commands/{name}.command').name
    # And generate the desired input for the `all` rule from that list of names.
    return expand('output/empirical-windows/npy/{name}.npy', name=window_names)


rule all:
    input:
        expand(
            'output/resistant-lines/{locus}-resistant-lines.txt',
            locus=['ace', 'chkov', 'cyp']
        ),
        'output/dgrp2/imputed-counts.frq.count',
        'output/dgrp2/imputed-heteroz.sites.pi',
        expand(
            'output/selection-scan/{chrom}-features{stats}.tsv',
            chrom=['2R', '3R'],
            stats=['', '-stats']
        ),
        # Get set of empirical windows from another rule and add their .npy as inputs to
        # this rule.
        wait_for_empirical_npy


rule selection_scan_features:
    input:
        genotypes = 'output/selection-scan/{chrom}.012'
    output:
        features = 'output/selection-scan/{chrom}-features.tsv',
        stats = 'output/selection-scan/{chrom}-features-stats.tsv',
        ms = temp('output/selection-scan/{chrom}.ms')
    conda: 'envs/simulate.yaml'
    benchmark: 'benchmarks/selection-scan-features/{chrom}.tsv'
    notebook: 'notebooks/prepare-drosophila/selection-scan.py.ipynb'


rule genotypes_2R:
    input: 'output/dgrp2/imputed.vcf.gz'
    output: temp('output/selection-scan/2R.012')
    conda: 'envs/simulate.yaml'
    shell:
        "vcftools --gzvcf {input} "
        "--chr 2R "
        "--mac 1 "
        "--012 "
        "--out \"output/selection-scan/2R\" "


rule genotypes_3R:
    input:
        vcf = 'output/dgrp2/imputed.vcf.gz',
        ace_pos = 'output/resistant-lines/ace-genotype.012.pos'
    output: temp('output/selection-scan/3R.012')
    conda: 'envs/simulate.yaml'
    shell:
        "vcftools --gzvcf {input.vcf} "
        "--chr 3R "
        "--mac 1 "
        "--exclude-positions \"{input.ace_pos}\" "
        "--012 "
        "--out \"output/selection-scan/3R\" "


rule empirical_window_features:
    input:
        genotypes = 'output/empirical-windows/genotypes/{window}.012',
        normalization_stats = config['stats-file-location']
    output:
        npy = 'output/empirical-windows/npy/{window}.npy',
        ms = temp('output/empirical-windows/ms/{window}.ms'),
        features = temp('output/empirical-windows/features/{window}.tsv'),
        stats = 'output/empirical-windows/features/{window}-stats.tsv',
    params:
        outdir = 'output/empirical-windows/npy'
    conda: 'envs/simulate.yaml'
    benchmark: 'benchmarks/empirical-window-features/{window}.tsv'
    notebook: 'notebooks/prepare-drosophila/empirical-window-features.py.ipynb'


rule extract_empirical_windows:
    input: 'output/empirical-windows/commands/{window}.command'
    output: temp('output/empirical-windows/genotypes/{window}.012')
    conda: 'envs/simulate.yaml'
    shell:
        "bash {input}"


checkpoint empirical_windows_commands:
    input:
        'output/dgrp2/imputed.vcf.gz',
        'output/resistant-lines/ace-genotype.012.pos',
        # Notebook as input re-runs the rule if we decide to change empirical windows.
        'notebooks/prepare-drosophila/empirical-windows-commands.py.ipynb'
    output: directory('output/empirical-windows/commands')
    conda: 'envs/simulate.yaml'
    benchmark: 'benchmarks/empirical-windows-commands.tsv'
    log:
        notebook = 'logs/prepare-drosophila/empirical-windows-commands.py.ipynb'
    notebook: 'notebooks/prepare-drosophila/empirical-windows-commands.py.ipynb'


rule resistant_lines:
    input:
        expand(
            'output/resistant-lines/{locus}-genotype.{ending}',
            locus=['ace', 'cyp'],
            ending=['012', '012.indv', '012.pos']
        ),
        jiggins = config['jiggins-data']
    output: 
        expand(
            'output/resistant-lines/{locus}-resistant-lines.txt',
            locus=['ace', 'chkov', 'cyp']
        ),
        'output/resistant-lines/dgrp1-dgrp2-lines-comparison.txt'
    conda: 'envs/simulate.yaml'
    benchmark: 'benchmarks/resistant-lines.tsv'
    log:
        notebook = 'logs/prepare-drosophila/resistant-lines.py.ipynb'
    notebook: 'notebooks/prepare-drosophila/resistant-lines.py.ipynb'


rule dgrp_stats:
    input: 'output/dgrp2/imputed.vcf.gz'
    output:
        'output/dgrp2/imputed-counts.frq.count',
        'output/dgrp2/imputed-heteroz.sites.pi'
    conda: 'envs/simulate.yaml'
    shell:
        "vcftools --gzvcf {input} "
        "--counts "
        "--out \"output/dgrp2/imputed-counts\" ;"
        "vcftools --gzvcf {input} "
        "--site-pi "
        "--out \"output/dgrp2/imputed-heteroz\" ;"


rule beagle_impute:
    input: 'output/dgrp2/dgrp2-snps-maxmissing15percent.recode.vcf.gz'
    output: 'output/dgrp2/imputed.vcf.gz'
    conda: 'envs/simulate.yaml'
    threads: 8
    params:
        beagle = config['beagle-location']
    shell:
        "java -Xmx16g -jar {params.beagle} gt=\"{input}\" out=\"output/dgrp2/imputed\" nthreads={threads}"


rule cyp_resistant_lines:
    input: config['dgrp2-location']
    output:
        multiext(
            'output/resistant-lines/cyp-genotype',
            '.012', '.012.indv', '.012.pos'
        )
    conda: 'envs/simulate.yaml'
    shell:
        "vcftools --gzvcf {input} "
        "--012 "
        "--out \"output/resistant-lines/cyp-genotype\" "
        "--snp 2R_8072884_INS "


rule ace_resistant_lines:
    input: config['dgrp2-location']
    output:
        multiext(
            'output/resistant-lines/ace-genotype',
            '.012', '.012.indv', '.012.pos'
        )
    conda: 'envs/simulate.yaml'
    shell:
        "vcftools --gzvcf {input} "
        "--012 "
        "--out \"output/resistant-lines/ace-genotype\" "
        "--snp 3R_9069054_SNP "
        "--snp 3R_9069408_SNP "
        "--snp 3R_9069721_SNP "


rule filter_dgrp:
    input: config['dgrp2-location']
    output: 'output/dgrp2/dgrp2-snps-maxmissing15percent.recode.vcf.gz'
    conda: 'envs/simulate.yaml'
    shell:
        "vcftools --gzvcf {input} "
        "--mac 1 "
        "--max-alleles 2 "
        "--min-alleles 2 "
        "--max-missing 0.85 "
        "--recode "
        "--remove-filtered-all "
        "--remove-indels "
        "--out \"output/dgrp2/dgrp2-snps-maxmissing15percent\" ;"
        "gzip output/dgrp2/dgrp2-snps-maxmissing15percent.recode.vcf ;"
