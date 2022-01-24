
rule all:
    input: "fig/overfitting-learning-curves.pdf"

rule overfitting_learning_curves:
    output: "fig/overfitting-learning-curves.pdf"
    conda: "envs/plotting.yaml"
    notebook: "notebooks/plotting/overfitting-learning-curves.r.ipynb"
