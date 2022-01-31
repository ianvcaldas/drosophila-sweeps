
configfile: '03_config.yaml'

rule all:
    input:
        expand("fig/overfitting-learning-curves_{training}.pdf",
               training=config["training_ids"]),
        "fig/subwindow-diagram.pdf",
        "fig/fixed-sweeps-validation.pdf"

rule fixed_sweeps_validation:
    input:
        parameters = "output/simulation-data-processed/balanced/log-sel-strength_main-fixedsweeps_validation.tsv",
        selstrength = "output/inferences-training/log-sel-strength_main-fixedsweeps_validation.tsv",
        sweepmode = "output/inferences-training/sweep-mode_main-fixedsweeps_validation.tsv",
        sweepmode_roc = "output/extra-information/sweep-mode_main-fixedsweeps_roc-curve.tsv"
    output: "fig/fixed-sweeps-validation.pdf"
    conda: "envs/plotting.yaml"
    notebook: "notebooks/plotting/main-validation.r.ipynb"

rule prepare_roc_curve:
    input:
        data="output/inferences-training/{target}_{training}_validation.tsv",
        labels="output/trained-models/{target}_{training}_labels.txt"
    output: "output/extra-information/{target}_{training}_roc-curve.tsv"
    conda: "envs/ml.yaml"
    notebook: "notebooks/plotting/prepare-roc-curve.py.ipynb"


rule subwindow_diagram:
    input: "output/extra-information/subwindow-sizes.tsv"
    output: "fig/subwindow-diagram.pdf"
    conda: "envs/plotting.yaml"
    notebook: "notebooks/plotting/subwindow-diagram.r.ipynb"

rule prepare_subdwindow_data:
    output: "output/extra-information/subwindow-sizes.tsv"
    conda: "envs/simulate.yaml"
    notebook: "notebooks/plotting/prepare-subwindow-data.py.ipynb"


rule overfitting_learning_curves:
    input:
        expand(
            "output/model-fitting/{target}_{{training}}_overfitting.tsv",
            target=config["inference_targets"]
        )
    output: "fig/overfitting-learning-curves_{training}.pdf"
    conda: "envs/plotting.yaml"
    notebook: "notebooks/plotting/overfitting-learning-curves.r.ipynb"
