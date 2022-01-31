
configfile: '03_config.yaml'

rule all:
    input:
        expand("fig/overfitting-learning-curves_{training}.pdf",
               training=config["training_ids"]),
        "fig/subwindow-diagram.pdf"

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
