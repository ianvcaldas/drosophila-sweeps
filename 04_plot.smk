
configfile: '03_config.yaml'

rule all:
    input:
        expand("fig/overfitting-learning-curves_{training}.pdf", training=config["training_ids"])

rule overfitting_learning_curves:
    input:
        expand(
            "output/model-fitting/{target}_{{training}}_overfitting.tsv",
            target=config["inference_targets"]
        )
    output: "fig/overfitting-learning-curves_{training}.pdf"
    conda: "envs/plotting.yaml"
    notebook: "notebooks/plotting/overfitting-learning-curves.r.ipynb"
