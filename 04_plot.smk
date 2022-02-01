
configfile: '03_config.yaml'

def get_feature_subsets(num_features):
    """Feature subsets are coded as a string, where each character is "0" or "1" for
    presence or absence of a feature."""
    no_features = ["0"]*num_features
    all_features = ["1"]*num_features
    results = []
    for ix in range(num_features):
        only_this_one = no_features[:]
        only_this_one[ix] = "1"
        without_this_one = all_features[:]
        without_this_one[ix] = "0"
        results.extend(["".join(only_this_one), "".join(without_this_one)])
    return results

rule all:
    input:
        expand("fig/overfitting-learning-curves_{training}.pdf",
               training=config["training_ids"]),
        "fig/subwindow-diagram.pdf",
        "fig/fixed-sweeps-validation.pdf",
        "fig/fixed-sweeps-validation-selection-brackets.pdf",
        "fig/fixed-sweeps-validation-secondary-models.pdf",
        "fig/partial-sweeps-validation.pdf"

rule partial_sweeps_validation:
    input:
        parameters = "output/simulation-data-processed/train-valid-split/main-partialsweeps_validation.tsv",
        selstrength = "output/inferences-training/log-sel-strength_main-partialsweeps_validation.tsv",
        sweepmode = "output/inferences-training/sweep-mode_main-partialsweeps_validation.tsv",
        sweepmode_roc =
        "output/metrics/sweep-mode_main-partialsweeps_roc-curve.tsv"
    output: "fig/partial-sweeps-validation.pdf"
    conda: "envs/plotting.yaml"
    notebook: "notebooks/plotting/partialsweeps-validation.r.ipynb"

rule fixed_sweeps_validation_secondary_models:
    input:
        parameters = "output/simulation-data-processed/train-valid-split/main-fixedsweeps_validation.tsv",
        rnm_vs_sgv = "output/inferences-training/rnm-vs-sgv_main-fixedsweeps_validation.tsv",
        rnm_vs_sgv_roc = "output/metrics/rnm-vs-sgv_main-fixedsweeps_roc-curve.tsv",
        hard_vs_soft = "output/inferences-training/hard-vs-soft_main-fixedsweeps_validation.tsv",
        hard_vs_soft_roc = "output/metrics/hard-vs-soft_main-fixedsweeps_roc-curve.tsv",
        feature_analysis = "output/metrics/main-fixedsweeps_metrics-with-features.tsv",
        feature_analysis_code = "output/metrics/main-fixedsweeps_feature-analysis-grid.tsv"
    output: "fig/fixed-sweeps-validation-secondary-models.pdf"
    conda: "envs/plotting.yaml"
    notebook: "notebooks/plotting/fixedsweeps-validation-secondary-models.r.ipynb"

rule fixed_sweeps_validation_selection_brackets:
    input:
        parameters = "output/simulation-data-processed/train-valid-split/main-fixedsweeps_validation.tsv",
        selstrength = "output/inferences-training/log-sel-strength_main-fixedsweeps_validation.tsv",
        sweepmode = "output/inferences-training/sweep-mode_main-fixedsweeps_validation.tsv"
    output: "fig/fixed-sweeps-validation-selection-brackets.pdf"
    conda: "envs/plotting.yaml"
    notebook: "notebooks/plotting/selection-brackets-validation.r.ipynb"

rule fixed_sweeps_validation:
    input:
        parameters = "output/simulation-data-processed/train-valid-split/main-fixedsweeps_validation.tsv",
        selstrength = "output/inferences-training/log-sel-strength_main-fixedsweeps_validation.tsv",
        sweepmode = "output/inferences-training/sweep-mode_main-fixedsweeps_validation.tsv",
        sweepmode_roc =
        "output/metrics/sweep-mode_main-fixedsweeps_roc-curve.tsv",
        feature_analysis = "output/metrics/main-fixedsweeps_metrics-with-features.tsv",
        feature_analysis_code = "output/metrics/main-fixedsweeps_feature-analysis-grid.tsv"
    output: "fig/fixed-sweeps-validation.pdf"
    conda: "envs/plotting.yaml"
    notebook: "notebooks/plotting/fixedsweeps-validation.r.ipynb"

rule prepare_metrics:
    input:
        feature_subsets=expand(
            "output/inferences-feature-analysis/{target}_{{training}}_features-{features}_validation.tsv",
            target=config["inference_targets"],
            features=get_feature_subsets(num_features=7)
        ),
        baseline=expand(
            "output/inferences-training/{target}_{{training}}_validation.tsv",
            target=config["inference_targets"]
        )
    output:
        feature_analysis = "output/metrics/{training}_metrics-with-features.tsv",
        feature_analysis_code = "output/metrics/{training}_feature-analysis-grid.tsv"
    conda: "envs/ml.yaml"
    notebook: "notebooks/plotting/prepare-metrics-with-features.py.ipynb"

rule prepare_roc_curve:
    input:
        data="output/inferences-training/{target}_{training}_validation.tsv",
        labels="output/trained-models/{target}_{training}_labels.txt"
    output:
        curve="output/metrics/{target}_{training}_roc-curve.tsv",
        auc="output/metrics/{target}_{training}_roc-curve-auc.tsv"
    conda: "envs/ml.yaml"
    notebook: "notebooks/plotting/prepare-roc-curve.py.ipynb"


rule subwindow_diagram:
    input: "output/metrics/subwindow-sizes.tsv"
    output: "fig/subwindow-diagram.pdf"
    conda: "envs/plotting.yaml"
    notebook: "notebooks/plotting/subwindow-diagram.r.ipynb"

rule prepare_subdwindow_data:
    output: "output/metrics/subwindow-sizes.tsv"
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
