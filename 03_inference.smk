
configfile: '03_config.yaml'
report: 'captions/inference.rst'

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
        data_reports = expand(
            "output/simulation-data-processed/info/{sim_id}_data-report.txt",
            sim_id=(config["training_ids"] + config["testing_ids"])
        ),
        testing_inferences = expand(
            "output/inferences-testing/{target}_{training}_{testing}.tsv",
            target=config["inference_targets"],
            training=config["training_ids"],
            testing=config["testing_ids"]
        ),
        partial_sweeps_testing_inferences = expand(
            "output/inferences-partial/{target}_main-fixedsweeps_main-partialsweeps.tsv",
            target=config["inference_targets"]
        ),
        training_inferences = expand(
            "output/inferences-training/{target}_{training}_{testing}.tsv",
            target=config["inference_targets"],
            training=config["training_ids"],
            testing=["training", "validation"]
        ),
        empirical_inferences = expand(
            "output/inferences-empirical/{target}_{training}_empirical.tsv",
            target=config["inference_targets"],
            training=config["training_ids"]
        ),
        feature_analysis = expand(
            "output/inferences-feature-analysis/{target}_{training}_features-{features}_{testing}.tsv",
            target=config["inference_targets"],
            # We only do feature analysis for fixed sweeps training dataset
            training=config["training_ids"][0],
            features=get_feature_subsets(num_features=7),
            testing=["training", "validation"]
        ),
        overfitting_reports = expand(
            "output/model-fitting/{target}_{training}_overfitting.tsv",
            target=config["inference_targets"],
            training=config["training_ids"]
        ),
        gradient_boost_inferences = expand(
            "output/inferences-gradientboost/{target}_{training}_{testing}.tsv",
            target=config["inference_targets"],
            training=config["training_ids"],
            testing=["training", "validation"]
        ),
        empirical_inference_replicates = expand(
            "output/inferences-empirical-replicates/{target}_{training}_empirical_replicate-{k}.tsv",
            target=config["inference_targets"],
            training=config["training_ids"],
            k=range(config["num_empirical_replicates"])
        ),
        learning_curve_per_sample_size = expand(
            "output/model-fitting-learning-curve/{target}_{training}_fit_{k}.tsv",
            target=config["inference_targets"],
            training=config["training_ids"],
            k=range(config["num_overfitting_replicates"])
        )


rule learning_curve:
    input:
        training = "output/simulation-data-processed/balanced/{target}_{training}_training.tsv",
        validation = "output/simulation-data-processed/balanced/{target}_{training}_validation.tsv",
        data = "output/simulation-data/{training}/data.tar",
        logdata  = "output/simulation-data/{training}/logdata.tar"
    output:
        fit_report = "output/model-fitting-learning-curve/{target}_{training}_fit_{k}.tsv"
    params:
        use_log_data = False,
        epochs = config["epochs_for_model_training"]
    conda: "envs/ml.yaml"
    notebook: "notebooks/inference/fit-neural-network-learning-curve.py.ipynb"


rule apply_replicate_model_to_empirical_data:
    input:
        fit_model = "output/trained-models-replicates/{target}_{training}_replicate-{k}.pth",
        model_object = "output/trained-models-replicates/{target}_{training}_replicate-{k}.pkl",
        model_labels = "output/trained-models-replicates/{target}_{training}_labels_replicate-{k}.txt",
        data = "output/empirical-windows/data.tar",
        logdata = "output/empirical-windows/logdata.tar"
    output:
        inferences = "output/inferences-empirical-replicates/{target}_{training}_empirical_replicate-{k}.tsv"
    params:
        application_type = "empirical"
    conda: "envs/ml.yaml"
    notebook: "notebooks/inference/apply-model.py.ipynb"


rule fit_model_replicate:
    input:
        training = "output/simulation-data-processed/balanced/{target}_{training}_training.tsv",
        validation = "output/simulation-data-processed/balanced/{target}_{training}_validation.tsv",
        data = "output/simulation-data/{training}/data.tar",
        logdata  = "output/simulation-data/{training}/logdata.tar"
    output:
        fit_model = "output/trained-models-replicates/{target}_{training}_replicate-{k}.pth",
        model_object = "output/trained-models-replicates/{target}_{training}_replicate-{k}.pkl",
        model_labels = "output/trained-models-replicates/{target}_{training}_labels_replicate-{k}.txt",
        training_inferences = "output/inferences-training-replicates/{target}_{training}_training_replicate-{k}.tsv",
        validation_inferences = "output/inferences-training-replicates/{target}_{training}_validation_replicate-{k}.tsv",
        fit_report = "output/model-fitting-replicates/{target}_{training}_fit_replicate-{k}.tsv"
    params:
        save_model = True,
        save_inferences = True,
        use_log_data = False,
        epochs = config["epochs_for_model_training"]
    conda: "envs/ml.yaml"
    notebook: "notebooks/inference/fit-neural-network.py.ipynb"


rule fit_gradient_boost:
    input:
        training = "output/simulation-data-processed/balanced/{target}_{training}_training.tsv",
        validation = "output/simulation-data-processed/balanced/{target}_{training}_validation.tsv",
        data = "output/simulation-data/{training}/data.tar",
        logdata  = "output/simulation-data/{training}/logdata.tar"
    output:
        training_inferences = "output/inferences-gradientboost/{target}_{training}_training.tsv",
        validation_inferences = "output/inferences-gradientboost/{target}_{training}_validation.tsv"
    params:
        num_gb_estimators = config["num_gradient_boosting_estimators"],
        use_log_data = False
    conda: "envs/ml.yaml"
    benchmark: "benchmarks/inference/fit-gradient-boost_{target}_{training}.tsv"
    log: "logs/inference/fit-gradient-boost_{target}_{training}.py.ipynb"
    notebook: "notebooks/inference/fit-gradient-boost.py.ipynb"


rule apply_model_to_empirical_data:
    input:
        fit_model = "output/trained-models/{target}_{training}.pth",
        model_object = "output/trained-models/{target}_{training}.pkl",
        model_labels = "output/trained-models/{target}_{training}_labels.txt",
        data = "output/empirical-windows/data.tar",
        logdata = "output/empirical-windows/logdata.tar"
    output:
        inferences = "output/inferences-empirical/{target}_{training}_empirical.tsv"
    params:
        application_type = "empirical"
    conda: "envs/ml.yaml"
    benchmark: "benchmarks/inference/apply-model-empirical_{target}_{training}.tsv"
    log: "logs/inference/apply-model-empirical_{target}_{training}.py.ipynb"
    notebook: "notebooks/inference/apply-model.py.ipynb"


rule apply_fixed_sweeps_model_to_partial_sweeps:
    input:
        fit_model = "output/trained-models/{target}_main-fixedsweeps.pth",
        model_object = "output/trained-models/{target}_main-fixedsweeps.pkl",
        model_labels = "output/trained-models/{target}_main-fixedsweeps_labels.txt",
        data = "output/simulation-data/main-partialsweeps/data.tar",
        sim_parameters = "output/simulation-data-processed/train-valid-split/main-partialsweeps_validation.tsv",
        logdata = "output/simulation-data/main-partialsweeps/logdata.tar",
    output:
        inferences = "output/inferences-partial/{target}_main-fixedsweeps_main-partialsweeps.tsv",
    params:
        application_type = "testing"
    conda: "envs/ml.yaml"
    benchmark: "benchmarks/apply-model-testing_{target}_main-fixedsweeps_main-partialsweeps.tsv"
    log: "logs/inference/apply-model-testing_{target}_main-fixedsweeps_main-partialsweeps.py.ipynb"
    notebook: "notebooks/inference/apply-model.py.ipynb"



rule apply_model_to_testing_data:
    input:
        fit_model = "output/trained-models/{target}_{training}.pth",
        model_object = "output/trained-models/{target}_{training}.pkl",
        model_labels = "output/trained-models/{target}_{training}_labels.txt",
        data = "output/simulation-data/{testing}/data.tar",
        sim_parameters = "output/simulation-data-processed/parameters/{testing}_parameters-clean.tsv",
        logdata  = "output/simulation-data/{testing}/logdata.tar"
    output:
        inferences = "output/inferences-testing/{target}_{training}_{testing}.tsv",
    params:
        application_type = "testing"
    conda: "envs/ml.yaml"
    benchmark: "benchmarks/apply-model-testing_{target}_{training}_{testing}.tsv"
    log: "logs/inference/apply-model-testing_{target}_{training}_{testing}.py.ipynb"
    notebook: "notebooks/inference/apply-model.py.ipynb"


rule aggregate_overfitting_replicates:
    input:
        expand(
            "output/model-fitting/{{target}}_{{training}}_replicate-{k}.tsv",
            k=range(config["num_overfitting_replicates"])
        )
    output:
        aggregated = "output/model-fitting/{target}_{training}_overfitting.tsv"
    conda: "envs/simulate.yaml"
    benchmark: "benchmarks/inference/aggregate-overfitting_{target}_{training}.tsv"
    script: "scripts/inference/aggregate-overfitting-replicates.py"


rule fit_model_feature_subset:
    input:
        training = "output/simulation-data-processed/balanced/{target}_{training}_training.tsv",
        validation = "output/simulation-data-processed/balanced/{target}_{training}_validation.tsv",
        data = "output/simulation-data/{training}/data.tar",
        logdata  = "output/simulation-data/{training}/logdata.tar"
    output:
        training_inferences = "output/inferences-feature-analysis/{target}_{training}_features-{features}_training.tsv",
        validation_inferences = "output/inferences-feature-analysis/{target}_{training}_features-{features}_validation.tsv",
        fit_report = "output/model-fitting/{target}_{training}_features-{features}_fit.tsv"
    params:
        save_model = False,
        save_inferences = True,
        use_log_data = False,
        epochs = config["epochs_for_model_training"]
    conda: "envs/ml.yaml"
    benchmark: "benchmarks/inference/fit-neural-network_{target}_{training}_features-{features}.tsv"
    log: "logs/inference/fit-neural-network_{target}_{training}_features-{features}.py.ipynb"
    notebook: "notebooks/inference/fit-neural-network.py.ipynb"


rule fit_model:
    input:
        training = "output/simulation-data-processed/balanced/{target}_{training}_training.tsv",
        validation = "output/simulation-data-processed/balanced/{target}_{training}_validation.tsv",
        data = "output/simulation-data/{training}/data.tar",
        logdata  = "output/simulation-data/{training}/logdata.tar"
    output:
        fit_model = "output/trained-models/{target}_{training}.pth",
        model_object = "output/trained-models/{target}_{training}.pkl",
        model_labels = "output/trained-models/{target}_{training}_labels.txt",
        training_inferences = "output/inferences-training/{target}_{training}_training.tsv",
        validation_inferences = "output/inferences-training/{target}_{training}_validation.tsv",
        fit_report = "output/model-fitting/{target}_{training}_fit.tsv"
    params:
        save_model = True,
        save_inferences = True,
        use_log_data = False,
        epochs = config["epochs_for_model_training"]
    conda: "envs/ml.yaml"
    benchmark: "benchmarks/inference/fit-neural-network_{target}_{training}.tsv"
    log: "logs/inference/fit-neural-network_{target}_{training}.py.ipynb"
    notebook: "notebooks/inference/fit-neural-network.py.ipynb"


rule balance_training_data:
    input:
        training = "output/simulation-data-processed/train-valid-split/{training}_training.tsv",
        validation = "output/simulation-data-processed/train-valid-split/{training}_validation.tsv"
    output:
        balanced_training = "output/simulation-data-processed/balanced/{target}_{training}_training.tsv",
        balanced_validation = "output/simulation-data-processed/balanced/{target}_{training}_validation.tsv"
    params:
        random_seed = 13
    conda: "envs/simulate.yaml"
    benchmark: "benchmarks/inference/balance-training-data_{target}_{training}.tsv"
    script: "scripts/inference/balance-training-data.py"

    
rule train_validation_split:
    input:
        sim_params = "output/simulation-data-processed/parameters/{training}_parameters-clean.tsv"
    output:
        training = "output/simulation-data-processed/train-valid-split/{training}_training.tsv",
        validation = "output/simulation-data-processed/train-valid-split/{training}_validation.tsv"
    params:
        random_seed = 13
    conda: "envs/ml.yaml"
    benchmark: "benchmarks/inference/train-validation-split_{training}.tsv"
    script: "scripts/inference/train-validation-split.py"


rule overfitting_simple_fit:
    input:
        training = "output/simulation-data-processed/balanced-overfitting/{target}_{training}_replicate-{k}_training.tsv",
        validation = "output/simulation-data-processed/balanced-overfitting/{target}_{training}_replicate-{k}_validation.tsv",
        data = "output/simulation-data/{training}/data.tar",
        logdata  = "output/simulation-data/{training}/logdata.tar"
    output:
        fit_report = "output/model-fitting/{target}_{training}_replicate-{k}.tsv"
    params:
        save_model = False,
        save_inferences = False,
        use_log_data = False,
        epochs = config["epochs_for_overfitting_analyses"]
    conda: "envs/ml.yaml"
    benchmark: "benchmarks/inference/overfitting-simple-fit_{target}_{training}_replicate-{k}.tsv"
    log: "logs/inference/overfitting-simple-fit_{target}_{training}_replicate-{k}.py.ipynb"
    notebook: "notebooks/inference/fit-neural-network.py.ipynb"


rule overfitting_balance_training_data:
    input:
        training = "output/simulation-data-processed/train-valid-split-overfitting/{training}_replicate-{k}_training.tsv",
        validation = "output/simulation-data-processed/train-valid-split-overfitting/{training}_replicate-{k}_validation.tsv"
    output:
        balanced_training = "output/simulation-data-processed/balanced-overfitting/{target}_{training}_replicate-{k}_training.tsv",
        balanced_validation = "output/simulation-data-processed/balanced-overfitting/{target}_{training}_replicate-{k}_validation.tsv"
    conda: "envs/simulate.yaml"
    benchmark: "benchmarks/inference/overfitting-balance-training-data_{target}_{training}_replicate-{k}.tsv"
    script: "scripts/inference/balance-training-data.py"


rule overfitting_train_validation_split:
    input:
        sim_params = "output/simulation-data-processed/parameters/{training}_parameters-clean.tsv"
    output:
        training = "output/simulation-data-processed/train-valid-split-overfitting/{training}_replicate-{k}_training.tsv",
        validation = "output/simulation-data-processed/train-valid-split-overfitting/{training}_replicate-{k}_validation.tsv"
    conda: "envs/ml.yaml"
    benchmark: "benchmarks/inference/overfitting-train-validation-split_{training}_replicate-{k}.tsv"
    script: "scripts/inference/train-validation-split.py"


rule clean_datasets:
    input:
        sim_params = "output/simulation-data/{sim_id}/parameters.tsv"
    output:
        cleaned_parameters = "output/simulation-data-processed/parameters/{sim_id}_parameters-clean.tsv",
        sweep_mode_report = "output/simulation-data-processed/info/{sim_id}_sweep-modes.txt",
        successful_report = "output/simulation-data-processed/info/{sim_id}_data-report.txt",
        failed_report = "output/simulation-data-processed/info/{sim_id}_failed-simulations-report.txt"
    params:
        random_seed = 13
    conda: "envs/simulate.yaml"
    benchmark: "benchmarks/inference/clean-dataset_{sim_id}.tsv"
    script: "scripts/inference/clean-dataset.py"


rule combine_simulation_tasks:
    output:
        data = "output/simulation-data/{sim}/data.tar",
        parameters = "output/simulation-data/{sim}/parameters.tsv",
        info = "output/simulation-data/{sim}/info.txt",
        features = "output/simulation-data/{sim}/features.tar.gz",
        logdata  = "output/simulation-data/{sim}/logdata.tar"
    conda: "envs/simulate.yaml"
    benchmark: "benchmarks/inference/combine-simulations_{sim}.tsv"
    log: "logs/inference/combine-simulations_{sim}.py.ipynb"
    notebook: "notebooks/inference/combine-simulations.py.ipynb"
