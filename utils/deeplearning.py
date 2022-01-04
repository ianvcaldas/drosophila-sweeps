import tarfile
from io import BytesIO

import tqdm
import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype

import torch
import torch.nn as nn


class SweepsDataset(torch.utils.data.Dataset):
    def __init__(
        self, data_tar, df, target_column, is_validation=False, feature_subset=None
    ):
        self.df = df
        self.target_column = target_column
        self.task = self.get_task()[0]
        self.labels = self.get_labels()
        self.feature_subset = feature_subset
        self.data = self.untar_data(data_tar)
        self.is_validation = is_validation

    def untar_data(self, tar):
        result = dict()
        with tarfile.open(tar) as data:
            uuids = self.df.uuid
            for uuid in tqdm.notebook.tqdm(uuids):
                name = uuid + ".npy"
                array_file = BytesIO()
                array_file.write(data.extractfile(name).read())
                array_file.seek(0)
                np_array = np.load(array_file).transpose([2, 0, 1])
                if self.feature_subset is not None:
                    np_array = np_array[self.feature_subset, :, :]
                tensor = torch.Tensor(np_array).div_(255)
                result[uuid] = tensor
        return result

    def get_task(self):
        """Returns the type of inference problem and number of labels."""
        labels_col = self.df[self.target_column]
        if pd.api.types.is_numeric_dtype(labels_col):
            return ("regression", 1)
        else:
            return ("classification", labels_col.nunique())

    def get_labels(self):
        if self.task == "classification":
            return sorted(self.df[self.target_column].unique())
        elif self.task == "regression":
            return self.target_column

    def __getitem__(self, ix):
        simulation_uid = self.df.iloc[ix].uuid
        item = self.data[simulation_uid]
        label = self.df.iloc[ix][self.target_column]
        if self.task == "classification":
            y = self.labels.index(label)
            y_tensor = torch.LongTensor([y]).squeeze()
        elif self.task == "regression":
            y = label
            y_tensor = torch.Tensor([y]).squeeze()
        return item, y_tensor

    def __len__(self):
        return self.df.shape[0]


class SimpleCNN2Layer(nn.Module):
    def __init__(self, input_dim, output_dim, in_channels=1):
        super().__init__()

        self.conv1 = nn.Conv2d(
            in_channels=in_channels,
            out_channels=128,
            stride=1,
            kernel_size=2,
            padding=1,
        )
        self.relu1 = nn.ReLU()
        self.maxpool1 = nn.MaxPool2d(kernel_size=2)

        self.conv2 = nn.Conv2d(
            in_channels=128, out_channels=64, stride=1, kernel_size=2, padding=1
        )
        self.relu2 = nn.ReLU()
        self.maxpool2 = nn.MaxPool2d(kernel_size=2)

        self.fc = nn.Linear(in_features=self.fc_dim(input_dim), out_features=output_dim)

    def fc_dim(self, indim):
        """Manually calculate the required dimension of the last fully connected layer of the network. This uses the parameters of layers set in the network description, so if the network architecture changes at all, this need to change as well!"""
        after_first_conv = self.outdim(
            indim, padding=1, dilation=1, kernel_size=2, stride=1
        )
        after_first_max = self.outdim(
            after_first_conv, padding=0, dilation=1, kernel_size=2, stride=2
        )
        after_sec_conv = self.outdim(
            after_first_max, padding=1, dilation=1, kernel_size=2, stride=1
        )
        after_sec_max = self.outdim(
            after_sec_conv, padding=0, dilation=1, kernel_size=2, stride=2
        )
        current_channels = 64
        fc_size = current_channels * after_sec_max ** 2
        return int(fc_size)

    def outdim(self, indim, padding, dilation, kernel_size, stride):
        """
        Calculates output dimension of a convolutional or maxpool layer. The formula is
        the same for both kinds of layer.
        References:
        https://pytorch.org/docs/stable/generated/torch.nn.Conv2d.html
        https://pytorch.org/docs/stable/generated/torch.nn.MaxPool2d.html
        """
        x = indim + 2 * padding - dilation * (kernel_size - 1) - 1
        return x / stride + 1

    def forward(self, x):
        out = self.maxpool1(self.relu1(self.conv1(x)))
        out = self.maxpool2(self.relu2(self.conv2(out)))
        out = out.view(out.size(0), -1)
        out = self.fc(out)
        return out


def get_inferences(model, loader, dataset_name):
    data = {"training": (loader.train_ds, 0), "validation": (loader.valid_ds, 1)}
    dataset_object, dataset_idx = data[dataset_name]
    predictions, targets = model.get_preds(ds_idx=dataset_idx)
    if dataset_object.task == "classification":
        predicted_ix = predictions.numpy().argmax(axis=1)
        true_ix = targets.numpy()
        labels = dataset_object.labels
        true_label = [labels[i] for i in true_ix]
        predicted_label = [labels[i] for i in predicted_ix]
        result = (
            pd.DataFrame.from_records(predictions.numpy(), columns=labels)
            .assign(
                uuid=dataset_object.df.uuid,
                true_ix=true_ix,
                predicted_ix=predicted_ix,
                true_label=true_label,
                predicted_label=predicted_label,
            )
            .set_index("uuid")
        )
    elif dataset_object.task == "regression":
        label = dataset_object.labels
        result = (
            pd.DataFrame.from_records(
                predictions.numpy(), columns=["predicted_" + label]
            )
            .assign(uuid=dataset_object.df.uuid, true=targets.numpy(),)
            .set_index("uuid")
            .rename({"true": "true_" + label}, axis="columns")
        )
    return result
