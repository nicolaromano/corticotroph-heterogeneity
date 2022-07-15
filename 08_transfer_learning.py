# Silence Tensorflow warnings
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

import argparse
from multiprocessing.sharedctypes import Value
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
# import KFold from scikit-learn
from sklearn.model_selection import KFold
from wandb.keras import WandbCallback
import wandb
from tensorflow import keras
from typing import List, Tuple
from tqdm import tqdm
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--epochs", type=int, default=50,
                    help="Number of epochs to train", required=True)
parser.add_argument("--batch_size", type=int, default=32,
                    help="Batch size", required=True)
parser.add_argument("--lr", type=float, default=0.001,
                    help="Learning rate", required=True)
parser.add_argument("--dropout", type=float, default=0.5,
                    help="Dropout rate (for all layers)", required=True)
parser.add_argument("--layers", type=int, default=3,
                    help="Number of hidden layers", required=True)
parser.add_argument("--units", type=int, nargs='+',
                    help="Number of units", required=True)
parser.add_argument("--ref_dataset", type=int,
                    default=0, help="Reference dataset")

args = parser.parse_args()

if args.units is None:
    args.units = [2**(args.layers + 5 - i) for i in range(args.layers)]

# Check that we have the correct number of units for the number of layers
if len(args.units) != args.layers:
    raise ValueError("Number of units must be equal to number of layers")


class cell_labeler():
    def __init__(self, epochs: int, batch_size: int, lr: float, dropout: float, layers: int, units: List[int],
                 ref_dataset: int):
        self.epochs = epochs
        self.batch_size = batch_size
        self.lr = lr
        self.dropout = dropout
        self.layers = layers
        self.units = units
        self.ref_dataset = ref_dataset

        datasets = pd.read_csv("datasets.csv")

        if self.ref_dataset not in datasets["study_id"]:
            raise ValueError(
                f"Dataset not found.\nThese are the available datasets: {datasets['study_id'].unique()}")

        study_id = datasets["study_id"][self.ref_dataset]

        wandb.init(project="scRNA-seq-classification",
                config={"batch_size": self.batch_size,
                        "epochs": self.epochs,
                        "lr": self.lr,
                        "dropout": self.dropout,
                        "layers": self.layers,
                        "units": self.units,
                        "ref_dataset": datasets["study_id"][self.ref_dataset]})

        wandb.log({"accuracy": 0.0, "loss": 0.0})

        filename_expr = f"expr_matrices/{study_id}.csv.gz"
        filename_labels = f"expr_matrices/{study_id}_clusters.csv"

        self.expr = pd.read_csv(filename_expr, index_col=0)
        self.clusters = pd.read_csv(filename_labels)

    def build_model(self) -> None:
        """
        Build a multi-class MLP classifier.
        """

        # Now build our MLP
        self.model = keras.Sequential()
        self.model.add(keras.layers.InputLayer(
            input_shape=(self.expr.shape[1],)))

        for i in range(self.layers):
            self.model.add(keras.layers.Dense(
                self.units[i], activation="relu", name=f"dense_{i}"))
            self.model.add(keras.layers.Dropout(
                self.dropout, name=f"dropout_{i}"))

        # Output layer
        # Note the +1 to take into account the "other" class
        self.n_clusters = self.clusters.unique().shape[0] + 1
        self.model.add(keras.layers.Dense(self.n_clusters,
                       activation="softmax", name="output"))

        self.model.compile(optimizer=keras.optimizers.Adam(lr=self.lr),
                    loss="categorical_crossentropy", metrics=["accuracy"])

    def prepare_training_data(self, n_augment: int, perc_shuffle_genes_aug: float) -> None:
        """
        Prepare the training data for the MLP. Generates a training and test set.

        param n_augment: Number of augmented data points to generate
        param perc_shuffle_genes_aug: Percentage of genes to shuffle in the augmented data
        """

        # We split into training (90%) and test (10%) data
        self.x_train, self.x_test, self.y_train, self.y_test = train_test_split(
            self.expr, self.clusters, test_size=0.1, random_state=42)

        # Scale the data with MinMaxScaler
        self.scaler = MinMaxScaler()
        self.x_train = self.scaler.fit_transform(self.x_train)
        self.x_test = self.scaler.transform(self.x_test)

        x_train_aug, y_train_aug = self.get_augmented_samples(
            n_augment=500, perc_shuffle_genes_aug=0.01)
        x_train_aug = self.scaler.transform(x_train_aug)

        self.x_train = np.vstack([self.x_train, x_train_aug])
        self.y_train = np.vstack([self.y_train, y_train_aug])

        # One-hot encode labels
        self.y_train = keras.utils.to_categorical(
            self.y_train, num_classes=self.n_clusters)
        self.y_test = keras.utils.to_categorical(
            self.y_test, num_classes=self.n_clusters)

    def get_augmented_samples(self, n: int = 1, perc_shuffle_genes: float = 0.01) -> Tuple[np.array, np.array]:
        """
        Creates augmented samples from the training data.

        param: n (int) - the number of augmented samples to create. Default is 1.
        param: perc_shuffle_genes (float) - the percentage of genes to shuffle (per cell). Default is 0.01.
        return: x_train_aug, y_train_aug (np.array, np.array) - the augmented training data
        """

        # Get a list of cluster labels
        y_train_aug = np.random.choice(
            np.unique(np.argmax(self.y_train, axis=1)), size=n)

        # Shuffle the data on a per-cluster/per-gene basis

        # Split x_train by cluster
        x_by_cluster = [self.x_train[np.where(
            np.argmax(self.y_train, axis=1) == i)] for i in range(self.n_clusters - 1)]

        for x in x_by_cluster:
            genes = np.random.choice(x.shape[1], int(
                perc_shuffle_genes * x.shape[1]), replace=False)

            # Randomly swaps the values of the cells gene-wise
            for g in genes:
                x[:, g] = np.random.permutation(x[:, g])

        # Pick shuffled data corresponding to the cluster labels
        x_train_aug = np.vstack([x_by_cluster[i][np.random.randint(
            0, x_by_cluster[i].shape[0]), :] for i in y_train_aug])

        # One-hot encode labels
        y_train_aug = keras.utils.to_categorical(
            y_train_aug, num_classes=self.n_clusters)

        return (x_train_aug, y_train_aug)

    def train_model(self) -> None:
        """
        Train the MLP model.
        """

        # Do 5-fold cross-validation
        kf = KFold(n_splits=5, shuffle=True)        
        
        val_acc = []
        for train_index, test_index in kf.split(self.x_train):
            self.build_model()
            self.prepare_training_data(n_augment=500, perc_shuffle_genes_aug=0.01)

            self.model.fit(self.x_train[train_index], self.y_train[train_index],
                           epochs=self.epochs, batch_size=self.batch_size,
                           validation_data=(self.x_train[test_index], self.y_train[test_index]))
            # Get validation accuracy
            val_acc.append(self.model.evaluate(self.x_train[test_index], self.y_train[test_index])[1])
        
        print(f"Validation accuracy: {np.mean(val_acc)}")

    def evaluate_model(self) -> float:
        """
        Evaluate the model.

        return: accuracy (float) - the accuracy of the model
        """

        # Evaluate the model, return accuracy
        return (self.model.evaluate(self.x_test, self.y_test)[2])

labeller = cell_labeler(epochs=args.epochs,
                        batch_size=args.batch_size,
                        lr=args.lr,
                        dropout=args.dropout,
                        layers=args.layers,
                        units=args.units,
                        ref_dataset=args.ref_dataset)

labeller.train_model()
labeller.evaluate_model()

wandb.finish()