# Silence Tensorflow warnings
from gc import callbacks
import os

from yaml import parse
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

import argparse
from multiprocessing.sharedctypes import Value
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
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
parser.add_argument("--reg_rate", type=float, default=0.0,
                    help="Regularization rate", required=False)
parser.add_argument("--reg_type", type=str, default="l2",
                    help="Type of regularization. One of l1, l2 (default), l1_l2", required=False)
parser.add_argument("--layers", type=int, default=3,
                    help="Number of hidden layers", required=True)
parser.add_argument("--units", type=int, nargs='+',
                    help="Number of units")
parser.add_argument("--dataset", type=str,
                    help="Id of the dataset used to train the network", required=True)
parser.add_argument("--run_name", type=str, help="Name of the Weights and Biases run", required=False)
parser.add_argument("--WB-log", action=argparse.BooleanOptionalAction,
                    help="Whether to log to Weights and Biases (--WB-log) or not (--no-WB-log).")

args = parser.parse_args()

# Check that we have the correct number of units for the number of layers
if args.units is not None and len(args.units) != args.layers:
    raise ValueError("Number of units must be equal to number of layers")

class cell_labeler():
    def __init__(self, epochs: int, batch_size: int, lr: float, dropout: float, layers: int, 
                units: List[int], dataset: str, run_name: None) -> None:
        """
        Initialize the cell_labeler class.

        param epochs: int - Number of epochs to train the model
        param batch_size: int - Batch size for training
        param lr: float - Learning rate
        param dropout: float - Dropout rate
        param layers: int - Number of hidden layers
        param units: List[int] - Number of units in each hidden layer
        param dataset: str - Id of the dataset used to train the model
        param run_name: str - Name of the run. If None (default), a random name is generated.
        """

        if units is None:
            units = [2**(layers + 5 - i) for i in range(layers)]

        self.epochs = epochs
        self.batch_size = batch_size
        self.lr = lr
        self.dropout = dropout
        self.layers = layers
        self.units = units
        self.reg_rate = args.reg_rate
        self.reg_type = args.reg_type
        self.dataset = dataset

        self.model = None
        self.x_train, self.x_test, self.y_train, self.y_test = None, None, None, None

        datasets = pd.read_csv("datasets.csv")

        if self.dataset not in datasets["study_id"].values:
            raise ValueError(
                f"Dataset {self.dataset} not found.\nThese are the available datasets: {datasets['study_id'].unique()}")

        if args.WB_log is True:
            wandb.init(project="scRNA-seq-classification",
                    config={"batch_size": self.batch_size,
                            "epochs": self.epochs,
                            "lr": self.lr,
                            "dropout": self.dropout,
                            "layers": self.layers,
                            "units": self.units,
                            "reg_rate": self.reg_rate,
                            "reg_type": self.reg_type,
                            "dataset": self.dataset})

            if run_name is not None:
                wandb.run.name = run_name 
                wandb.run.save()

            wandb.log({"accuracy": 0.0, 
                    "loss": 0.0})

        filename_expr = f"expr_matrices/{self.dataset}.csv.gz"
        filename_labels = f"expr_matrices/{self.dataset}_clusters.csv"

        self.expr = pd.read_csv(filename_expr, index_col=0)        
        self.expr = self.expr.T
        self.clusters = pd.read_csv(filename_labels, index_col=0)
        # Note the +1 to take into account the "other" class
        self.n_clusters = self.clusters["Cluster"].unique().shape[0] + 1
        self.clusters_names = [f"Cluster {i}" for i in range(self.n_clusters-1)] + ["Other"]

    def __build_model(self) -> None:
        """
        Build a multi-class MLP classifier.
        """

        # Now build our MLP
        self.model = keras.Sequential()
        self.model.add(keras.layers.InputLayer(
            input_shape=(self.expr.shape[1],)))

        if self.reg_type == "l1":
            reg = keras.regularizers.l1(self.reg_rate)
        elif self.reg_type == "l2":
            reg = keras.regularizers.l2(self.reg_rate)
        elif self.reg_type == "l1_l2":
            reg = keras.regularizers.l1_l2(self.reg_rate) 

        for i in range(self.layers):
            self.model.add(keras.layers.Dense(
                self.units[i], activation="relu",
                
                name=f"dense_{i}"))
            self.model.add(keras.layers.Dropout(
                self.dropout, name=f"dropout_{i}"))

        # Output layer
        self.model.add(keras.layers.Dense(self.n_clusters,
                       kernel_regularizer=reg,
                       activation="softmax", name="output"))

        self.model.compile(optimizer=keras.optimizers.Adam(lr=self.lr),
                    loss="categorical_crossentropy", metrics=["accuracy"])

    def __prepare_training_data(self, n_augment: int, perc_shuffle_genes_aug: float) -> None:
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

        # One-hot encode labels
        self.y_train = keras.utils.to_categorical(
            self.y_train, num_classes=self.n_clusters)
        self.y_test = keras.utils.to_categorical(
            self.y_test, num_classes=self.n_clusters)

        x_train_aug, y_train_aug = self.__get_augmented_samples(
            n=n_augment, perc_shuffle_genes=perc_shuffle_genes_aug)
        x_train_aug = self.scaler.transform(x_train_aug)

        self.x_train = np.vstack([self.x_train, x_train_aug])
        self.y_train = np.vstack([self.y_train, y_train_aug])


    def __get_augmented_samples(self, n: int = 1, perc_shuffle_genes: float = 0.01) -> Tuple[np.array, np.array]:
        """
        Creates augmented samples from the training data.

        param: n (int) - the number of augmented samples to create. Default is 1.
        param: perc_shuffle_genes (float) - the percentage of genes to shuffle (per cell). Default is 0.01.
        return: x_train_aug, y_train_aug (np.array, np.array) - the augmented training data
        """
        # Get a list of cluster labels
        y_train_aug = np.random.choice(np.unique(np.argmax(self.y_train, axis=1)), size=n)
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

        self.__prepare_training_data(n_augment=500, perc_shuffle_genes_aug=0.01)

        for train_index, test_index in kf.split(self.x_train):
            self.__build_model()

            # TODO add early stopping and checkpointing
            model_callbacks = []

            if (args.WB_log):
                model_callbacks.append(WandbCallback(save_model=False))

            self.model.fit(self.x_train[train_index], self.y_train[train_index],
                           epochs=self.epochs, batch_size=self.batch_size,
                           validation_data=(self.x_train[test_index], self.y_train[test_index]),
                           callbacks=[model_callbacks])
            # Get validation accuracy
            val_acc.append(self.model.evaluate(self.x_train[test_index], self.y_train[test_index])[1])
        
        print(f"Validation accuracy: {np.mean(val_acc)}")
        # Push to W&B
        wandb.log({"Validation accuracy": np.mean(val_acc)})

    def evaluate_model(self) -> Tuple[float, float]:
        """
        Evaluate the model.

        return: accuracy (float) - the accuracy of the model
        """

        # Evaluate the model
        test_loss, test_accuracy = self.model.evaluate(self.x_test, self.y_test)
        y_pred = self.model.predict(self.x_test)

        # Push to W&B
        if (args.WB_log):
            wandb.log({"Test loss": test_loss, "Test accuracy": test_accuracy,
                    "conf_mat" : wandb.plot.confusion_matrix(probs=None,
                        y_true=np.argmax(self.y_test, axis = 1), preds=np.argmax(y_pred, axis = 1),
                        class_names=self.clusters_names)})
        
        return (test_loss, test_accuracy)

labeller = cell_labeler(epochs=args.epochs,
                        batch_size=args.batch_size,
                        lr=args.lr,
                        dropout=args.dropout,
                        layers=args.layers,
                        units=args.units,
                        dataset=args.dataset,
                        run_name=f"confusion_mtx_{args.dataset}")

labeller.train_model()
labeller.evaluate_model()

if (args.WB_log):
    wandb.finish()