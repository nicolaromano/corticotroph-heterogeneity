# Silence Tensorflow warnings
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

import pandas as pd
import numpy as np
from tqdm import tqdm
from typing import List, Tuple
from tensorflow import keras
import wandb
from wandb.keras import WandbCallback
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
from multiprocessing.sharedctypes import Value
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--epochs", type=int, default=50,
                    help="Number of epochs to train")
parser.add_argument("--batch_size", type=int, default=32, help="Batch size")
parser.add_argument("--lr", type=float, default=0.001, help="Learning rate")
parser.add_argument("--dropout", type=float, default=0.5,
                    help="Dropout rate (for all layers)")
parser.add_argument("--layers", type=int, default=3,
                    help="Number of hidden layers")
parser.add_argument("--units", type=int, nargs='+',
                    default=[512, 256, 32], help="Number of units")
parser.add_argument("--ref_dataset", type=int,
                    default=0, help="Reference dataset")

args = parser.parse_args()

# Check that we have the correct number of units for the number of layers
if len(args.units) != args.layers:
    raise ValueError("Number of units must be equal to number of layers")

datasets = pd.read_csv("datasets.csv")

if args.ref_dataset not in datasets["study_id"]:
    raise ValueError(
        f"Dataset not found.\nThese are the available datasets: {datasets['study_id'].unique()}")

study_id = datasets["study_id"][args.ref_dataset]
print(f"Using dataset {study_id} as reference")

wandb.init(project="scRNA-seq-classification",
           config={"batch_size": args.batch_size,
                   "epochs": args.epochs,
                   "lr": args.lr,
                   "dropout": args.dropout,
                   "layers": args.layers,
                   "units": args.units,
                   "ref_dataset": datasets["study_id"][args.ref_dataset]})

wandb.log({"accuracy": 0.0, "loss": 0.0})

# Find all csv.gz files in the expr_matrices directory
filenames_expr = [f for f in os.listdir(
    "expr_matrices") if f.endswith("M.csv.gz")]
filenames_clusters = [f for f in os.listdir(
    "expr_matrices") if f.endswith("M_clusters.csv")]

print("Reading expression matrices...")
expr = [pd.read_csv(f"expr_matrices/{f}") for f in tqdm(filenames_expr)]

# Now intersect the gene ids
common_genes = []
for item in expr:
    item.rename(columns={item.columns[0]: "gene_id"}, inplace=True)
    item.set_index("gene_id", inplace=True)
    if len(common_genes) == 0:
        common_genes = item.index
    else:
        common_genes = common_genes.intersection(item.index)

print(f"Training on {len(common_genes)} genes, common to all datasets.")

for i in range(len(expr)):
    expr[i] = expr[i].loc[common_genes]
    # Convert to Numpy array and transpose
    expr[i] = expr[i].values.T

print("Reading clusters...")
clusters = [pd.read_csv(f"expr_matrices/{f}")
            for f in tqdm(filenames_clusters)]


def build_model(n_clusters: int, n_layers: int = 3, 
                n_nodes: List[int] = [512, 256, 16], 
                learning_rate = 0.001) -> keras.Model:
    """
    Build a 3-layer multi-class MLP classifier.

    param: n_clusters (int) - the number of clusters (possible classes)
    param: n_layers (int) - the number of hidden layers
    param: n_nodes (List[int]) - the number of nodes in each hidden layer

    return: model (keras.Model) - the model    
    """

    # Now build our MLP
    model = keras.Sequential()
    model.add(keras.layers.InputLayer(input_shape=(len(common_genes))))
    model.add(keras.layers.Dense(512, activation="relu",
              kernel_regularizer=keras.regularizers.l1(0.001), name="dense_1"))
    model.add(keras.layers.Dense(256, activation="relu",
              kernel_regularizer=keras.regularizers.l1(0.001), name="hidden_2"))
    model.add(keras.layers.Dense(32, activation="relu",
              kernel_regularizer=keras.regularizers.l1(0.001), name="hidden_3"))
    # Output layer
    # Note the +1 to take into account the "other" class
    model.add(keras.layers.Dense(n_clusters + 1,
              activation="softmax", name="output"))

    model.compile(optimizer=keras.optimizers.Adam(lr=learning_rate),
                  loss="categorical_crossentropy", metrics=["accuracy"])
    # model.summary()

    return model

def prepare_training_data(expr: List, clusters: List, dataset_id: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Prepare the training data for the MLP.

    param: expr (list) - the expression matrices
    param: clusters (list) - the list of cluster identities
    param: dataset_id (int) - the id of the dataset we are using as reference
    return: x_train (np.array) - the training data
    return: y_train (np.array) - the training labels
    """

    # We will use the expression matrices as the training data
    # The labels will be the cluster labels

    # We get the expression matrix and labels for the dataset we are using as reference...
    expr = expr[dataset_id]
    clusters = clusters[dataset_id]["Cluster"].values
    # ... and split into training and test data
    x_train, x_test, y_train, y_test = train_test_split(expr, clusters, test_size=0.1, random_state=42)

    # Scale the data with MinMaxScaler
    scaler = MinMaxScaler()
    x_train = scaler.fit_transform(x_train)
    x_test = scaler.transform(x_test)

    # One-hot encode labels
    y_train = keras.utils.to_categorical(y_train, num_classes=len(np.unique(clusters)) + 1)
    y_test = keras.utils.to_categorical(y_test, num_classes=len(np.unique(clusters)) + 1)

    return (x_train, y_train, x_test, y_test)

def get_augmented_samples(x_train: np.array, y_train: np.array, n: int = 1, perc_shuffle_genes: float = 0.01) -> Tuple[np.array, np.array]:
    """
    Creates augmented samples from the training data.

    param: x_train (np.array) - the training data
    param: y_train (np.array) - the training labels
    param: n (int) - the number of augmented samples to create. Default is 1.
    param: perc_shuffle_genes (float) - the percentage of genes to shuffle (per cell). Default is 0.01.
    return: x_train_aug, y_train_aug (np.array, np.array) - the augmented training data
    """

    # Get a list of cluster labels
    y_train_aug = np.random.choice(np.unique(np.argmax(y_train, axis=1)), size = n)

    # Shuffle the data on a per-cluster/per-gene basis

    # Split x_train by cluster
    x_by_cluster = [x_train[np.where(np.argmax(y_train, axis=1)==i)] for i in range(y_train.shape[1] - 1)]

    for x in x_by_cluster:
        genes = np.random.choice(x.shape[1], int(
            perc_shuffle_genes * x.shape[1]), replace=False)

        # Randomly swaps the values of the cells gene-wise
        for g in genes:
            x[:, g] = np.random.permutation(x[:, g])

    # Pick shuffled data corresponding to the cluster labels    
    x_train_aug = np.vstack([x_by_cluster[i][np.random.randint(0, x_by_cluster[i].shape[0]), :] for i in y_train_aug])

    # One-hot encode labels
    y_train_aug = keras.utils.to_categorical(y_train_aug, num_classes=y_train.shape[1])
    
    return (x_train_aug, y_train_aug)

ref_dataset = 0

dataset_name = filenames_expr[ref_dataset].split(".")[0]
print(f"Using dataset {dataset_name} as reference")

# Update in W&B
wandb.config.dataset_name = dataset_name

model  = build_model(len(clusters[ref_dataset]["Cluster"].unique()))
x_train, y_train, x_test, y_test = prepare_training_data(expr, clusters, ref_dataset)
x_train_aug, y_train_aug = get_augmented_samples(x_train, y_train, n=500, perc_shuffle_genes=0.01)

x_train = np.vstack([x_train, x_train_aug])
y_train = np.vstack([y_train, y_train_aug])

print(model.summary())

history = model.fit(x_train, y_train, epochs=epochs, batch_size=batch_size, validation_data=(x_test, y_test), callbacks=[WandbCallback()])