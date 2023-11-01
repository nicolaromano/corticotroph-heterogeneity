"""
label_transfer_predictions.py

Generates predictions for all datasets using each of the trained models.
"""

import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
from tensorflow_addons.metrics import F1Score
import tensorflow as tf
from tensorflow import keras
import random
import numpy as np
import pandas as pd
from glob import glob
import os
# Silence TensorFlow warnings
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'


class LabelTransferPredictor:
    """
    The class used for generating predictions using the trained models.
    """

    def __init__(self, data_dir: str = "exported_matrices/", output_dir: str = "label_transfer_model_output/predictions/",
                 model_dir: str = "label_transfer_model_output/",
                 verbose: bool = True, random_seed: int = 12345) -> None:
        """
        The constructor of the class

        Parameters
        ----------
        data_dir : str
            The path to the directory where the data is stored. For each dataset, there should be a file called
            <dataset_name>_expression.csv and <dataset_name>_clusters.csv
        output_dir : str, optional, default = "model_output"
            The path to the directory where the trained models and training history will be saved.
        num_features : int, optional, default = 11103
            The number of features in the dataset. This can be safely left as the default value as all datasets have
            the same number of features.
        verbose : bool, optional, default = True
            Whether to print output messages or not.
        random_seed : int, optional, default = 12345
            The random seed used e.g. for splitting the data into training, validation and test sets and for
            initializing the weights of the models.

        Returns
        -------
        None            
        """

        # We save verbose into a class variable so that we can directly use it in other methods
        self.verbose = verbose

        self.data_dir = data_dir
        self.output_dir = output_dir

        # Ensure directories end with a slash
        self.data_dir = self.data_dir if self.data_dir.endswith(
            '/') else f'{self.data_dir}/'
        self.output_dir = self.output_dir if self.output_dir.endswith(
            '/') else f'{self.output_dir}/'
        self.model_dir = model_dir if model_dir.endswith(
            '/') else f'{model_dir}/'

        if (self.verbose):
            print(f"Data directory: {self.data_dir}")
            print(f"Output directory: {self.output_dir}")

        # Create the output directory if it doesn't exist
        if (not os.path.exists(self.output_dir)):
            os.makedirs(self.output_dir)

        self.random_seed = random_seed
        # Set the random seed for numpy, random and tensorflow
        np.random.seed(self.random_seed)
        random.seed(self.random_seed)
        tf.random.set_seed(self.random_seed)

        self.datasets = pd.read_csv('datasets.csv', index_col=None)
        # Remove Ho
        self.datasets = self.datasets[~self.datasets['study_id'].isin(
            ['Ho2020M', 'Ho2020F'])]
        # Load the data
        self.expression_data = {}
        self.Seurat_labels = {}
        self.umaps = {}
        self.models = {}
        self._load_expression_data()
        self._load_models()

        self.predictions = {}

    def _load_expression_data(self) -> None:
        """
        Loads the expression data, Seurat labels and UMAP coordinates for each dataset.

        Parameters
        ----------

        None

        Returns
        -------

        None, the expression data, Seurat labels are saved in the self.expression_data and self.Seurat_labels
        """

        for dataset in self.datasets['study_id'].unique():
            if (self.verbose):
                print(f"Loading expression data for dataset {dataset}")

            filename = f"{self.data_dir}{dataset}_expression.csv"
            expr_data = pd.read_csv(filename, index_col=0)
            expr_data = expr_data.T
            self.expression_data[dataset] = expr_data

            if (self.verbose):
                print(f"Loading Seurat labels for dataset {dataset}")

            Seurat_labels = pd.read_csv(
                f"{self.data_dir}{dataset}_clusters.csv", index_col=None)
            self.Seurat_labels[dataset] = Seurat_labels

            if (self.verbose):
                print(f"Loading UMAP coordinates for dataset {dataset}")

            umap = pd.read_csv(
                f"{self.data_dir}{dataset}_umap.csv", index_col=None)
            self.umaps[dataset] = umap

    def _load_models(self) -> None:
        """
        Loads the trained models from the model directory.

        Parameters
        ----------

        None

        Returns
        -------

        None, the models are saved in the self.models dictionary.                
        """

        for dataset in self.datasets['study_id'].unique():
            if (self.verbose):
                print(f"Loading model for dataset {dataset}")

            # The file is saved as <dataset_name>_best_<date>-<time>_ep<num_epochs>-f1sc<f1_score>.hdf5
            filename = glob(f"{self.model_dir}/{dataset}_best*")

            if not filename:
                raise FileNotFoundError(
                    f"No model files found for dataset {dataset}")

            if (len(filename) > 1):
                print(
                    f"Multiple model files found for dataset {dataset}. Using the first one ({filename[0]})")

            self.models[dataset] = keras.models.load_model(
                filename[0], custom_objects={'f1_score': F1Score})

        if (self.verbose):
            print(f"Loaded {len(self.models)} models.")

    def predict_all(self, threshold: float = 0.7):
        """
        Generates predictions for all datasets using each of the trained models.

        Parameters
        ----------
        threshold : float, optional, default = 0.7
            The threshold used for determining whether a cell belongs to a cluster or not. If the probability of a cell
            is greater than the threshold, the cell is assigned to the cluster. Otherwise, the cell is assigned to the
            "Unknown" cluster (-1)

        Returns
        -------
        None
        """

        for dataset in self.expression_data.keys():
            sex = dataset[-1]
            all_predictions = {}
            filename = f"{self.output_dir}{dataset}_predictions.csv"

            for model_dataset in self.models.keys():
                # Only predict using models trained on the same sex
                if (model_dataset[-1] != sex):
                    continue

                if (self.verbose):
                    print(
                        f"Predicting labels for dataset {dataset} using model trained on dataset {model_dataset}")

                # Get the predictions
                predictions = self.models[model_dataset].predict(
                    self.expression_data[dataset])

                predictions = [np.argmax(prediction) if np.max(
                    prediction) > threshold else -1 for prediction in predictions]

                if (all_predictions == {}):
                    all_predictions = {'Cell': self.Seurat_labels[dataset]['Barcode'],
                                       'Seurat_cluster': self.Seurat_labels[dataset]['Cluster'],
                                       f'Predicted_cluster_{model_dataset}': predictions}
                else:
                    all_predictions[f'Predicted_cluster_{model_dataset}'] = predictions

            # Save the predictions to file
            all_predictions = pd.DataFrame(all_predictions)
            self.predictions[dataset] = all_predictions
            all_predictions.to_csv(filename, index=False, header=True)

    def plot_predictions_umap(self, dataset: str, save: bool = False) -> None:
        """
        Plots the predictions on a UMAP plot.

        Parameters
        ----------

        dataset : str
            The name of the dataset to plot.
            
        save : bool, optional, default = False. If True, the plot is saved to a file named <dataset>_predictions.pdf

        Returns
        -------

        None                
        """

        if (dataset not in self.predictions.keys()):
            print(
                f"Dataset {dataset} not found in predictions. Have you run the predict_all() method yet?")
            print("Available datasets:")
            for key in self.predictions.keys():
                print(key)
            return

        color_palette = ["#FF595E", "#FF924C", "#FFCA3A", "#C5CA30",
                         "#8AC926", "#52A675", "#1982C4", "#4267AC", "#6A4C93"]

        nplots = self.predictions[dataset].shape[1] - 1 # -1 because of the Cell column
        nr, nc = 2, nplots//2 if nplots % 2 == 0 else nplots//2 + 1

        _, ax = plt.subplots(nr, nc, figsize=(5*nplots//2, 10))

        # First plot is the Seurat clusters
        ax[0, 0].scatter(self.umaps[dataset]['UMAP_1'], self.umaps[dataset]['UMAP_2'],
                         s=10, c=self.Seurat_labels[dataset]['Cluster'].apply(lambda x: color_palette[x]))
        ax[0, 0].set_title("Seurat clusters")

        for i, a in enumerate(ax.flatten()[1:]):
            if i < (nplots - 1):
                model_dataset = list(self.predictions[dataset].columns)[
                    i+2].split('_')[2:]
                # To avoid Ruf_Zamojski to be split into two words
                model_dataset = '_'.join(model_dataset)

                predictions = self.predictions[dataset][f'Predicted_cluster_{model_dataset}']
                # Color the unassigned cells lightgray
                colors = [color_palette[prediction] if prediction != -1
                        else 'lightgray' for prediction in predictions]
                a.scatter(self.umaps[dataset]['UMAP_1'], self.umaps[dataset]['UMAP_2'],
                        s=10, c=colors)
                a.set_title(
                    f"Predictions from {model_dataset}\n% unassigned: {round(100 * np.sum(predictions == -1) / len(predictions), 2)}%")
                a.set_xlabel("$UMAP_1$")
                a.set_ylabel("$UMAP_2$")
            else:
                a.axis('off')

        if (save):
            plt.savefig(f"{self.output_dir}{dataset}_predictions.pdf",
                        bbox_inches='tight')
        else:
            plt.tight_layout()

    def compute_saliency_map(self, dataset: str, input_data: np.ndarray, class_index: int) -> np.ndarray:
        """
        Compute a saliency map for the specified class index for the given input data.

        Parameters
        ----------
        model : keras.Model
            The trained Keras model
        input_data : np.ndarray
            The input data (gene expression)
        class_index : int
            The class index (cluster label)

        Returns
        -------
        np.ndarray : The saliency map
        """

        input_tensor = tf.convert_to_tensor(input_data)

        with tf.GradientTape() as tape:
            tape.watch(input_tensor)
            # Forward pass through the model
            output = self.models[dataset](input_tensor, training=False)
            class_output = output[:, class_index]

        # Compute the gradient of the output with respect to the input
        gradients = tape.gradient(class_output, input_tensor)

        return gradients.numpy()


    def get_saliency(self, dataset:str) -> dict:
        """
        Calculates the saliency maps for the dataset. This is done by calculating the gradient of the model output with
        respect to the input. The saliency map is then calculated by taking the gradient and averaging over the samples for each cluster.
        We get a value per input (gene)
        
        Parameters
        ----------
        
        dataset : str
            The name of the dataset to calculate the saliency maps for.
            
        Returns
        -------
        
        saliency_maps : dict
            A dictionary containing the saliency maps for each cluster. The keys are the cluster labels and the values are
            the average and sd of saliency for each gene (pd.DataFrame)
        """
        
        saliency_maps = {}
        
        for cluster in self.Seurat_labels[dataset]['Cluster'].unique():
            cluster_indices = self.Seurat_labels[dataset]['Cluster'] == cluster
            input_data = self.expression_data[dataset][cluster_indices.values]
            
            saliency_map = [self.compute_saliency_map(dataset, input_data[i:i+1], cluster) for i in range(input_data.shape[0])]
            
            saliency_map_avg = np.mean(saliency_map, axis=0)[0]
            saliency_map_sd = np.std(saliency_map, axis=0)[0]

            saliency_map_df = pd.DataFrame(
                {'Saliency_mean': saliency_map_avg, 'Saliency_sd': saliency_map_sd, 'Cluster': cluster}, index=input_data.columns)
            saliency_map_df = saliency_map_df.sort_values(
                by='Saliency_mean', ascending=False)
            
            saliency_map_df.index.rename("Gene", inplace=True)
            saliency_maps[cluster] = saliency_map_df
            
        return saliency_maps