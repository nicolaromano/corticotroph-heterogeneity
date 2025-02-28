"""
label_transfer_predictions.py

Generates predictions for all datasets using each of the trained models.
"""

from model_trainer import LabelTransferTrainer
import pandas as pd
import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
# Silence TensorFlow warnings
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'


class LabelTransferPredictor (LabelTransferTrainer):
    """
    The class used for generating predictions using the trained models.
    """

    def __init__(self, hyperparams_file: str = "hyperparameters.csv",
                 data_dir: str = "exported_matrices", output_dir: str = "model_output",
                 num_features=11103,
                 verbose: bool = True, random_seed: int = 12345) -> None:

        super().__init__(hyperparams_file=hyperparams_file, data_dir=data_dir, output_dir=output_dir,
                         num_features=num_features, verbose=verbose, random_seed=random_seed)

        # Load the trained models
        super().load_best_models()

    def _load_all_data(self, dataset: str) -> pd.DataFrame:
        """
        Loads the expression data and labels data for a single dataset.

        Parameters
        ----------
        dataset : str
            The name of the dataset to load the data for.

        Returns
        -------
        pd.DataFrame
            A dataframe containing the expression data for the dataset and the labels data for the dataset, with cell barcodes as the index.
        """

        expr_data = pd.read_csv(
            f"{self.data_dir}{dataset}_expression.csv", index_col=0)
        expr_data = expr_data.T
        expr_data.reset_index(drop=True, inplace=True)
        labels_data = pd.read_csv(
            f"{self.data_dir}{dataset}_clusters.csv", index_col=None)
        expr_data['Cluster'] = labels_data['Cluster']
        expr_data['Barcode'] = labels_data['Barcode']

        return expr_data

    def _load_umap(self, dataset: str) -> pd.DataFrame:
        """
        Loads the UMAP coordinates for a single dataset.

        Parameters
        ----------
        dataset : str
            The name of the dataset to load the UMAP coordinates for.

        Returns
        -------
        pd.DataFrame
            A DataFrame containing the UMAP coordinates for each cell in the dataset.
        """

        filename = f"{self.data_dir}{dataset}_umap.csv"
        return pd.read_csv(filename, index_col=None)

    def get_data(self, dataset: str) -> pd.DataFrame:
        """
        Returns the expression data and labels data for a single dataset.

        Parameters
        ----------
        dataset : str
            The name of the dataset to return the data for.

        Returns
        -------
        pd.DataFrame
            A dataframe containing the expression data for the dataset and the labels data for the dataset, with cell barcodes as the index.
        """

        return self._load_all_data(dataset)

    def predict_dataset(self, dataset: str, threshold: float = 0.7) -> pd.DataFrame:
        """
        Generates predictions for a single dataset using each of the trained models.

        Parameters
        ----------
        dataset : str
            The name of the dataset to generate predictions for.
        threshold : float, optional, default = 0.7
            The threshold used for determining whether a cell belongs to a cluster or not. If the probability of a cell
            is greater than the threshold, the cell is assigned to the cluster. Otherwise, the cell is assigned to the
            "Unknown" cluster (-1)

        Returns
        -------
        pd.DataFrame
            A DataFrame containing the predictions for each cell in the dataset. The columns are the Seurat cluster
            labels, the cell barcodes and the predicted cluster labels for each model.
        """

        expr_data = self.get_data(dataset)

        all_predictions = {}
        sex = dataset[-1]

        # Only predict using models trained on the same sex
        for model_dataset in [k for k in self.models.keys() if k[-1] == sex]:
            predictions = self.models[model_dataset].predict(
                expr_data.drop(columns=['Cluster', 'Barcode']))
            predictions = [np.argmax(prediction) if np.max(
                prediction) > threshold else -1 for prediction in predictions]
            all_predictions[f'Predicted_cluster_{model_dataset}'] = predictions

        all_predictions = pd.DataFrame(all_predictions)
        all_predictions['Cell'] = expr_data.Barcode
        all_predictions['Seurat_cluster'] = expr_data.Cluster

        return all_predictions

    def plot_predictions_on_UMAP(self, dataset: str, predictions: pd.DataFrame,
                                 display: bool = True, save: bool = False) -> None:
        """
        Plots the predictions on a UMAP plot.

        Parameters
        ----------
        dataset : str
            The name of the dataset to plot.
        predictions : pd.DataFrame
            The predictions DataFrame for the dataset.
        display : bool, optional, default = True
            If True, the plot is displayed
        save : bool, optional, default = False
            If True, the plot is saved to a file named <dataset>_predictions.pdf

        Returns
        -------
        None
        """

        umap = self._load_umap(dataset)

        color_palette = ["#FF595E", "#FF924C", "#FFCA3A", "#C5CA30",
                         "#8AC926", "#52A675", "#1982C4", "#4267AC", "#6A4C93"]

        n_plots = predictions.shape[1] - 1  # -1 because of the Cell column
        rows, cols = 2, n_plots//2 if n_plots % 2 == 0 else n_plots//2 + 1

        fig, ax = plt.subplots(rows, cols, figsize=(4*cols, 4*rows))

        for i, col in enumerate(predictions.drop('Cell', axis=1).columns):
            a = ax.flat[i]

            a.scatter(umap['UMAP_1'], umap['UMAP_2'], s=10, c=predictions[col].apply(
                lambda x: color_palette[x] if x != -1 else 'lightgray'))
            dataset_from = col.split('_')[-1]
            if col == "Seurat_cluster":
                a.set_title("Seurat Clusters")
            else:
                a.set_title(
                    f"Predictions from {dataset_from}\n% unassigned: {round(100 * np.sum(predictions[col] == -1) / len(predictions), 2)}%")
            a.set_xlabel("$UMAP_1$")
            a.set_ylabel("$UMAP_2$")
            fig.suptitle(f"Predictions for {dataset}", fontsize=16)

        # Hide the remaining axes
        for i in range(n_plots, rows*cols):
            ax.flat[i].axis('off')

        if (display):
            plt.tight_layout()
            plt.show()
        else:
            plt.close()  # Close the plot so it doesn't display when saving

        if (save):
            plt.savefig(f"{self.output_dir}{dataset}_predictions.pdf",
                        bbox_inches='tight')

    def _compute_saliency_map(self, dataset: str, input_data: np.ndarray, class_index: int) -> np.ndarray:
        """
        Compute a saliency map for the specified class index for the given input data.

        Parameters
        ----------
        dataset : str
            The name of the dataset to calculate the saliency maps for.
        input_data : np.ndarray
            The input data (gene expression)
        class_index : int
            The class index (cluster label)

        Returns
        -------
        np.ndarray : The saliency map
        """

        input_tensor = tf.convert_to_tensor(input_data, dtype=tf.float32)
            
        with tf.GradientTape() as tape:
            tape.watch(input_tensor)
            # Forward pass through the model
            output = self.models[dataset](input_tensor, training=False)
            class_output = tf.reduce_sum(output[:, class_index])
        # Compute the gradient of the output with respect to the input
        gradients = tape.gradient(class_output, input_tensor)            

        return gradients.numpy()

    def get_saliency(self, dataset: str) -> (dict, dict):
        """
        Calculates the saliency maps for a single dataset. This is done by calculating the gradient of the model output with
        respect to the input. The saliency map is then calculated by taking the gradient and averaging over the samples for each cluster.
        Saliency values are normalized using z-scores across genes per cluster.

        Parameters
        ----------
        dataset : str
            The name of the dataset to calculate the saliency maps for.

        Returns
        -------
        saliency_maps_summary : dict
            A dictionary containing the normalized saliency maps for each cluster. The keys are the cluster labels, and the values are
            the normalized average and SD of saliency for each gene (pd.DataFrame).
        all_saliency_maps : dict
            A dictionary containing the normalized saliency maps for each cluster. The keys are the cluster labels, and the values are
            the normalized saliency maps for each sample and gene (pd.DataFrame).
        """

        saliency_maps_summary = {}
        all_saliency_maps = {}
        expr_data = self._load_all_data(dataset)

        for cluster in expr_data['Cluster'].unique():
            cluster_indices = expr_data['Cluster'] == cluster
            input_data = expr_data.drop(columns=['Cluster', 'Barcode'])[
                cluster_indices.values]

            # Compute saliency maps for all samples in the cluster
            saliency_map = [self._compute_saliency_map(
                dataset, input_data[i:i+1], cluster) for i in range(input_data.shape[0])]

            # Convert saliency maps to NumPy array for aggregation
            saliency_map = np.array(saliency_map)

            # Compute average and standard deviation of saliency for each gene
            saliency_map_avg = np.mean(saliency_map, axis=0)[0]
            saliency_map_sd = np.std(saliency_map, axis=0)[0]

            # Normalize using z-score across genes (per cluster)
            mean_per_gene = saliency_map_avg.mean()
            std_per_gene = saliency_map_avg.std()

            normalized_saliency_map = (saliency_map - mean_per_gene) / std_per_gene

            # Compute normalized average and SD for summary
            normalized_saliency_mean = (saliency_map_avg - mean_per_gene) / std_per_gene
            normalized_saliency_sd = saliency_map_sd / std_per_gene  # Normalize SD consistently

            # Create DataFrame for normalized saliency map summary
            saliency_map_df = pd.DataFrame(
                {
                    'Saliency_mean': normalized_saliency_mean,
                    'Saliency_sd': normalized_saliency_sd,
                    'Cluster': cluster
                },
                index=input_data.columns
            )
            saliency_map_df = saliency_map_df.sort_values(
                by='Saliency_mean', ascending=False)

            saliency_map_df.index.rename("Gene", inplace=True)

            # Save summary and per-sample normalized maps
            saliency_maps_summary[cluster] = saliency_map_df
            all_saliency_maps[cluster] = pd.DataFrame(
                normalized_saliency_map[:, 0],  # Convert normalized array to DataFrame
                columns=input_data.columns,
                index=[f"Sample_{i}" for i in range(normalized_saliency_map.shape[0])]
            )

        return saliency_maps_summary, all_saliency_maps

