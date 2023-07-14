"""
10_label_transfer_predictions.py

Generates predictions for all datasets using each of the trained models.
"""

import os
# Silence TensorFlow warnings
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

from glob import glob
import pandas as pd
import numpy as np
import random
from tensorflow import keras
import tensorflow as tf
from tensorflow_addons.metrics import F1Score

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
        self.datasets = self.datasets[~self.datasets['study_id'].isin(['Ho2020M', 'Ho2020F'])]
        # Load the data
        self.expression_data = {}
        self.Seurat_labels = {}

        for dataset in self.datasets['study_id'].unique():
            if (self.verbose):
                print(f"Loading expression data for dataset {dataset}")
                
            filename = f"{self.data_dir}{dataset}_expression.csv"
            expr_data = pd.read_csv(filename, index_col=0)
            expr_data = expr_data.T
            self.expression_data[dataset] = expr_data
            
            if (self.verbose):
                print(f"Loading Seurat labels for dataset {dataset}")
                
            Seurat_labels = pd.read_csv(f"{self.data_dir}{dataset}_clusters.csv", index_col=0)
            Seurat_labels = Seurat_labels['Cluster'].values
            Seurat_labels = keras.utils.to_categorical(Seurat_labels)
            self.Seurat_labels[dataset] = Seurat_labels

        # Load the models
        self.models = {}
        for dataset in self.datasets['study_id'].unique():
            if (self.verbose):
                print(f"Loading model for dataset {dataset}")
                
            # The file is saved as <dataset_name>_best_<date>-<time>_ep<num_epochs>-f1sc<f1_score>.hdf5            
            filename = glob(f"{self.model_dir}/{dataset}_best*")
            
            if not filename:
                raise FileNotFoundError(f"No model files found for dataset {dataset}")
            
            if (len(filename) > 1):
                print(f"Multiple model files found for dataset {dataset}. Using the first one ({filename[0]})")
                
            self.models[dataset] = keras.models.load_model(filename[0], custom_objects={'f1_score': F1Score})
        
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
        
        for model_dataset in self.models.keys():
            for dataset in self.expression_data.keys():
                # Skip if the model dataset is the same as the dataset
                if (model_dataset == dataset):
                    continue
                
                if (self.verbose):
                    print(f"Predicting labels for dataset {dataset} using model trained on dataset {model_dataset}")
                
                # Get the predictions
                predictions = self.models[model_dataset].predict(self.expression_data[dataset])
                
                predictions = [np.argmax(prediction) if np.max(prediction) > threshold else -1 for prediction in predictions]
                
                # save the predictions to file
                filename = f"{self.output_dir}{dataset}_predictions_model_{model_dataset}.csv"
                
                pd.DataFrame(predictions).to_csv(filename, index=False, header=False)
                
                