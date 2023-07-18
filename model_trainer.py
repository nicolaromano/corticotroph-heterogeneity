"""
Defines the LabelTransferTrainer class, used to train the ANN models for label transfer.

The class uses the hyperparameters defined in the hyperparameters.csv file to build a MLP model for each dataset.
This class is used by the 09_train_models.ipynb script to train the models.
"""
from typing import Tuple
import time
import tensorflow as tf
from tensorflow import keras
from keras import backend as K
from keras.regularizers import L1
import tensorflow_addons as tfa
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, KFold
import random
import os


class LabelTransferTrainer:
    """
    The class used to train the ANN models
    """

    def __init__(self, hyperparams_file: str = "hyperparameters.csv",
                 data_dir: str = "exported_matrices", output_dir: str = "model_output",
                 num_features=11103,
                 verbose: bool = True, random_seed: int = 12345) -> None:
        """
        The constructor of the class

        Parameters
        ----------
        hyperparams_file : str
            The path to the hyperparameters file, defaults to "hyperparameters.csv". 
            This is a CSV file with the following columns:
            - dataset: the name of the dataset, e.g. Cheung2018M. This will also be used as the name of the model file.            
            - batch_size: the batch size to use during training
            - learning_rate: the learning rate to use during training
            - num_nodes: the number of nodes to use in the hidden layer (all models have a single hidden layer)
            - num_classes: the number of output classes
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

        self.hyperparams = pd.read_csv(hyperparams_file, index_col=False)
        # We save verbose into a class variable so that we can directly use it in other methods
        self.verbose = verbose

        self.data_dir = data_dir
        self.output_dir = output_dir

        # Ensure directories end with a slash
        self.data_dir = self.data_dir if self.data_dir.endswith(
            '/') else f'{self.data_dir}/'
        self.output_dir = self.output_dir if self.output_dir.endswith(
            '/') else f'{self.output_dir}/'

        if (self.verbose):
            print(f"Data directory: {self.data_dir}")
            print(f"Output directory: {self.output_dir}")

        self.num_features = num_features

        # Create the output directory if it doesn't exist
        if (not os.path.exists(self.output_dir)):
            os.makedirs(self.output_dir)

        if (self.verbose):
            print(
                f"Loaded hyperparameters for {len(self.hyperparams)} models.")

        self.models = {}
        self.training_histories = {}

        self.random_seed = random_seed

        # Set the random seed for numpy, random and tensorflow
        np.random.seed(self.random_seed)
        random.seed(self.random_seed)
        tf.random.set_seed(self.random_seed)

        self._setup_models()

    def _create_model(self, learning_rate: float, num_nodes: int,
                      reg_factor: float, bn_momentum: float, num_classes: int) -> keras.Sequential:
        """
        Creates a single model.
        The model consists of an input layer, batch normalization, then a single hidden layer with ReLU activation 
        and L1 regularization, followed by a softmax output layer.

        Parameters
        ----------

        learning_rate : float
            The learning rate
        num_nodes : int
            The number of nodes in the hidden layer.
        reg_factor : float
            The regularization factor. We are using L1 regularization for all models.
        bn_momentum : float
            The momentum for the batch normalization layer.
        num_classes : int
            The number of output classes.

        Returns
        -------

        model : keras.Sequential
            The Keras model.
        """

        model = keras.Sequential()
        model.add(keras.layers.InputLayer(input_shape=(
            self.num_features,), name='input_layer'))
        model.add(keras.layers.BatchNormalization(
            name='batchnorm_layer', momentum=bn_momentum))
        model.add(keras.layers.Dense(num_nodes, activation='relu',
                  kernel_regularizer=L1(reg_factor), name='hidden_layer'))
        model.add(keras.layers.Dense(
            num_classes, activation='softmax', name='output_layer'))

        optimizer = keras.optimizers.Adam(learning_rate=learning_rate)
        model.compile(optimizer=optimizer,
                      loss='categorical_crossentropy',
                      metrics=[tfa.metrics.F1Score(num_classes=num_classes, average='macro')])

        return model

    def _setup_models(self) -> None:
        """
        Sets up the models to be trained. Models are stored in the self.models list.

        Parameters
        ----------
        None

        Returns
        -------
        None        
        """

        if (self.verbose):
            print("Setting up models...")

        for _, h in self.hyperparams.iterrows():
            # Clear the session to avoid memory leaks
            # See https://www.tensorflow.org/api_docs/python/tf/keras/backend/clear_session
            tf.keras.backend.clear_session()

            model = self._create_model(learning_rate=h['learning_rate'],
                                       num_nodes=h['num_nodes'], reg_factor=h['reg_factor'],
                                       bn_momentum=h['bn_momentum'], num_classes=h['num_classes'])

            self.models[h['dataset']] = model

            if (self.verbose):
                print(f"Created model {h['dataset']}.")

    def save_models_plots(self) -> None:
        """
        Plots the models to a <dataset_name>_model.png file in the output directory.

        Parameters
        ----------

        None

        Returns
        -------

        None, but saves the model plots to the output directory.
        """

        for name, model in self.models.items():
            keras.utils.plot_model(model, to_file=os.path.join(
                self.output_dir, f"{name}_model.png"), show_shapes=True)

    def _load_data(self, dataset: str, test_size: float = 0.9) -> Tuple[pd.DataFrame, pd.DataFrame, np.ndarray, np.ndarray]:
        """
        Loads the data for a single dataset.

        Parameters
        ----------

        dataset : str
            The name of the dataset, e.g. Cheung2018M.
        test_size : float, optional, default = 0.9
            The size of the test set. Defaults to 0.9, i.e. 90% of the data is used for training and 10% for testing.

        Returns
        -------

        (expr_data_train, expr_data_test, labels_data_train, labels_data_test) : Tuple[pd.DataFrame, pd.DataFrame, np.ndarray, np.ndarray]
            The expression data and labels for the training and test sets.
        """

        expr_data = pd.read_csv(os.path.join(
            self.data_dir, f"{dataset}_expression.csv"), index_col=0)
        expr_data = expr_data.T
        labels_data = pd.read_csv(os.path.join(
            self.data_dir, f"{dataset}_clusters.csv"), index_col=0)
        labels_data = labels_data['Cluster'].values
        labels_data = keras.utils.to_categorical(labels_data)

        expr_data_train, expr_data_test, labels_data_train, labels_data_test = train_test_split(
            expr_data, labels_data, test_size=test_size, random_state=self.random_seed)

        return (expr_data_train, expr_data_test,
                labels_data_train, labels_data_test)

    def _train_model(self,
                     model: keras.Sequential,
                     dataset: str,
                     train_expr: pd.DataFrame, train_labels: pd.DataFrame,
                     val_expr: pd.DataFrame, val_labels: pd.DataFrame,
                     epochs, batch_size) -> Tuple[keras.Sequential, dict]:
        """
        Trains a single model.

        Parameters
        ----------
        model : keras.Sequential
            The model to train.
        dataset : str
            The name of the dataset. This will be used to save the model and training history.
        train_expr, val_expr, test_expr : pd.DataFrame
            The expression data for the training, validation and test sets.
        train_labels, val_labels, test_labels : pd.DataFrame
            The labels for the training, validation and test sets.
        epochs : int, optional
            The number of epochs to train the model.
        batch_size : int, optional
            The batch size to use during training.

        Returns
        -------
        (model, history) : Tuple[keras.Sequential, dict]
            The trained model at the best epoch and the training history (full)
        """
        start_time = time.time()
        print(f"Training model {dataset} for {epochs} epochs...")

        date = time.strftime("%Y%m%d-%H%M%S")
        # Save model as <dataset_name>_best_<date>-<time>_ep<num_epochs>-f1sc<f1_score>.hdf5
        checkpoint = keras.callbacks.ModelCheckpoint(f"{self.output_dir}/{dataset}_best_{date}_ep" +
                                                     "{epoch:02d}-f1sc{val_f1_score:.2f}.hdf5",
                                                     monitor='val_f1_score',
                                                     save_best_only=True,
                                                     initial_value_threshold=0.5,
                                                     mode='max')

        history = model.fit(train_expr, train_labels,
                            epochs=epochs,
                            batch_size=batch_size,
                            validation_data=(val_expr, val_labels),
                            verbose=self.verbose,
                            callbacks=[checkpoint])

        # Find all the files in the output directory with the pattern
        # f"{self.output_dir}/{dataset}_best_{date}_ep_{epoch:02d}-f1_{val_f1_score:.2f}.hdf5"
        filenames = [f for f in os.listdir(
            self.output_dir) if f.startswith(f"{dataset}_best_")]

        # Find the file with the highest F1 score
        scores = [f.split('-')[2].split('_')[1] for f in filenames]
        scores = [float(s[:-5]) for s in scores]

        for i, f in enumerate(filenames):
            if (scores[i] < max(scores)):
                os.remove(os.path.join(self.output_dir, f))

        # Load the best model
        model = keras.models.load_model(
            os.path.join(self.output_dir, filenames[np.argmax(scores)]),
            custom_objects={'F1Score': tfa.metrics.F1Score})

        # Save history to file
        pd.DataFrame(history.history).to_csv(
            os.path.join(self.output_dir, f"{dataset}_history.csv"))

        end_time = time.time()

        if (self.verbose):
            print(
                f"Training model {dataset} took {end_time - start_time} seconds.")

        return (model, history.history)

    def train_all_models(self) -> None:
        """
        Trains all the models. The models are stored in the self.models dictionary.

        Parameters
        ----------

        None

        Returns
        -------

        None, trains the models

        """

        for i, h in self.hyperparams.iterrows():
            dataset = h['dataset']

            # Note: this is deterministic because we set the random seed
            expr_train, expr_val, _, labels_train, labels_val, _ = self._load_data(dataset)

            model = self.models[dataset]

            # Train the model. This returns the trained model at the best epoch and the training history (full)
            self.models[dataset], self.training_histories[dataset] = self._train_model(model=model, dataset=dataset,
                                                                                       train_expr=expr_train, train_labels=labels_train,
                                                                                       val_expr=expr_val, val_labels=labels_val,
                                                                                       epochs=h['epochs'], batch_size=h['batch_size'])

    def train_single_model(self, dataset: str) -> dict:
        """
        Trains a single model.
        This is a convenience method that can be used to train a single model, without having to train all models or 
        having to pass the expression data and labels to the _train_model method.

        Parameters
        ----------

        dataset : str
            The name of the dataset to train.

        Returns
        -------

        history : dict
            The training history of the model.        
        """

        expr_train, expr_val, _, labels_train, labels_val, _ = self._load_data(
            dataset)

        hyperparameters = self._get_hyperparameters(dataset)

        return self._train_model(model=self.models[dataset], dataset=dataset,
                                 train_expr=expr_train, train_labels=labels_train,
                                 val_expr=expr_val, val_labels=labels_val,
                                 epochs=hyperparameters['epochs'],
                                 batch_size=hyperparameters['batch_size'])

    def evaluate_models(self) -> dict:
        """
        Evaluates all the models on the test set.

        Parameters
        ----------

        None

        Returns
        -------

        metrics : dict
            A dictionary with the metrics for each model.
        """

        metrics = {}
        for i, row in self.hyperparams.iterrows():
            _, _, expr_test, _, _, labels_test = self._load_data(
                row['dataset'])
            metrics[row['dataset']] = self.models[row['dataset']].evaluate(
                expr_test, labels_test, batch_size=row['batch_size'])

        return metrics

    def evaluate_single_model(self, dataset: str) -> dict:
        """
        Evaluates all the models on the test set.

        Parameters
        ----------

        dataset : str
            The name of the dataset to evaluate.

        Returns
        -------

        metrics : dict
            A dictionary with the metrics.
        """

        h = self._get_hyperparameters(dataset)

        _, _, expr_test, _, _, labels_test = self._load_data(dataset)

        return self.models[dataset].evaluate(
            expr_test, labels_test, batch_size=h['batch_size'])

    def _get_hyperparameters(self, dataset: str) -> pd.DataFrame:
        """
        Returns the model for a given dataset.

        Parameters
        ----------
        dataset: str
            The name of the dataset.

        Returns
        -------
        hyperparameters : pd.DataFrame
            The hyperparameters for the given dataset.        
        """

        assert (dataset in self.hyperparams['dataset'].values)

        return self.hyperparams[self.hyperparams['dataset'] == dataset].to_dict(orient='records')[0]

    def tune_hyperparam(self, dataset: str, hyperparam: str, values: list, n_folds: int = 5) -> dict:
        """
        Tests the effect of changing an hyperparameter on the performance of a single model

        Parameters
        ----------

        dataset : str
            The name of the dataset to use.
        hyperparam : str
            The name of the hyperparameter to test. This should be one of the following:
            learning_rate, num_nodes, reg_factor, bn_momentum, batch_size            
        values : list
            The values to test.
        n_folds : int, optional, default = 5
            The number of folds to use for cross-validation.

        Returns
        -------

        The F1 score for each value of the hyperparameter.
        """
        
        assert hyperparam in {
            'learning_rate',
            'num_nodes',
            'reg_factor',
            'bn_momentum',
            'batch_size',
        }

        hyper = self._get_hyperparameters(dataset)

        tune_res = {}

        expr_train, expr_test, labels_train, labels_test = self._load_data(
            dataset, test_size=0.9)

        for v in values:
            K.clear_session()

            num_nodes = hyper['num_nodes'] if hyperparam != 'num_nodes' else v
            reg_factor = hyper['reg_factor'] if hyperparam != 'reg_factor' else v
            bn_momentum = hyper['bn_momentum'] if hyperparam != 'bn_momentum' else v            
            learning_rate = hyper['learning_rate'] if hyperparam != 'learning_rate' else v
            batch_size = hyper['batch_size'] if hyperparam != 'batch_size' else v
                
            model = self._create_model(learning_rate=learning_rate,
                                       num_nodes=num_nodes,
                                       reg_factor=reg_factor,
                                       bn_momentum=bn_momentum,
                                       num_classes=hyper['num_classes'])

            print(f"Training model {dataset} for {hyper['epochs']} epochs, with {hyperparam} = {v}...")

            kfold = KFold(n_splits=n_folds, shuffle=True, random_state=self.random_seed)
            for i, (train, val) in enumerate(kfold.split(expr_train, labels_train)):
                print(f"Fold {i+1}...")
                            
                history = model.fit(expr_train.iloc[train], labels_train[train],
                        validation_data=(expr_train.iloc[val], labels_train[val]),
                        epochs=hyper['epochs'], 
                        batch_size=hyper['batch_size'], 
                        workers=4,
                        verbose=False)

                # We get the best F1 score on the validation set.
                # Note that we completely ignore the test set here, as we don't 
                # want to use it for model selection, since that would bias the results.
                # Also, when training the final model we will use the highest scoring
                # model on the validation set (saved through checkpoints), so we
                # use the max here.                
                tune_res[f"{v}_{i}"] = max(history.history['val_f1_score'])

        return tune_res

