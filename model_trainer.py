"""
Defines the LabelTransferTrainer class, used to train the ANN models for label transfer.

The class uses the hyperparameters defined in the hyperparameters.csv file to build a MLP model for each dataset.
This class is used by the 09_train_models.ipynb script to train the models.
"""
import time
import tensorflow as tf
import keras
from keras import backend as K
from keras.regularizers import L1, L2, L1L2
from keras.callbacks import ModelCheckpoint, LearningRateScheduler
from keras.metrics import F1Score, AUC
# import tensorflow_addons as tfa
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, KFold
from sklearn.preprocessing import MinMaxScaler
import random
import os
from math import exp

os.environ["KERAS_BACKEND"] = "tensorflow"


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

    def _create_model(self, learning_rate: float, num_layers: int, num_nodes: int,
                      reg_factor: float, reg_type: str, dropout_rate: float,
                      bn_momentum: float, num_classes: int) -> keras.Sequential:
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
        num_layers : int
            The number of BN + hidden layers to use.
        reg_factor : float
            The regularization factor.
        reg_type : str
            The type of regularization to use. Can be 'l1', 'l2' or 'l1l2'.
        dropout_rate : float
            The dropout rate to use. If set to 0, no dropout layer will be added.
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
        model.add(keras.layers.InputLayer(shape=(
            self.num_features,), name='input_layer'))

        for i in range(num_layers):
            model.add(keras.layers.BatchNormalization(
                name=f'batchnorm_layer_{i}', momentum=bn_momentum))
            kernel_regularizer = L1 if reg_type == 'l1' else L2 if reg_type == 'l2' else L1L2
            model.add(keras.layers.Dense(num_nodes, activation='relu',
                      kernel_regularizer=kernel_regularizer(reg_factor), name=f'hidden_layer_{i}'))
            if (dropout_rate > 0):
                model.add(keras.layers.Dropout(
                    dropout_rate, name=f'dropout{i}'))

        model.add(keras.layers.Dense(
            num_classes, activation='softmax', name='output_layer'))

        optimizer = keras.optimizers.Adam(learning_rate=learning_rate)
        model.compile(optimizer=optimizer,
                      loss='categorical_crossentropy',
                      metrics=['accuracy', F1Score(average='micro', name='micro_f1_score'),
                               F1Score(average='macro',
                                       name='macro_f1_score'),
                               AUC(num_thresholds=100, name='ROCAUC')])
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

            model = self._create_model(learning_rate=h['learning_rate'], num_layers=h['num_layers'],
                                       num_nodes=h['num_nodes'], reg_factor=h['reg_factor'],
                                       reg_type=h['reg_type'], dropout_rate=h['dropout_rate'],
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

    def _load_data(self, dataset: str, test_size: float = 0.1) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Loads the data for a single dataset.

        Parameters
        ----------

        dataset : str
            The name of the dataset, e.g. Cheung2018M.
        test_size : float, optional, default = 0.1
            The size of the test set. Defaults to 0.1, i.e. 10% of the data is used for training and 90% for testing.

        Returns
        -------

        (expr_data_train, expr_data_test, labels_data_train, labels_data_test) : tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
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
            expr_data, labels_data, train_size=1.0-test_size, random_state=self.random_seed)

        # Min-max scaling
        scaler = MinMaxScaler()
        expr_data_train = scaler.fit_transform(expr_data_train)
        expr_data_test = scaler.transform(expr_data_test)

        return (expr_data_train, expr_data_test,
                labels_data_train, labels_data_test)

    def clear_saved_models(self) -> None:
        """
        Clears the saved models and training history from the output directory.

        Parameters
        ----------

        None

        Returns
        -------

        None
        """

        for f in os.listdir(self.output_dir):
            if (f.endswith('.keras') or f.endswith('.csv')) and "_best_" in f:
                os.remove(os.path.join(self.output_dir, f))
            
            if (f.endswith('history.csv')):
                os.remove(os.path.join(self.output_dir, f))

    def load_best_models(self, clear_others: bool = False) -> pd.DataFrame:
        """
        Loads the best models from the output directory.

        Parameters
        ----------
        clear_others : bool, optional, default = False
            Whether to clear the other (non-best) models from the output directory or not.

        Returns
        -------
        best_models_dict : pd.DataFrame
            A DataFrame with the best models for each dataset. The columns are:
            - filename: the filename of the model
            - f1: the F1 score of the model
            - date: the date the model was trained
        """
        best_models_dict = {}
        all_models = [f for f in os.listdir(self.output_dir) if (
            f.endswith('.keras') and "_best_" in f)]

        studies = set([m.split("_best_")[0]
                      for m in all_models])  # Get the name of the studies
        for s in studies:
            models = [m for m in all_models if m.startswith(s)]
            # Get F1 and date from the model name
            f1 = [m.split("f1_")[1].replace(".keras", "") for m in models]
            date = [m.split("_best_")[1].split("_ep")[0] for m in models]
            # Get the model with the highest F1 score. If there are multiple, get the latest one
            # Sort the models by date, with latest first
            models = [m for _, m in sorted(zip(date, models), reverse=True)]
            best_model = models[np.argmax([float(f) for f in f1])]
            self.models[s] = keras.models.load_model(
                os.path.join(self.output_dir, best_model))
            self.training_histories[s] = pd.read_csv(
                os.path.join(self.output_dir, f"{s}_{date[np.argmax([float(f) for f in f1])]}_history.csv"))
            if clear_others:
                for m in models:
                    if m != best_model:
                        os.remove(os.path.join(self.output_dir, m))
                        os.remove(os.path.join(
                            self.output_dir, f"{s}_{date[np.argmax([float(f) for f in f1])]}_history.csv"))

            if self.verbose:
                print(f"Loaded best model for {s} from {best_model}")
            best_models_dict[s] = {
                "filename": best_model, "f1": f1[np.argmax([float(f) for f in f1])], "date": date[np.argmax([float(f) for f in f1])]}

        return pd.DataFrame(best_models_dict).T

    def _train_model(self,
                     model: keras.Sequential,
                     dataset: str,
                     train_expr: pd.DataFrame, train_labels: pd.DataFrame,
                     epochs, batch_size) -> tuple[keras.Sequential, dict]:
        """
        Trains a single model.

        Parameters
        ----------
        model : keras.Sequential
            The model to train.
        dataset : str
            The name of the dataset. This will be used to save the model and training history.
        train_expr : pd.DataFrame
            The expression data for the training set. A 20% validation split will be used.
        train_labels : pd.DataFrame
            The labels for the training set.
        epochs : int, optional
            The number of epochs to train the model.
        batch_size : int, optional
            The batch size to use during training.

        Returns
        -------
        (model, history) : tuple[keras.Sequential, dict]
            The trained model at the best epoch and the training history (full)
        """
        start_time = time.time()
        print(f"Training model {dataset} for {epochs} epochs...")

        date = time.strftime("%Y%m%d-%H%M%S")
        # Save model as <dataset_name>_best_<date>-<time>_ep<num_epochs>-f1sc<macro_f1_score>.keras
        checkpoint = ModelCheckpoint(f"{self.output_dir}{dataset}_best_{date}_ep" +
                                     "{epoch:02d}-f1_{val_macro_f1_score:.2f}.keras",
                                     monitor='val_macro_f1_score',
                                     save_best_only=True,
                                     initial_value_threshold=0.3,
                                     mode='max')

        decay_rate = self._get_hyperparameters(dataset)['lr_exp_decay']
        initial_lr = self._get_hyperparameters(dataset)['learning_rate']

        def lr_exp_decay(epoch):
            return initial_lr * np.exp(-decay_rate * epoch)

        lr_scheduler_callback = LearningRateScheduler(
            lambda epoch: lr_exp_decay(epoch))

        history = model.fit(train_expr, train_labels,
                            validation_split=0.2,
                            epochs=epochs,
                            batch_size=batch_size,
                            verbose="auto",
                            callbacks=[checkpoint, lr_scheduler_callback])

        # Find all the files in the output directory with the pattern
        # f"{self.output_dir}/{dataset}_best_{date}_ep_{epoch:02d}-f1_{macro_f1_score:.2f}.keras"
        filenames = [f for f in os.listdir(
            self.output_dir) if f.startswith(f"{dataset}_best_")]

        if (len(filenames) == 0):
            raise Exception(
                f"No files found in {self.output_dir} for dataset {dataset}.")

        # Find the file with the highest F1 score
        scores = [f.split('f1_')[1] for f in filenames]
        scores = [float(f.split('.keras')[0]) for f in scores]

        for i, f in enumerate(filenames):
            if (scores[i] < max(scores)):
                os.remove(os.path.join(self.output_dir, f))

        if (self.verbose):
            print(
                f"Loading best model for {dataset} from {filenames[np.argmax(scores)]}...")

        # Load the best model
        model = keras.models.load_model(filepath=os.path.join(
            self.output_dir, filenames[np.argmax(scores)]))

        # Save history to file
        pd.DataFrame(history.history).to_csv(
            os.path.join(self.output_dir, f"{dataset}_{date}_history.csv"), index=False)

        end_time = time.time()

        if (self.verbose):
            print(
                f"Training model {dataset} took {end_time - start_time} seconds.")

        return (model, history.history)

    def train_all_models(self, reset_models: bool = False) -> None:
        """
        Trains all the models. The models are stored in the self.models dictionary.

        Parameters
        ----------

        reset_models : bool, optional, default = False
            Whether to reset the models before training them. If set to True, the models will be re-created from scratch.

        Returns
        -------

        None, trains the models

        """
        if (reset_models):
            self._setup_models()

        for i, h in self.hyperparams.iterrows():
            dataset = h['dataset']

            # Note: this is deterministic because we set the random seed
            expr_train, _, labels_train, _ = self._load_data(dataset)

            model = self.models[dataset]

            # Train the model. This returns the trained model at the best epoch and the training history (full)
            self.models[dataset], self.training_histories[dataset] = self._train_model(model=model, dataset=dataset,
                                                                                       train_expr=expr_train, train_labels=labels_train,
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

        expr_train, _, labels_train, _ = self._load_data(dataset)

        hyperparameters = self._get_hyperparameters(dataset)

        return self._train_model(model=self.models[dataset], dataset=dataset,
                                 train_expr=expr_train, train_labels=labels_train,
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
            _, expr_test, _, labels_test = self._load_data(
                row['dataset'])
            metrics[row['dataset']] = self.models[row['dataset']].evaluate(
                expr_test, labels_test, batch_size=row['batch_size'])
            metrics_names = ['accuracy', 'micro_f1_score',
                             'macro_f1_score', 'ROCAUC']
            metrics[row['dataset']] = dict(
                zip(metrics_names, metrics[row['dataset']]))

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

        _, expr_test, _, labels_test = self._load_data(dataset)

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

        assert (
            dataset in self.hyperparams['dataset'].values), f"Dataset {dataset} not found. Possible values are {self.hyperparams['dataset'].values}"

        return self.hyperparams[self.hyperparams['dataset'] == dataset].to_dict(orient='records')[0]

    def tune_hyperparam(self, dataset: str, hyperparam: str, values: list, n_folds: int = 5) -> pd.DataFrame:
        """
        Tests the effect of changing an hyperparameter on the performance of a single model

        Parameters
        ----------

        dataset : str
            The name of the dataset to use.
        hyperparam: str
            The name of the hyperparameter to test. This should be one of the following:
            learning_rate, num_layers, num_nodes, reg_factor, reg_type, dropout_rate, bn_momentum, batch_size
        values : list
            The values to test for the hyperparameter
        n_folds : int, optional, default = 5
            The number of folds to use for cross-validation.

        Returns
        -------

        tune_res : pd.DataFrame containing the results of the hyperparameter tuning.
        """

        valid_hyperparameters = ['learning_rate', 'num_layers', 'num_nodes', 'reg_factor',
                                 'reg_type', 'dropout_rate', 'bn_momentum', 'batch_size']

        assert hyperparam in valid_hyperparameters

        hyper = self._get_hyperparameters(dataset)

        expr_train, _, labels_train, _ = self._load_data(
            dataset, test_size=0.9)

        tune_res = pd.DataFrame()

        for v in values:
            num_layers = hyper['num_layers'] if hyperparam != 'num_layers' else v
            num_nodes = hyper['num_nodes'] if hyperparam != 'num_nodes' else v
            reg_factor = hyper['reg_factor'] if hyperparam != 'reg_factor' else v
            reg_type = hyper['reg_type'] if hyperparam != 'reg_type' else v
            bn_momentum = hyper['bn_momentum'] if hyperparam != 'bn_momentum' else v
            dropout_rate = hyper['dropout_rate'] if hyperparam != 'dropout_rate' else v
            learning_rate = hyper['learning_rate'] if hyperparam != 'learning_rate' else v
            batch_size = hyper['batch_size'] if hyperparam != 'batch_size' else v

            print(
                f"Training model {dataset} for {hyper['epochs']} epochs, with {hyperparam} = {v}...")

            kfold = KFold(n_splits=n_folds, shuffle=True,
                          random_state=self.random_seed)

            for i, (train, val) in enumerate(kfold.split(expr_train, labels_train)):
                print(f"Fold {i+1}", end="...")

                # Reset model
                K.clear_session()
                model = self._create_model(learning_rate=learning_rate,
                                           num_layers=num_layers,
                                           num_nodes=num_nodes,
                                           reg_factor=reg_factor,
                                           reg_type=reg_type,
                                           dropout_rate=dropout_rate,
                                           bn_momentum=bn_momentum,
                                           num_classes=hyper['num_classes'])

                history = model.fit(expr_train[train], labels_train[train],
                                    validation_data=(
                                        expr_train[val], labels_train[val]),
                                    epochs=hyper['epochs'],
                                    batch_size=batch_size,
                                    verbose=self.verbose)

                print("max F1 score:", max(
                    history.history['val_macro_f1_score']))

                tune_res = pd.concat([tune_res, pd.DataFrame({
                    hyperparam: [v] * hyper['epochs'],
                    "Fold": [i+1] * hyper['epochs'],
                    "Epoch": range(1, hyper['epochs']+1),
                    "Loss": history.history['loss'],
                    "Macro F1 score": history.history['macro_f1_score'],
                    "Micro F1 score": history.history['micro_f1_score'],
                    "ROCAUC": history.history['ROCAUC'],
                    "Val loss": history.history['val_loss'],
                    "Val Macro F1 score": history.history['val_macro_f1_score'],
                    "Val Micro F1 score": history.history['val_micro_f1_score'],
                    "Val ROCAUC": history.history['val_ROCAUC']
                })])

                # Enforce types
                tune_res = tune_res.astype({
                    hyperparam: float if hyperparam not in {'reg_type'} else str,
                    "Fold": int,
                    "Epoch": int,
                    "Loss": float,
                    "Macro F1 score": float,
                    "Micro F1 score": float,
                    "ROCAUC": float,
                    "Val loss": float,
                    "Val Macro F1 score": float,
                    "Val Micro F1 score": float,
                    "Val ROCAUC": float
                })

        return tune_res

    def tune_two_hyperparameters(self, dataset: str,
                                 hyperparam1: str, hyperparam2: str, values1: list, values2: list, n_folds: int = 5) -> pd.DataFrame:
        """
        Tests the effect of changing two hyperparameters on the performance of a single model

        Parameters
        ----------

        dataset: str
            The name of the dataset to use.
        hyperparam1, 2: str
            The name of the hyperparameters to test. This should be one of the following:
            learning_rate, num_layers, num_nodes, reg_factor, reg_type, dropout_rate, bn_momentum, batch_size
        values1, 2 : list
            The values to test for the hyperparameters
        n_folds : int, optional, default = 5
            The number of folds to use for cross-validation.

        Returns
        -------

        tune_res : pd.DataFrame containing the results of the hyperparameter tuning.
        """

        valid_hyperparameters = ['learning_rate', 'num_layers', 'num_nodes', 'reg_factor',
                                 'reg_type', 'dropout_rate', 'bn_momentum', 'batch_size']

        if not hyperparam1 in valid_hyperparameters:
            raise ValueError(f"Invalid hyperparameter {hyperparam1}")
        if not hyperparam2 in valid_hyperparameters:
            raise ValueError(f"Invalid hyperparameter {hyperparam2}")

        if hyperparam1 == hyperparam2:
            raise ValueError(
                "The two hyperparameters to tune should be different")

        hyper = self._get_hyperparameters(dataset)
        expr_train, _, labels_train, _ = self._load_data(
            dataset, test_size=0.9)

        tune_res = pd.DataFrame()

        for v1 in values1:
            for v2 in values2:
                num_layers = v1 if hyperparam1 == 'num_layers' else v2 if hyperparam2 == 'num_layers' else hyper[
                    'num_layers']
                num_nodes = v1 if hyperparam1 == 'num_nodes' else v2 if hyperparam2 == 'num_nodes' else hyper[
                    'num_nodes']
                reg_factor = v1 if hyperparam1 == 'reg_factor' else v2 if hyperparam2 == 'reg_factor' else hyper[
                    'reg_factor']
                reg_type = v1 if hyperparam1 == 'reg_type' else v2 if hyperparam2 == 'reg_type' else hyper[
                    'reg_type']
                bn_momentum = v1 if hyperparam1 == 'bn_momentum' else v2 if hyperparam2 == 'bn_momentum' else hyper[
                    'bn_momentum']
                dropout_rate = v1 if hyperparam1 == 'dropout_rate' else v2 if hyperparam2 == 'dropout_rate' else hyper[
                    'dropout_rate']
                learning_rate = v1 if hyperparam1 == 'learning_rate' else v2 if hyperparam2 == 'learning_rate' else hyper[
                    'learning_rate']
                batch_size = v1 if hyperparam1 == 'batch_size' else v2 if hyperparam2 == 'batch_size' else hyper[
                    'batch_size']

                print(
                    f"Training model {dataset} for {hyper['epochs']} epochs, with {hyperparam1} = {v1} and {hyperparam2} = {v2}")

                kfold = KFold(n_splits=n_folds, shuffle=True,
                              random_state=self.random_seed)

                for i, (train, val) in enumerate(kfold.split(expr_train, labels_train)):
                    print(f"Fold {i+1}", end="...")

                    # Reset model
                    K.clear_session()
                    model = self._create_model(learning_rate=learning_rate,
                                               num_layers=num_layers,
                                               num_nodes=num_nodes,
                                               reg_factor=reg_factor,
                                               reg_type=reg_type,
                                               dropout_rate=dropout_rate,
                                               bn_momentum=bn_momentum,
                                               num_classes=hyper['num_classes'])

                    history = model.fit(expr_train[train], labels_train[train],
                                        validation_data=(
                                            expr_train[val], labels_train[val]),
                                        epochs=hyper['epochs'],
                                        batch_size=batch_size,
                                        verbose=self.verbose)

                    print("max macro F1 score:", max(
                        history.history['val_macro_f1_score']))

                    tune_res = pd.concat([tune_res, pd.DataFrame({
                        hyperparam1: [v1] * hyper['epochs'],
                        hyperparam2: [v2] * hyper['epochs'],
                        "Fold": [i+1] * hyper['epochs'],
                        "Epoch": range(1, hyper['epochs']+1),
                        "Loss": history.history['loss'],
                        "Macro F1 score": history.history['macro_f1_score'],
                        "Micro F1 score": history.history['micro_f1_score'],
                        "ROCAUC": history.history['ROCAUC'],
                        "Val loss": history.history['val_loss'],
                        "Val Macro F1 score": history.history['val_macro_f1_score'],
                        "Val Micro F1 score": history.history['val_micro_f1_score'],
                        "Val ROCAUC": history.history['val_ROCAUC']
                    })])

                    # Enforce types
                    tune_res = tune_res.astype({
                        hyperparam1: float if hyperparam1 not in {'reg_type'} else str,
                        hyperparam2: float if hyperparam2 not in {'reg_type'} else str,
                        "Fold": int,
                        "Epoch": int,
                        "Loss": float,
                        "Macro F1 score": float,
                        "Micro F1 score": float,
                        "ROCAUC": float,
                        "Val loss": float,
                        "Val Macro F1 score": float,
                        "Val Micro F1 score": float,
                        "Val ROCAUC": float
                    })

        return tune_res
