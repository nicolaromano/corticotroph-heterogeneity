"""
Defines the LabelTransferTrainer class, used to train the XGBoost models for label transfer.

The class uses the hyperparameters defined in the hyperparameters.csv file to build a XGBoost model for each dataset.
This class is used by the 09_train_models.ipynb script to train the models.
"""

import time
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score
import xgboost as xgb
import random
import os


class LabelTransferTrainer:
    """
    The class used to train the XGBoost models
    """

    def __init__(self, hyperparams_file: str = "hyperparameters.csv",
                 data_dir: str = "exported_matrices", output_dir: str = "label_transfer_model_output",
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
            - learning_rate: the learning rate eta to use during training
            - max_depth: the maximum depth of the trees
            - n_boost_rounds: the number of boosting rounds            
            - gamma: the gamma parameter, used for regularization
            - num_classes: the number of output classes
        data_dir : str
            The path to the directory where the data is stored. For each dataset, there should be a file called
            <dataset_name>_expression.csv and <dataset_name>_clusters.csv
        output_dir : str, optional, default = "xgboost_model_output"
            The path to the directory where the trained models and training history will be saved.
        num_features : int, optional, default = 11103
            The number of features in the dataset. This can be safely left as the default value as all datasets have
            been saved with the same number of features.
        verbose : bool, optional, default = True
            Whether to print output messages or not.
        random_seed : int, optional, default = 12345
            The random seed used e.g. for splitting the data into training and test sets.

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

        self.random_seed = random_seed
        np.random.seed(self.random_seed)
        random.seed(self.random_seed)

        # Create the XGBoost models using the hyperparameters stored in hyperparams.csv        
        self._create_models()

    def _create_models(self) -> None:
        """
        Creates the XGBoost models for each dataset.

        Parameters
        ----------

        None

        Returns
        -------

        None, saves the models in the self.models dictionary.
        """

        self.models = {}
        
        for i, row in self.hyperparams.iterrows():
            dataset = row['dataset']
            params = {
                'device': 'cuda',
                'nthread': 4,
                'objective': 'multi:softmax',
                'num_class': row['num_classes'],
                'max_depth': row['max_depth'],
                'learning_rate': row['learning_rate'],
                'gamma': row['gamma'],
                'eval_metric': ['mlogloss', 'auc']
            }
            self.models[dataset] = xgb.XGBClassifier(**params)        

    def _load_data(self, dataset: str, validation_size: float = 0.2, test_size: float = 0.1) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Loads the data for a single dataset.

        Parameters
        ----------

        dataset : str
            The name of the dataset, e.g. Cheung2018M.
        validation_size : float, optional, default = 0.2. The size of the validation set. Defaults to 0.2, i.e. 20% of the training (not total!) data is used for validation.
        test_size : float, optional, default = 0.1
            The size of the test set. Defaults to 0.1, i.e. 10% of the data is used for training and 90% for testing.

        Returns
        -------

        (expr_data_train, expr_data_val, expr_data_test, labels_data_train, labels_data_val, labels_data_test) : tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
            The expression data and labels for the training, validation and test sets.
        """

        expr_data = pd.read_csv(os.path.join(
            self.data_dir, f"{dataset}_expression.csv"), index_col=0)
        expr_data = expr_data.T

        labels_data = pd.read_csv(os.path.join(
            self.data_dir, f"{dataset}_clusters.csv"), index_col=0)
        labels_data = labels_data['Cluster'].values

        expr_data_train, expr_data_test, labels_data_train, labels_data_test = train_test_split(
            expr_data, labels_data, train_size=1.0-test_size, random_state=self.random_seed)
        expr_data_train, expr_data_val, labels_data_train, labels_data_val = train_test_split(
            expr_data_train, labels_data_train, test_size=validation_size, random_state=self.random_seed)

        if self.verbose:
            print(
                f"Loaded data for {dataset}.\nTrain: {len(expr_data_train)}\nValidation: {len(expr_data_val)}\nTest: {len(expr_data_test)}")
            
        return (expr_data_train, expr_data_val, expr_data_test,
                labels_data_train, labels_data_val, labels_data_test)

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
            if f.endswith('.json') and "_f1_" in f:
                os.remove(os.path.join(self.output_dir, f))

    def load_best_models(self) -> None:
        """
        Loads the best models from the output directory.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        all_models = [f for f in os.listdir(self.output_dir) if (
            f.endswith('.json') and "_f1_" in f)]

        # TODO change this
        studies = {m.split("_best_")[0] for m in all_models}
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
                os.path.join(self.output_dir, f"{s}_history.csv"))
            if self.verbose:
                print(f"Loaded best model for {s} from {best_model}")

    def _train_model(self,
                     dataset: str) -> xgb.Booster:
        """
        Trains a single XGBoost model.

        Parameters
        ----------
        dataset : str
            The name of the dataset. This will be used to save the model and training history.

        Returns
        -------
        xgb_model : xgb.Booster. The trained XGBoost model.
        """

        start_time = time.time()
        print(
            f"Training model {dataset} for {self.hyperparams.loc[self.hyperparams['dataset'] == dataset, 'n_boost_rounds'].values[0]} rounds...")

        # Note: this is deterministic because we set the random seed
        train_expr, val_expr, _, train_labels, val_labels, _ = self._load_data(
            dataset)
        hp = self._get_hyperparameters(dataset)
        dtrain = xgb.DMatrix(train_expr, label=train_labels)
        dval = xgb.DMatrix(val_expr, label=val_labels)

        params = {
            'device': 'cuda',
            'nthread': 4,
            'objective': 'multi:softmax',
            'num_class': hp['num_classes'],
            'max_depth': hp['max_depth'],
            'learning_rate': hp['learning_rate'],
            'gamma': hp['gamma'],
            'eval_metric': ['mlogloss', 'auc']
        }

        evallist = [(dtrain, 'train'), (dval, 'eval')]
        self.models[dataset] = xgb.train(
            params, dtrain, hp['n_boost_rounds'], evallist)

        # Save model as <dataset_name>_<date>-<time>_f1_<macro_f1_score>.json
        date = time.strftime("%Y%m%d-%H%M%S")
        f1 = f1_score(val_labels, self.models[dataset].predict(
            dval), average='macro')
        self.models[dataset].save_model(os.path.join(
            self.output_dir, f"{dataset}_{date}_f1_{f1}.json"))

        end_time = time.time()

        if (self.verbose):
            print(
                f"Training XGBoost model for {dataset} took {end_time - start_time} seconds.")

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

            self._train_model(dataset=dataset)

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

            metrics[row['dataset']] = self.evaluate_single_model(
                row['dataset'])

        return metrics

    def evaluate_single_model(self, dataset: str) -> dict:
        """
        Evaluates a single model on the test set.

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

        dtest = xgb.DMatrix(expr_test, label=labels_test)
        y_pred = self.models[dataset].predict(dtest)
        metrics = {}
        metrics['accuracy'] = np.mean(labels_test == y_pred)
        metrics['micro_f1_score'] = f1_score(
            labels_test, y_pred, average='micro')
        metrics['macro_f1_score'] = f1_score(
            labels_test, y_pred, average='macro')

        return metrics

    def _get_hyperparameters(self, dataset: str) -> dict:
        """
        Returns the model for a given dataset.

        Parameters
        ----------
        dataset: str
            The name of the dataset.

        Returns
        -------
        hyperparameters : dict
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
            learning_rate, max_depth, n_boost_rounds, gamma
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

        expr_train, expr_val, _, labels_train, labels_val, _ = self._load_data(
            dataset, validation_size=0.2, test_size=0.1)

        tune_res = pd.DataFrame()

        for v in values:
            params = {
                'device': 'cuda',
                'objective': 'multi:softmax',
                'num_class': hyper['num_classes'],
                'max_depth': hyper['max_depth'] if hyperparam != 'max_depth' else v,
                'learning_rate': hyper['learning_rate'] if hyperparam != 'learning_rate' else v,
                'gamma': hyper['gamma'] if hyperparam != 'gamma' else v,
                'eval_metric': ['mlogloss', 'auc']
            }

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
