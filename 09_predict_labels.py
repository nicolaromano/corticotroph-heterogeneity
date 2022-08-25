# Silence Tensorflow warnings
import os

from sklearn.metrics import confusion_matrix
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

import argparse
import seaborn as sns
import matplotlib.pyplot as plt
from tensorflow import keras
import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix

parser = argparse.ArgumentParser()
parser.add_argument('--ref_dataset', type=str, required=True, help='Reference dataset (model trained on this must exist.')
parser.add_argument('--target_dataset', type=str, required=True, help='Target dataset (labels will be predicted on this).')

args = parser.parse_args()

model = keras.models.load_model(f"best_models/{args.ref_dataset}_weights.h5")

filename_expr = f"expr_matrices/{args.target_dataset}.csv.gz"
filename_labels = f"expr_matrices/{args.target_dataset}_clusters.csv"

expr = pd.read_csv(filename_expr, index_col=0)        
expr = expr.T
real_clusters = pd.read_csv(filename_labels, index_col=0)

predicted_clusters = np.argmax(model.predict(expr), axis=1)

cm = confusion_matrix(y_true=real_clusters, y_pred=predicted_clusters, normalize="true")

sns.heatmap(cm)
plt.title(f"Transfer from {args.ref_dataset} to {args.target_dataset}")
plt.show()

print(cm)
