# -*- coding: utf-8 -*-
"""
Created on Thu May  2 16:59:02 2024

@author: armelsoubeiga
"""

from tslearn.datasets import UCR_UEA_datasets
import pandas as pd
import numpy as np
import os

data_loader = UCR_UEA_datasets()
data_loader.list_multivariate_datasets()
l = ['UWaveGestureLibrary','StandWalkJump']

os.chdir('F:/phd/1-Longitudinal/Multiview/Experimental/Data/time-series')
for nom_dataset in l:
    X_train, y_train, X_test, y_test = data_loader.load_dataset(nom_dataset)

    id_col = np.repeat(np.arange(X_train.shape[0]), X_train.shape[1])
    time_col = np.tile(np.arange(X_train.shape[1]), X_train.shape[0])
    label_col = np.repeat(y_train, X_train.shape[1])

    df = pd.DataFrame({
        'id': id_col,
        'time': time_col,
        'label': label_col,
    })

    for i in range(X_train.shape[2]):
        df[f'dimension_{i+1}'] = X_train[:, :, i].reshape(-1)

    df.to_excel(f"{nom_dataset}.xlsx", index=False)
