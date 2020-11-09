# -*- coding: utf-8 -*-
"""Yesterday method. Simple.
"""
# Author: Yinchen Wu <yinchen@uchicago.edu>

from __future__ import division
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import random
from arch import arch_model
import pandas as pd
import math
import pmdarima as pm
from pmdarima import model_selection
import os
import dis
import statistics
from sklearn import metrics
import sklearn
from ptsa.models.detectorA import DetectorA



import warnings
from collections import defaultdict

from ptsa.utils.utility import _get_sklearn_version

if _get_sklearn_version() > 20:
    from inspect import signature
else:
    from sklearn.externals.funcsigs import signature

import abc
import six

import numpy as np
from numpy import percentile
from scipy.special import erf
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import roc_auc_score
from sklearn.utils import deprecated
from sklearn.utils.validation import check_is_fitted
from sklearn.utils.multiclass import check_classification_targets
import scipy.stats

from .sklearn_base import _pprint
from ptsa.utils.utility import precision_n_scores


class Yesterday(DetectorA):
    """A straighforward method that uses past data to estimate
    ----------
    Day : int, optional(default = 288)
        Number of samples per day
    Period : int, optional(default = 1)
        Number of days between the past reference data and current data
    neighborhood : int, string 'all', optional (default = 'all')
        The number of samples to fit for one subsequence. The timeseries may have
        seaonal shifts, hence if neighborhood = N, it assumes the N nearest points around 
        the sequence thats being estimated come from same distribution. If N = 'all', it 
        treats it as a consistent process.
    window: int, optional (default = 20)
        The length of the window to detect the given anomolies 
    contamination : float in (0., 0.55), optional (default=0.1)
        The amount of contamination of the data set, i.e. the proportion
        of outliers in the data set. Used when fitting to define the threshold
        on the decision function.

    Attributes
    ----------
    estimators_ : The ARIMA class that's being generated at the end. 
    decision_scores_ : numpy array of shape (n_samples,)
        The outlier scores of the training data.
        The higher, the more abnormal. Outliers tend to have higher
        scores. This value is available once the detector is
        fitted.
    threshold_ : float
        The threshold is based on ``contamination``. It is the
        ``n_samples * contamination`` most abnormal samples in
        ``decision_scores_``. The threshold is calculated for generating
        binary outlier labels.
    labels_ : int, either 0 or 1
        The binary labels of the training data. 0 stands for inliers
        and 1 for outliers/anomalies. It is generated by applying
        ``threshold_`` on ``decision_scores_``.
    """
    def __init__(self,  window = 20,  day = 288, period =1, contamination = 0.1, neighborhood = 'all'):

        self.window = window
        self.day = day
        self.period = period
        self.contamination = contamination
        self.neighborhood = neighborhood



    def fit(self, X, y=None):
        """Fit detector. y is ignored in unsupervised methods.
        Parameters
        ----------
        X : numpy array of shape (n_samples, )
            The input samples.
        y : Ignored
            Not used, present for API consistency by convention.
        Returns
        -------
        self : object
            Fitted estimator.
        """
        # validate inputs X and y (optional)
        self._set_n_classes(y)

        self.X_train_ = X
        self.n_train_ = len(X)
        self.decision_scores_ = np.zeros([self.n_train_, 1])
        
        self.n_initial_ = min(500, int(0.1 * self.n_train_))
        self.X_initial_ = X[:self.n_initial_]                         
        
        day = self.day
        period = self.period
        window = self.window
        data = X
        estimation = X.copy()
        
        difference = day * period 
        
        start = max(difference, self.n_initial_)

        for i in range(start, self.n_train_):
            estimation[i] = X[i - difference]
            
        self.estimation = estimation
        self.estimator = None
        return self
        
        
    
    def decision_function(self, X= False, measure = None):
        """Derive the decision score based on the given distance measure
        Parameters
        ----------
        X : numpy array of shape (n_samples, )
            The input samples.
        measure : object
            object for given distance measure with methods to derive the score
        Returns
        -------
        self : object
            Fitted estimator.
        """
        neighborhood = self.neighborhood
        if neighborhood == 'all':
            neighborhood = self.n_train_
            self.neighborhood = neighborhood
        if type(X) != bool:
            self.X_train_ = X
        estimation = self.estimation
        window = self.window
        n_train_ = self.n_train_
        score = np.zeros(n_train_)
        
        measure.detector = self
        measure.set_param()
        i = 0
        while i + window + self.n_initial_ < (n_train_):
            index = i+self.n_initial_
            score[index: index+window] = measure.measure(self.X_train_[index:index+window], estimation[index:index+window], index)
            i += window
        if i + self.n_initial_ < n_train_:
            score[i + self.n_initial_: ] = measure.measure(self.X_train_[i + self.n_initial_:], estimation[i + self.n_initial_:], i + self.n_initial_)
        self.decision_scores_ = score
        return self

    def predict_proba(self, X, method='linear', measure = None):
        """Predict the probability of a sample being outlier. Two approaches
        are possible:
        1. simply use Min-max conversion to linearly transform the outlier
           scores into the range of [0,1]. The model must be
           fitted first.
        2. use unifying scores, see :cite:`kriegel2011interpreting`.
        Parameters
        ----------
        X : numpy array of shape (n_samples, n_features)
            The input samples.
        method : str, optional (default='linear')
            probability conversion method. It must be one of
            'linear' or 'unify'.
        Returns
        -------
        outlier_probability : numpy array of shape (n_samples,)
            For each observation, tells whether or not
            it should be considered as an outlier according to the
            fitted model. Return the outlier probability, ranging
            in [0,1].
        """

        check_is_fitted(self, ['decision_scores_', 'threshold_', 'labels_'])
        train_scores = self.decision_scores_

        self.fit(X)
        self.decision_function(measure = measure)
        test_scores = self.decision_scores_

        probs = np.zeros([X.shape[0], int(self._classes)])
        if method == 'linear':
            scaler = MinMaxScaler().fit(train_scores.reshape(-1, 1))
            probs[:, 1] = scaler.transform(
                test_scores.reshape(-1, 1)).ravel().clip(0, 1)
            probs[:, 0] = 1 - probs[:, 1]
            return probs

        elif method == 'unify':
            # turn output into probability
            pre_erf_score = (test_scores - self._mu) / (
                    self._sigma * np.sqrt(2))
            erf_score = erf(pre_erf_score)
            probs[:, 1] = erf_score.clip(0, 1).ravel()
            probs[:, 0] = 1 - probs[:, 1]
            return probs
        else:
            raise ValueError(method,
                             'is not a valid probability conversion method')