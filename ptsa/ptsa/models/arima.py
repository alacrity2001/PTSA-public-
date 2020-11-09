# -*- coding: utf-8 -*-
"""Anomoly detector using ARIMA estimation + GARCH value  
"""
# Author: Yinchen Wu <yinchen@uchicago.edu>

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
#from ptsa.models.poly import Poly
#from ptsa.models.cheb import Cheb
import sklearn

from ptsa.models.detectorA import DetectorA
from ..utils.utility import precision_n_scores

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
#from ptsa.models.poly import Poly
#from ptsa.models.cheb import Cheb
import sklearn

from ptsa.models.detectorA import DetectorA
from ptsa.utils.utility import precision_n_scores

class ARIMA(DetectorA):
    """An elementary method to detect anomolies using ARIMA approxiamtion. 
    A ARIMA process of certain order is fitted to the given timeseries dataset.
    A anomoly score is derived on each point based on the dissimiarlity measure used.
    Parameters
    ----------
     p_start, q_start,  max_p, max_q, d : int, optional (default=1, 1, 5, 5, 0)
        Parameters to initialize the auto arima process for the initial data. 
        It will pick the most likely combination to fit the data. 
    neighborhood : int, string 'all', optional (default = 'all')
        The number of samples to fit for one subsequence. The timeseries may have
        seaonal shifts, hence if neighborhood = N, it assumes the N nearest points around 
        the sequence thats being estimated come from same distribution. If N = 'all', it 
        treats it as a consistent process.
    window: int, optional (default = 20)
        The length of the window to detect the given anomalies 
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
    def __init__(self, clean = False, window = 20, max_lag = 30000,  p_start=1, q_start=1,  max_p=5, max_q = 5, d= 0, contamination = 0.1, neighborhood = 'all'):

        self.window = window
        self.p_start = p_start
        self.max_p = max_p
        self.max_q = max_q
        self.d = d
        self.q_start = q_start
        self.contamination = contamination
        self.neighborhood = neighborhood
        self.clean = clean
        self.max_lag = max_lag



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
        max_lag = self.max_lag
        self.X_train_ = X
        self.n_train_ = len(X)
        self.decision_scores_ = np.zeros([self.n_train_, 1])
        
        self.n_initial_ = min(500, int(0.1 * self.n_train_))
        self.X_initial_ = X[:self.n_initial_]    
        
        if max_lag < self.n_initial_:
            max_lag = int(self.n_initial_ * 1.1)

        window = self.window
        data = X
        
        #clean data
        clean = self.clean
        if clean:
            for i in range(self.n_initial_):
                k = 1
                while data[i] == 0:
                    data[i] = data[i + k]
                    k += 1

        p_start = self.p_start
        q_start = self.q_start
        max_p = self.max_p
        max_q = self.max_q
        d = self.d
        
        train, test = model_selection.train_test_split(data, train_size=self.n_initial_)
        arima = pm.auto_arima(train, start_p=p_start, start_q=q_start, d=d, max_p=max_p, max_q=max_q,
                      out_of_sample_size=10, suppress_warnings=True,
                      stepwise=True, error_action='ignore')
        
        self.start_arima = arima
        estimation = np.zeros(self.n_train_)

        i = 0
        j = 0
        while i + window < (len(test)):
            if i % (10 * window) == 0:
                print('i',i)
                print('j',j)
                
            #restart process
            if j > max_lag:
                j = 0
                train = data[i - self.n_initial_:i]
                arima = pm.auto_arima(train, start_p=p_start, start_q=q_start, d=d, max_p=max_p, max_q=max_q,
                      out_of_sample_size=10, suppress_warnings=True,
                      stepwise=True, error_action='ignore')
            pred = arima.predict(n_periods=window, return_conf_int=False)
        
            estimation[i+self.n_initial_: i + self.n_initial_ + window] = pred
            arima.update(test[i: i + window])

            i += window
            j += window

        pred = arima.predict(n_periods=(len(test) - i), return_conf_int=False)
        estimation[i+self.n_initial_: ] = pred
        estimation[:self.n_initial_] = train

        self.estimation = estimation
        self.estimator = arima
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
            if i % (10 * window) == 0:
                print(i)
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