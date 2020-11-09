# -*- coding: utf-8 -*-
"""SSA method based on the paper "Online Anomaly Detection for Sensor Systems: a Simple and
Efficient Approach" 
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

class SSA:
    
    def __init__(self, ep=3, method='GARCH', rf_method = 'all', n=720, a=0.2):
        self.ep = ep
        self.n = n
        self.method = method
        self.rf_method = rf_method
        self.a = a

    def reference_time(self, index, a, ground_truth = False):
        """Obtain the reference timeseries.
        Parameters
        ----------
        Index : integer
            the current datapoint being tested.
        ground_truth : boolean 
            If True, return X as the timeseries
        a: float, integer, or numpy array 
        weights to obtain the 
        Returns
        -------
        X1: numpy array of shape (n, )
        the reference timeseries
        X2: numpy array of shape (n, )
        the test timeseries
        """
        n  = self.n
        X  = self.X_train_
        if ground_truth == True:
            
            return X[index - 2*n: index - n], X[index - n: index]
        else:
            if type(a) == float:
                return X[index - 2*n: index - n]*(1-a) + a*X[index - n: index] , X[index - n: index]
            else:
                num = math.floor(index/n)
                A = a[:num]
                A = A / np.sum(A)
                rf = np.zeros(n)
                for i in range(len(A)):
                    rf += A[i] * X[index - (i+1)*n: index - i * n]
                    return rf,  X[index - n: index]

    def Linearization(self, X2, e=1):
        """Obtain the linearized curve.
        Parameters
        ----------
        X2 : numpy array of shape (n, )
            the time series curve to be fitted
        e: float, integer, or numpy array 
        weights to obtain the 
        Returns
        -------
        fit: parameters for the fitted linear curve
        """
        i = 0
        fit = {}
        fit['index'] = []
        fit['rep'] = []
        while i < len(X2):
            fit['index'].append(i)
            fit['Y'+str(i)]= X2[i]
            fit['rep'].append(np.array([i, X2[i]]))
            if i+1 >= len(X2):
                    break
            k = X2[i+1]-X2[i]
            b = -i*(X2[i+1]-X2[i])+X2[i]
            fit['reg' +str(i)]= np.array([k, b])
            i += 2
            if i >= len(X2):
                break
            d = np.abs(X2[i]- (k * i+b))
            while d < e:
                i +=1 
                if i >= len(X2):
                    break
                d = np.abs(X2[i]- (k * i+b)) 
        return fit

    def SSA(self, X2, X3, e = 1):
        """Obtain the SSA similarity score.
        Parameters
        ----------
        X2 : numpy array of shape (n, )
            the reference timeseries
        X3 : numpy array of shape (n, )
            the tested timeseries
        e: float, integer, or numpy array 
        weights to obtain the 
        Returns
        -------
        score: float, the higher the more dissimilar are the two curves 
        """       
        #linearization of data X2 and X3
        fit = self.Linearization(X2, e=e)
        fit2 = self.Linearization(X3, e=e)
    
        #line alinement 
        Index = []
        test_list = fit['index'] + fit2['index']
        [Index.append(x) for x in test_list if x not in Index]
        Y = 0
    
        #Similarity Computation
        for i in Index:
            if i in fit['index'] and i in fit2['index']:
                Y += abs(fit['Y'+str(i)]-fit2['Y'+str(i)])

            elif i in fit['index']:
                J = np.max(np.where(np.array(fit2['index']) < i ))
                index = fit2['index'][J]
                k = fit2['reg'+str(index)][0]
                b = fit2['reg'+str(index)][1]
                value = abs(k * i + b - fit['Y'+str(i)])
                Y += value
            elif i in fit2['index']:
                J = np.max(np.where(np.array(fit['index']) < i ))
                index = fit['index'][J]
                k = fit['reg'+str(index)][0]
                b = fit['reg'+str(index)][1]
                value = abs(k * i + b - fit2['Y'+str(i)])
                Y += value
        score = Y/len(Index)
        return score     

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
        n = self.n
        self.decision_scores_ = np.zeros(len(X))
        self.raw_decision_scores_ = np.zeros(len(X))
        a = self.a
        ep = self.ep

        if self.rf_method == 'all': 
            rf = np.zeros(self.n)
            num = math.floor(self.n_train_ / self.n)
            for i in range(num):
                rf += X[i * self.n: (i+1) * self.n]/num
            for i in range(2*n, len(X)):
                X1 = X[i-n:i]
                score = self.SSA(X1, rf)
                self.decision_scores_[i] = min(score/(2*ep), 1)
                self.raw_decision_scores_[i] = score
            
        elif self.rf_method == 'alpha':
            for i in range(2*n, len(X)):
                X1, X2 = self.reference_time(index = i, a=a)
                score = 0
                score = self.SSA(X1, X2)
                self.decision_scores_[i] = min(score/(2*ep), 1)
                self.raw_decision_scores_[i] = score
                
        else:
            raise ValueError(self.method, "is not a valid reference timeseries method")

        # flip the scores
        #self.decision_scores_ = self.decision_scores_.ravel() * -1
        #self._process_decision_scores()
        return self

    def _set_n_classes(self, y):
        """Set the number of classes if `y` is presented for multi-class outlier detection. If 
        y is not presented, nummber of classes is set to 2.
        Parameters
        ----------
        y : numpy array of shape (n_samples,)
            Ground truth.
        Returns
        -------
        self
        """

        self._classes = 2  # default as binary classification
        if y is not None:
            sklearn.utils.multiclass.check_classification_targets(y)
            self._classes = len(np.unique(y))
            warnings.warn(
                "y should not be presented in unsupervised learning.")
        return self   
