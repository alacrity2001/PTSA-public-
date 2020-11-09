
import keras 
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import numpy as np
import time
import pandas as pd
import os

class LSTM:
#LSTM based
    def __init__(self, window = 30, support = 200, contamination = 0.1, epochs = 10):
        self.support = support
        self.window = window   
        self.contamination = contamination
        self.epochs = epochs

    def fit(self, X, y=None, ratio = 0.15):
        # validate inputs X and y (optional)
        window = self.window
        support = self.support
        TIME_STEPS =  self.window 
        epochs = self.epochs
        self.X_train_ = X
        self.n_train_ = len(X)
        
        x_train = X[:int(0.15 * len(X))]
        x_test = X[int(0.15 * len(data)):]
        if y is None:
            y_train = np.zeros(x_train.shape)
        else:
            y_train = y[:int(0.15 * len(X))]


        X_train, Y_train = self.create_dataset(
          x_train,
          y_train,
          TIME_STEPS,
          support
        )
        
        X_test, Y_test = self.create_dataset(
          x_test,
          np.zeros(x_test.shape),
          TIME_STEPS,
          support
        )

        model = keras.Sequential()
        model.add(keras.layers.LSTM(
            units=64,
            input_shape=(X_train.shape[1], X_train.shape[2])))
                  
        model.add(keras.layers.Dense(window))
        model.compile(loss='mean_squared_error', optimizer='adam')
        model.fit(X_train, Y_train, epochs=10, batch_size=32, verbose=2)
        
        

        prediction = model.predict(X_test)
        self.X_test = X_test
        self.Y_test = Y_test
        self.estimation = prediction
        self.estimator = model
        self.n_initial = X_train.shape[0]
        
        return self
        
    def create_dataset(self, X, y, time_steps=1, support = 10):
        Xs, ys = [], []
        for i in range(len(X) - time_steps - support):
            v = X[i:(i + support)]
            w = X[i+support : i+time_steps + support]
            Xs.append(v)
            ys.append(w)
            
            
        return np.array(Xs), np.array(ys)
    
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
        if type(X) != bool:
            self.X_train_ = X
        n_train_ = self.n_train_
        self.neighborhood = n_train_
        autoencoder = self.estimator
        Y_test = self.Y_test
        window = self.window
        score = np.zeros(n_train_)
        measure.detector = self
        measure.set_param()
        estimation = self.estimation
        n_initial = self.n_initial

        i = 0
        m = 0
        for i in range(estimation.shape[0]):
            

            score[i - estimation.shape[0]] = measure.measure(Y_test[i], estimation[i], n_train_ - estimation.shape[0] + i)
        
        print('done!')
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


class LSTM_Auto:
#Autoencoder using LSTM    
    def __init__(self, window = 30, latent_space = 64, contamination = 0.1, epochs = 10):
        self.latent_space = latent_space
        self.window = window   
        self.contamination = contamination
        self.epochs = epochs

    def fit(self, X, y=None, ratio = 0.15):
        # validate inputs X and y (optional)
        units = self.latent_space
        TIME_STEPS =  self.window 
        epochs = self.epochs
        self.X_train_ = X
        self.n_train_ = len(X)
        
        x_train = X[:int(0.15 * len(X))]
        x_test = X[int(0.15 * len(data)):]
        if y is None:
            y_train = np.zeros(x_train.shape)
        else:
            y_train = y[:int(0.15 * len(X))]


        X_train, Y_train = self.create_dataset(
          x_train,
          y_train,
          TIME_STEPS
        )
        
        X_test, Y_test = self.create_dataset(
          x_test,
          np.zeros(x_test.shape),
          TIME_STEPS
        )
        
        model = keras.Sequential()
        model.add(keras.layers.LSTM(
            units=units,
            input_shape=(X_train.shape[1], X_train.shape[2])
        ))
        model.add(keras.layers.Dropout(rate=0.2))
        model.add(keras.layers.RepeatVector(n=X_train.shape[1]))
        model.add(keras.layers.LSTM(units=units, return_sequences=True))
        model.add(keras.layers.Dropout(rate=0.2))
        model.add(
          keras.layers.TimeDistributed(
            keras.layers.Dense(units=X_train.shape[2])
          )
        )
        model.compile(loss='mae', optimizer='adam')
        
        history = model.fit(
        X_train, Y_train,
        epochs=epochs,
        batch_size=32,
        validation_split=0.1,
        shuffle=False
    )

        prediction = model.predict(X_test)
        self.X_test = X_test
        self.estimation = prediction
        self.estimator = model
        self.n_initial = X_train.shape[0]
        
        return self
        
    def create_dataset(self, X, y, time_steps=1):
        Xs, ys = [], []
        for i in range(len(X) - time_steps):
            v = X[i:(i + time_steps)]
            Xs.append(v)
            ys.append(y[i + time_steps])
            
            
        return np.array(Xs), np.array(ys)
    
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
        if type(X) != bool:
            self.X_train_ = X
        n_train_ = self.n_train_
        self.neighborhood = n_train_
        autoencoder = self.estimator
        X_test = self.X_test
        window = self.window
        score = np.zeros(n_train_)
        measure.detector = self
        measure.set_param()
        estimation = self.estimation
        n_initial = self.n_initial

        i = 0
        m = 0
        for i in range(estimation.shape[0]):
            

            score[i + n_initial + window] = measure.measure(X_test[i], estimation[i], i + n_initial + window)
        
        print('done!')
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

class ANN:
#ANN based
    def __init__(self, window = 30, support = 200, contamination = 0.1, epochs = 10):
        self.support = support
        self.window = window   
        self.contamination = contamination
        self.epochs = epochs

    def fit(self, X, y=None, ratio = 0.15):
        # validate inputs X and y (optional)
        window = self.window
        support = self.support
        TIME_STEPS =  self.window 
        epochs = self.epochs
        self.X_train_ = X
        self.n_train_ = len(X)
        
        x_train = X[:int(0.15 * len(X))]
        x_test = X[int(0.15 * len(data)):]
        if y is None:
            y_train = np.zeros(x_train.shape)
        else:
            y_train = y[:int(0.15 * len(X))]


        X_train, Y_train = self.create_dataset(
          x_train,
          y_train,
          TIME_STEPS,
          support
        )
        
        X_test, Y_test = self.create_dataset(
          x_test,
          np.zeros(x_test.shape),
          TIME_STEPS,
          support
        )

        model = keras.Sequential()
        model.add(keras.layers.Dense(
            units=100,
            activation='relu', 
            input_dim=support))
                  
        model.add(keras.layers.Dense(window))
        model.compile(loss='mean_squared_error', optimizer='adam')
        

        model.fit(X_train, Y_train, epochs=epochs, batch_size=32, verbose=2)

        

        prediction = model.predict(X_test)
        self.X_test = X_test
        self.Y_test = Y_test
        self.estimation = prediction
        self.estimator = model
        self.n_initial = X_train.shape[0]
        
        return self
        
    def create_dataset(self, X, y, time_steps=1, support = 10):
        Xs, ys = [], []
        for i in range(len(X) - time_steps - support):
            v = X[i:(i + support)]
            w = X[i+support : i+time_steps + support]
            Xs.append(v)
            ys.append(w)
        x = np.array(Xs)
        y = np.array(ys)
        x = x.reshape((x.shape[0], x.shape[1]))
        y = y.reshape((y.shape[0], y.shape[1]))
        return x, y
    
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
        if type(X) != bool:
            self.X_train_ = X
        n_train_ = self.n_train_
        self.neighborhood = n_train_
        autoencoder = self.estimator
        Y_test = self.Y_test
        window = self.window
        score = np.zeros(n_train_)
        measure.detector = self
        measure.set_param()
        estimation = self.estimation
        n_initial = self.n_initial

        i = 0
        m = 0
        for i in range(estimation.shape[0]):
            

            score[i - estimation.shape[0]] = measure.measure(Y_test[i], estimation[i], n_train_ - estimation.shape[0] + i)
        
        print('done!')
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