# -*- coding: utf-8 -*-
"""Clustering Based Local Outlier Factor (CBLOF)
"""
# Author: Yue Zhao <zhaoy@cmu.edu>
#         Shangwen Huang <https://github.com/shangwen777>
# License: BSD 2 clause

from __future__ import division
from __future__ import print_function

import warnings
import numpy as np
from scipy.spatial.distance import cdist
from sklearn.cluster import KMeans
from sklearn.utils.validation import check_is_fitted
from sklearn.utils.validation import check_array
from sklearn.utils.estimator_checks import check_estimator

from .detectorB import DetectorB
from ..utils.utility import check_parameter
from ..utils.stat_models import pairwise_distances_no_broadcast

__all__ = ['CBLOF']


class CBLOF(DetectorB):
    r"""The CBLOF operator calculates the outlier score based on cluster-based
    local outlier factor.
    CBLOF takes as an input the data set and the cluster model that was
    generated by a clustering algorithm. It classifies the clusters into small
    clusters and large clusters using the parameters alpha and beta.
    The anomaly score is then calculated based on the size of the cluster the
    point belongs to as well as the distance to the nearest large cluster.
    Use weighting for outlier factor based on the sizes of the clusters as
    proposed in the original publication. Since this might lead to unexpected
    behavior (outliers close to small clusters are not found), it is disabled
    by default.Outliers scores are solely computed based on their distance to
    the closest large cluster center.
    By default, kMeans is used for clustering algorithm instead of
    Squeezer algorithm mentioned in the original paper for multiple reasons.
    See :cite:`he2003discovering` for details.
    Parameters
    ----------
    n_clusters : int, optional (default=8)
        The number of clusters to form as well as the number of
        centroids to generate.
    contamination : float in (0., 0.5), optional (default=0.1)
        The amount of contamination of the data set,
        i.e. the proportion of outliers in the data set. Used when fitting to
        define the threshold on the decision function.
    clustering_estimator : Estimator, optional (default=None)
        The base clustering algorithm for performing data clustering.
        A valid clustering algorithm should be passed in. The estimator should
        have standard sklearn APIs, fit() and predict(). The estimator should
        have attributes ``labels_`` and ``cluster_centers_``.
        If ``cluster_centers_`` is not in the attributes once the model is fit,
        it is calculated as the mean of the samples in a cluster.
        If not set, CBLOF uses KMeans for scalability. See
        https://scikit-learn.org/stable/modules/generated/sklearn.cluster.KMeans.html
    alpha : float in (0.5, 1), optional (default=0.9)
        Coefficient for deciding small and large clusters. The ratio
        of the number of samples in large clusters to the number of samples in
        small clusters.
    beta : int or float in (1,), optional (default=5).
        Coefficient for deciding small and large clusters. For a list
        sorted clusters by size `|C1|, \|C2|, ..., |Cn|, beta = |Ck|/|Ck-1|`
    use_weights : bool, optional (default=False)
        If set to True, the size of clusters are used as weights in
        outlier score calculation.
    check_estimator : bool, optional (default=False)
        If set to True, check whether the base estimator is consistent with
        sklearn standard.
        .. warning::
            check_estimator may throw errors with scikit-learn 0.20 above.
    random_state : int, RandomState or None, optional (default=None)
        If int, random_state is the seed used by the random
        number generator; If RandomState instance, random_state is the random
        number generator; If None, the random number generator is the
        RandomState instance used by `np.random`.
    n_jobs : integer, optional (default=1)
        The number of jobs to run in parallel for both `fit` and `predict`.
        If -1, then the number of jobs is set to the number of cores.
    Attributes
    ----------
    clustering_estimator_ : Estimator, sklearn instance
        Base estimator for clustering.
    cluster_labels_ : list of shape (n_samples,)
        Cluster assignment for the training samples.
    n_clusters_ : int
        Actual number of clusters (possibly different from n_clusters).
    cluster_sizes_ : list of shape (n_clusters_,)
        The size of each cluster once fitted with the training data.
    decision_scores_ : numpy array of shape (n_samples,)
        The outlier scores of the training data.
        The higher, the more abnormal. Outliers tend to have higher scores.
        This value is available once the detector is fitted.
    cluster_centers_ : numpy array of shape (n_clusters_, n_features)
        The center of each cluster.
    small_cluster_labels_ : list of clusters numbers
        The cluster assignments belonging to small clusters.
    large_cluster_labels_ : list of clusters numbers
        The cluster assignments belonging to large clusters.
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

    def __init__(self, n_clusters=8, contamination=0.1,
                 clustering_estimator=None, alpha=0.9, beta=5,
                 use_weights=False, check_estimator=False, random_state=None,
                 n_jobs=1):
        super(CBLOF, self).__init__(contamination=contamination)
        self.n_clusters = n_clusters
        self.clustering_estimator = clustering_estimator
        self.alpha = alpha
        self.beta = beta
        self.use_weights = use_weights
        self.check_estimator = check_estimator
        self.random_state = random_state
        self.n_jobs = n_jobs

    # noinspection PyIncorrectDocstring
    def fit(self, X, y=None):
        """Fit detector. y is ignored in unsupervised methods.
        Parameters
        ----------
        X : numpy array of shape (n_samples, n_features)
            The input samples.
        y : Ignored
            Not used, present for API consistency by convention.
        Returns
        -------
        self : object
            Fitted estimator.
        """

        # validate inputs X and y (optional)
        X = check_array(X)
        self._set_n_classes(y)
        n_samples, n_features = X.shape

        # check parameters
        # number of clusters are default to 8
        self._validate_estimator(default=KMeans(
            n_clusters=self.n_clusters,
            random_state=self.random_state,
            n_jobs=self.n_jobs))

        self.clustering_estimator_.fit(X=X, y=y)
        # Get the labels of the clustering results
        # labels_ is consistent across sklearn clustering algorithms
        self.cluster_labels_ = self.clustering_estimator_.labels_
        self.cluster_sizes_ = np.bincount(self.cluster_labels_)

        # Get the actual number of clusters
        self.n_clusters_ = self.cluster_sizes_.shape[0]

        if self.n_clusters_ != self.n_clusters:
            warnings.warn("The chosen clustering for CBLOF forms {0} clusters"
                          "which is inconsistent with n_clusters ({1}).".
                          format(self.n_clusters_, self.n_clusters))

        self._set_cluster_centers(X, n_features)
        self._set_small_large_clusters(n_samples)

        self.decision_scores_ = self._decision_function(X,
                                                        self.cluster_labels_)

        self._process_decision_scores()
        return self

    def decision_function(self, X):
        """Predict raw anomaly score of X using the fitted detector.
        The anomaly score of an input sample is computed based on different
        detector algorithms. For consistency, outliers are assigned with
        larger anomaly scores.
        Parameters
        ----------
        X : numpy array of shape (n_samples, n_features)
            The training input samples. Sparse matrices are accepted only
            if they are supported by the base estimator.
        Returns
        -------
        anomaly_scores : numpy array of shape (n_samples,)
            The anomaly score of the input samples.
        """
        check_is_fitted(self, ['decision_scores_', 'threshold_', 'labels_'])
        X = check_array(X)
        labels = self.clustering_estimator_.predict(X)
        return self._decision_function(X, labels)

    def _validate_estimator(self, default=None):
        """Check the value of alpha and beta and clustering algorithm.
        """
        check_parameter(self.alpha, low=0, high=1, param_name='alpha',
                        include_left=False, include_right=False)

        check_parameter(self.beta, low=1, param_name='beta',
                        include_left=False)

        if self.clustering_estimator is not None:
            self.clustering_estimator_ = self.clustering_estimator
        else:
            self.clustering_estimator_ = default

        # make sure the base clustering algorithm is valid
        if self.clustering_estimator_ is None:
            raise ValueError("clustering algorithm cannot be None")

        if self.check_estimator:
            check_estimator(self.clustering_estimator_)

    def _set_cluster_centers(self, X, n_features):
        # Noted not all clustering algorithms have cluster_centers_
        if hasattr(self.clustering_estimator_, 'cluster_centers_'):
            self.cluster_centers_ = self.clustering_estimator_.cluster_centers_
        else:
            # Set the cluster center as the mean of all the samples within
            # the cluster
            warnings.warn("The chosen clustering for CBLOF does not have"
                          "the center of clusters. Calculate the center"
                          "as the mean of the clusters.")
            self.cluster_centers_ = np.zeros([self.n_clusters_, n_features])
            for i in range(self.n_clusters_):
                self.cluster_centers_[i, :] = np.mean(
                    X[np.where(self.cluster_labels_ == i)], axis=0)

    def _set_small_large_clusters(self, n_samples):
        # Sort the index of clusters by the number of samples belonging to it
        size_clusters = np.bincount(self.cluster_labels_)

        # Sort the order from the largest to the smallest
        sorted_cluster_indices = np.argsort(size_clusters * -1)

        # Initialize the lists of index that fulfill the requirements by
        # either alpha or beta
        alpha_list = []
        beta_list = []

        for i in range(1, self.n_clusters_):
            temp_sum = np.sum(size_clusters[sorted_cluster_indices[:i]])
            if temp_sum >= n_samples * self.alpha:
                alpha_list.append(i)

            if size_clusters[sorted_cluster_indices[i - 1]] / size_clusters[
                sorted_cluster_indices[i]] >= self.beta:
                beta_list.append(i)

            # Find the separation index fulfills both alpha and beta
        intersection = np.intersect1d(alpha_list, beta_list)

        if len(intersection) > 0:
            self._clustering_threshold = intersection[0]
        elif len(alpha_list) > 0:
            self._clustering_threshold = alpha_list[0]
        elif len(beta_list) > 0:
            self._clustering_threshold = beta_list[0]
        else:
            raise ValueError("Could not form valid cluster separation. Please "
                             "change n_clusters or change clustering method")

        self.small_cluster_labels_ = sorted_cluster_indices[
                                     self._clustering_threshold:]
        self.large_cluster_labels_ = sorted_cluster_indices[
                                     0:self._clustering_threshold]

        # No need to calculate small cluster center
        # self.small_cluster_centers_ = self.cluster_centers_[
        #     self.small_cluster_labels_]

        self._large_cluster_centers = self.cluster_centers_[
            self.large_cluster_labels_]

    def _decision_function(self, X, labels):
        # Initialize the score array
        scores = np.zeros([X.shape[0], ])

        small_indices = np.where(
            np.isin(labels, self.small_cluster_labels_))[0]
        large_indices = np.where(
            np.isin(labels, self.large_cluster_labels_))[0]

        if small_indices.shape[0] != 0:
            # Calculate the outlier factor for the samples in small clusters
            dist_to_large_center = cdist(X[small_indices, :],
                                         self._large_cluster_centers)

            scores[small_indices] = np.min(dist_to_large_center, axis=1)

        if large_indices.shape[0] != 0:
            # Calculate the outlier factor for the samples in large clusters
            large_centers = self.cluster_centers_[labels[large_indices]]

            scores[large_indices] = pairwise_distances_no_broadcast(
                X[large_indices, :], large_centers)

        if self.use_weights:
            # Weights are calculated as the number of elements in the cluster
            scores = scores * self.cluster_sizes_[labels]

        return scores.ravel()