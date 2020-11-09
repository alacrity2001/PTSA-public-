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
from functools import partial
from multiprocessing import Pool
from version import __version__
import random as rn
from collections import Counter
import warnings
from numpy import percentile


from ..utils.utility import *


class hiForest(object):
    '''hybrid isolation forest method 
    '''

    def buildTree(self, i):
        ix = rn.sample(range(self.nobjs), self.sample)
        X_p = self.X[ix]
        return(hiTree(X_p, 0, self.limit))

    def __init__(self,X, ntrees,  sample_size, limit=None, nCore=1, contamination = 1):
        self.ntrees = ntrees
        self.X = X
        self.nobjs = len(X)
        self.sample = sample_size
        self.Trees = []
        self.limit = limit
        self.nCore=nCore
        self.X_in=None
        self.contamination = contamination

        if limit is None:
            self.limit = int(np.ceil(1.2*np.log2(self.sample)))
        self.c = c_factor(self.sample)
        self.Trees = []
        for i in range(self.ntrees):
            ix = rn.sample(range(self.nobjs), self.sample)
            X_p = X[ix]
            self.Trees.append(hiTree(X_p, 0, self.limit))
            
    def fit(self, X = None, a2 = 0.8, a1 = 0.4):
        """Calculate the average abnormality score of each point by their expected 
        number of path in a tree using the hybrid method. The algorithm also calcuates
        the score from the centrod.
        ----------
        X_in: input X, if None, used the dataset set up in the object
        a0, a1: parameters that weights the 3 scores. 

        Returns
        -------
        self : object
            Fitted estimator with decision score
        """
        if X == None:
            X = self.X
        Sc0 = np.zeros(X_train.shape[0])
        Sc1 = np.zeros(X_train.shape[0])
        Sc2 = np.zeros(X_train.shape[0])
        l = np.zeros(X_train.shape[0])
        for i in range(X_train.shape[0]):
            Sc0[i], l,  Sc1[i], Sc2[i] = hit.computeAggScore(X_train[i])

        score = a2 * (a1 * self.normalize(Sc0) + (1-a1)* self.normalize(Sc1)) + (1- a2) * self.normalize(Sc2)
        hiForest.raw_Sc0 = Sc0
        hiForest.raw_Sc1 = Sc1
        hiForest.raw_Sc2 = Sc2
        hiForest._decision_score = score
        return self
            
    
    def normalize(self, x):
        '''normalize data
        '''
        mini = np.min(x)
        maxi = np.max(x)
        if (maxi - mini) != 0:
            return (x - mini)/(maxi - mini)
        else:
            return np.zeros(len(x))

    def computeScore_paths(self, X_in = None):
        """Calculate the average abnormality score of each point by their expected 
        number of path in a tree 
        this is the standard Iforest method
        ----------
        X_in: input X, if None, used the dataset set up in the object

        Returns
        -------
        self : object
            Fitted estimator with decision score
        """
        if X_in is None:
            X_in = self.X
        S = np.zeros(len(X_in))
        for i in  range(len(X_in)):
            h_temp = 0
            for j in range(self.ntrees):
                h_temp += PathFactor(X_in[i],self.Trees[j]).path*1.0
            Eh = h_temp/self.ntrees
            S[i] = 2.0**(-Eh/self.c)
            self.path_scores_ = S
        return self

    def computeScore(self, i):
        """Calculate the abnormality score of individual point by their expected 
        number of path in a tree
        ----------
        i: integer, index of the point to be calculated in the dataset 

        Returns
        -------
        self : object
            Fitted estimator with decision score
        """
        h_temp = 0
        for j in range(self.ntrees):
            h_temp += PathFactor(self.X_in[i], self.Trees[j]).path * 1.0
        Eh = h_temp / self.ntrees
        return 2.0 ** (-Eh / self.c)
    
    def _process_decision_scores(self):
        """Internal function to calculate key attributes:
        - threshold_: used to decide the binary label
        - labels_: binary labels of training data
        Returns
        -------
        self
        """

        self.threshold_ = percentile(self._decision_score,
                                     100 * (1 - self.contamination))
        self.labels_ = (self._decision_score > self.threshold_).astype(
            'int').ravel()

        # calculate for predict_proba()

        self._mu = np.mean(self._decision_score)
        self._sigma = np.std(self._decision_score)

        return self


    def computeScore_pathsPool(self, X_in = None):
        pool = Pool(self.nCore)
        if X_in is None:
            X_in = self.X
        self.X_in=X_in
        print(np.shape(X_in), self.nCore)
        L=len(X_in)
        tab = list(range(L))
        print(tab)
        S=pool.map(self.computeScore, tab)
        return S

    def computeScore_paths_single(self, x):
        S = np.zeros(self.ntrees)
        for j in range(self.ntrees):
            path =  PathFactor(x,self.Trees[j]).path*1.0
            S[j] = 2.0**(-1.0*path/self.c)
        return S

    def computeScore_paths_single_with_labs(self, x):
        S = np.zeros(self.ntrees)
        labs=[]
        for j in range(self.ntrees):
            pf=PathFactor(x,self.Trees[j])
            path =  pf.path*1.0
            S[j] = 2.0**(-1.0*path/self.c)
            labs.append(pf.labs)
        return S, labs


    def computeAggScore(self, x):
        S = np.zeros(self.ntrees)
        labsCount = Counter([])
        ldist=[]
        ldist_a=[]
        for j in range(self.ntrees):
            pf=PathFactor(x,self.Trees[j])
            path =  pf.path*1.0
            S[j] = 2.0**(-1.0*path/self.c)
            labsCount=labsCount+pf.labs
            if(len(pf.ldist)>0):
                ldist.append(np.mean(pf.ldist))
            if(len(pf.ldist_a)>0):
                ldist_a.append(np.mean(pf.ldist_a, axis=0))
        meanDist=0
        if(len(ldist)>0):
            meanDist=np.mean(ldist)
        meanDist_r = 0
        if(len(ldist_a)>0):
            meanDist_a = np.mean(ldist_a, axis=0)
            if(meanDist_a>0):
                meanDist_r=meanDist/(meanDist_a)
                #meanDist_r = 1.0 / (meanDist_a)
        return np.mean(S), labsCount, meanDist, meanDist_r
    



    def addAnomaly(self, x, lab):
        for j in range(self.ntrees):
            pf=PathFactor(x,self.Trees[j])
            pf.addAnomaly(x, lab, self.Trees[j].root)

    def computeAnomalyCentroid(self):
        for j in range(self.ntrees):
            self.Trees[j].root.computeAnomalyCentroid()

    def getAverageBucketSize(self):
        szb=0
        nbb=0
        for j in range(self.ntrees):
            s,n=self.Trees[j].root.getAverageBucketSize()
            szb+=s
            nbb+=n
        out=0
        if nbb>0:
            out=szb/nbb
        return out, nbb/self.ntrees



class Node(object):
    def __init__(self, X, q, p, e, left, right, node_type = '' ):
        self.e = e
        self.size = len(X)
        self.X = X # to be removed
        self.q = q
        self.p = p
        self.left = left
        self.right = right
        self.ntype = node_type
        self.C = None
        self.Ca = None
        self.labs=[]
        self.Xanomaly=[]
        if(node_type == 'exNode' and self.size>0):
            self.C = np.mean(X, axis=0)


    def computeAnomalyCentroid(self):
        if self.ntype == 'exNode':
            if(len(self.Xanomaly)>0):
                self.Ca = np.mean(self.Xanomaly, axis=0)
        else:
            self.left.computeAnomalyCentroid()
            self.right.computeAnomalyCentroid()

    def getAverageBucketSize(self):
        if self.ntype == 'exNode':
            return self.size, 1
        else:
            s1, n1 = self.left.getAverageBucketSize()
            s2, n2 = self.right.getAverageBucketSize()
            return s1+s2, n1+n2

class hiTree(object):

    """
    Unique entries for X
    """

    def __init__(self,X,e,l):
        self.e = e # depth
        self.X = X #save data for now
        self.size = len(X) #  n objects
        self.Q = np.arange(np.shape(X)[1], dtype='int') # n dimensions
        self.l = l # depth limit
        self.p = None
        self.q = None
        self.exnodes = 0
        self.labs=[]
        self.root = self.make_tree(X,e,l)
        

    def make_tree(self,X,e,l):
        self.e = e
        if e >= l or len(X) <= 1:
            left = None
            right = None
            self.exnodes += 1
            return Node(X, self.q, self.p, e, left, right, node_type = 'exNode' )
        else:
            self.q = rn.choice(self.Q)
            self.p = rn.uniform(X[:,self.q].min(),X[:,self.q].max())
            w = np.where(X[:,self.q] < self.p,True,False)
            return Node(X, self.q, self.p, e,\
            left=self.make_tree(X[w],e+1,l),\
            right=self.make_tree(X[~w],e+1,l),\
            node_type = 'inNode' )

    def get_node(self, path):
        node = self.root
        for p in path:
            if p == 'L' : node = node.left
            if p == 'R' : node = node.right
        return node


class PathFactor(object):
    def __init__(self,x,hitree):
        self.path_list=[]
        self.labs = []
        self.ldist = []
        self.ldist_a = []
        self.x = x
        self.e = 0
        self.path = self.find_path(hitree.root)


    def find_path(self,T):
        if T.ntype == 'exNode':
            self.labs = Counter(T.labs)
            if not (T.C is None):
                self.ldist.append(EuclideanDist(self.x, T.C))
            if not (T.Ca is None):
                self.ldist_a.append(EuclideanDist(self.x, T.Ca))
            sz=T.size
            if(sz==0):
                sz+=1
            for key in self.labs:
                self.labs[key] /= sz
            if T.size == 1:
                return self.e
            else:
                self.e = self.e + c_factor(T.size)
                return self.e
        else:
            a = T.q
            self.e += 1
            if self.x[a] < T.p:
                self.path_list.append('L')
                return self.find_path(T.left)
            else:
                self.path_list.append('R')
                return self.find_path(T.right)

    def addAnomaly(self, x, lab, T):
        if T.ntype == 'exNode':
            T.labs.append(lab)
            T.Xanomaly.append(x)
        else:
            a = T.q
            if self.x[a] < T.p:
                return self.addAnomaly(x, lab, T.left)
            else:
                return self.addAnomaly(x, lab, T.right)