Below is one of the projects I did in my RA at Uchicago computer science department that is open to public. Many others entail private data and private 
algirithms that are build for paper publishing. They are not made public until the paper is reviewed and published in the journal. 

Python TimeSeries Anomaly Detection (PTSA)

PTSA(in deveopment) is a comprehensive and scalable **Python toolkit** for **detecting anomalies** in
timeseries data. I already included 20 algorithms in the package, and more will come as I keep developing the package. 
There is a python notebook demo about how to use the package. 




**API Demo**\ :


.. code-block:: python


    # train the Polynomial detector
   from ptsa.models.poly import Poly
   model = Poly(window = 200, power = 0)
   model.fit(data)

    # get outlier scores
   from ptsa.models.distance import Fourier
   measure = Fourier()
   measure.detector = model
   measure.set_param()
   model.decision_function(measure=measure)
   
----


Installation
^^^^^^^^^^^^

It is recommended to use **pip** for installation. The required package is in requirement.txt. 

.. code-block:: bash

   git clone https://github.com/alacrity2001/AnomalyDetection
   cd code/ptsa
   pip install .
   
----
   
 **Required Dependencies**\ :

* Python 3.5, 3.6, or 3.7
* combo>=0.0.8
* joblib
* numpy>=1.13
* numba>=0.35
* pandas>=0.25
* scipy>=0.19.1
* scikit_learn>=0.19.1
* six
* statsmodels
* suod


**Optional Dependencies (see details below)**\ :


* keras (optional, required for AutoEncoder, ANN, LSTM)
* matplotlib (optional, required for running examples)
* pandas (optional, required for reading demo)
* tensorflow (optional, required for AutoEncoder, other backend works)
* pmdarima (optional, required for ARIMA)
* arch (optional, required for ARCH measure)
* tsfresh (optional, required for tf feature extraction)
* hurst (optional, tf feature extraction)

**(Note) Individual Detection Algorithms** :

The ptsa/models has all the scripts. 

The machine learning models for iforest, cblof, cof, knn, loci, lof, mcd, pca, ocsvm are which I modifed from PYOD libaray.  
The RPIF.py is hybrid isolation forest which is in https://github.com/pfmarteau/HIF 

I coded all the others includind the feature transformation, dissimilarity measure, etc. 
All the dissimilarity measure is quite standard, except SSA, which I referenced from the paper  http://bourbon.usc.edu/iml/projects/nsf-dc-0917340/anomaly.pdf

utils.utility are some functions I gathered from sklearn library and other authors of other alogirthms. 

ptsa.others.normats.norma is the norma algorithm whose author didn't pusblish his paper yet. The code is directly copied from
their folder with some debuggings for some errors encountered when applying our data. 

ptsa.utils.metrics has the class metricor I coded as a convenient tool to evaluate algorithms. The script I coded to compute 
ranged recall/precision is referenced from the paper https://arxiv.org/pdf/1803.03639.pdf 
