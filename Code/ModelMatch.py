# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Last modified 16 July 2014 by KO & PBS
# 
# This notebook contains functions to implement model-based matching, including stratifying samples, computing test statistics, and conducting permutation tests.

# <codecell>

import math
import numpy as np
import scipy
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
from permute import *
import sklearn.ensemble

# <codecell>

def binByQuantiles(predicted, treatment, nbins = 4, verbosity = 0):
    ''' 
    Stratify observations by their predicted response values.
    Use the quantiles of the sample to return nbins groups, each containing
    approximately the same number of individuals.
    
    Inputs: 
    predicted = pandas Series containing model-predicted response
    treatment = pandas Series containing treatment values
    nbins = number of groups
    verbosity = flag for printing and debugging. 
        0 for no printed output, 1 for regular verbosity, 2 for verbosity when treatment is binary
    
    Dependencies: numpy (as np), scipy.stats, pandas as pd
    '''
    n = len(predicted)
    q = (np.arange(nbins+1))*100.0/nbins
    quantiles = scipy.stats.scoreatpercentile(predicted, q)
    groups = pd.Series(np.zeros(n))
    # potential problem - how to handle ties in the quantiles?
    for i in np.arange(nbins)+1:
        groups[predicted>quantiles[i]] += 1
    if verbosity>0:
        print 'Quantiles used for binning: ', quantiles
        print '\nNumber of observations assigned to each bin:'
        if verbosity == 2: 
            print pd.crosstab(groups,treatment,rownames=['Bin'],colnames=['Treatment Group'])
        else:
            print groups.value_counts()
    return groups

# <codecell>

def modelMatch(predicted, response, conditions, bin_method = "quantile", nbins = 4,
               testStatistic = "difference_in_means", iterations = 1000, verbosity = 0):
    '''
    This function is a wrapper for the binning and stratified permutation test methods.
    
    Inputs: 
    predicted = response predicted using all variables EXCEPT the treatment
    response = Series of measured responses, one entry per individual
    conditions = Series of treatment levels, one entry per individual
    bin_method = a string to specify binning method. For now, the only option is "quantile"
    nbins = number of bins, default 4
    testStatistic = a string to specify the test statistic to be used in the stratified permutation test.
        For now, the options are "difference_in_means", "pearson_r"
    iterations = number of iterations for simulating the permutation distribution
    verbosity = flag for printing and debugging. 
        0 for no printed output, 1 for some verbosity, 2 for maximum verbosity
        
    Outputs:
    pLeft, pRight, pBoth = permutation test p-values
    tst = test statistic
    dist = (empirical) permutation distribution of test statistic
    '''
    binning = dict( \
             quantile = lambda u: binByQuantiles(u, treatment = conditions, nbins = nbins, verbosity = verbosity)
            )
    try:
        stratify = binning[bin_method]
    except KeyError:
        raise ValueError("Unrecognized binning method: " + bin_method) 
    
    stats = dict( \
             difference_in_means = lambda u: stratifiedPermutationTest(u, conditions, response, iterations=iterations, testStatistic=stratifiedPermutationTestMean, verbosity=verbosity),
             pearson_r = lambda u: stratifiedPermutationTest(u, conditions, response-predicted, iterations=iterations, testStatistic=stratifiedPermutationTestPearsonr, verbosity=verbosity) \
             )
    try:
        tst = stats[testStatistic]
    except KeyError:
        raise ValueError("Unrecognized test statistic: " + testStatistic)  
        
    groups = stratify(predicted)
    pLeft, pRight, pBoth, tst, dist = tst(groups)
    return pLeft, pRight, pBoth, tst, dist

# <markdowncell>




