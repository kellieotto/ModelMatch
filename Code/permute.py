# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Last modified 21 July 2014 by KO
# 
# This notebook implements a variety of permutation tests, including stratified permutation tests.
# It also implements exact confidence intervals for binomial p and hypergeometric parameters, by inverting tests.
# 
# <hr />

# <codecell>

import math
import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt

# <codecell>

def bisect(lo, hi, tol, fun):
    mid = (lo+hi)/2.0
    while (hi-lo)/2 > tol:
        if fun(mid) == 0.0:
            return mid
        elif fun(lo)*fun(mid) < 0.0:
            hi = mid
        else:
            lo = mid
        mid = (lo+hi)/2.0
    return mid
        
def binoLowerCL(n, x, cl = 0.975, inc=0.000001, p = None):
    "Lower confidence level cl confidence interval for Binomial p, for x successes in n trials"
    if p is None:
            p = float(x)/float(n)
    lo = 0.0
    if (x > 0):
            f = lambda q: cl - scipy.stats.binom.cdf(x-1, n, q)
            lo = bisect(0.0, p, inc, f)
    return lo

def binoUpperCL(n, x, cl = 0.975, inc=0.000001, p = None):
    "Upper confidence level cl confidence interval for Binomial p, for x successes in n trials"
    if p is None:
            p = float(x)/float(n)
    hi = 1.0
    if (x < n):
            f = lambda q: scipy.stats.binom.cdf(x, n, q) - (1-cl)
            hi = bisect(p, 1.0, inc, f)
    return hi

def permuTestMean(x, y, reps = 10**5, stat = 'mean', side = 'greater', CI =  False, CL = 0.95):
    """
       One-sided or two-sided, two-sample permutation test for equality of two 
       means, with p-value estimated by simulated random sampling with reps replications.
       
       Tests the hypothesis that x and y are a random partition of x,y
       against the alternative that x comes from a population with mean
           (a) greater than that of the population from which y comes, if side = 'greater_than'
           (b) less than that of the population from which y comes, if side = 'less_than'
           (c) different from that of the population from which y comes, if side = 'both'
       
       If stat == 'mean', the test statistic is (mean(x) - mean(y))
       (equivalently, sum(x), since those are monotonically related)
       
       If stat == 't', the test statistic is the two-sample t-statistic--but the p-value 
       is still estimated by the randomization, approximating the permutation distribution.
       The t-statistic is computed using scipy.stats.ttest_ind
       
       If CI == 'upper', computes an upper confidence bound on the true
       p-value based on the simulations by inverting Binomial tests.
       
       If CI == 'lower', computes a lower confidence bound on the true
       p-value based on the simulations by inverting Binomial tests.
       
       If CI == 'both', computes lower and upper confidence bounds on the true
       p-value based on the simulations by inverting Binomial tests.
       
       CL is the confidence limit for the confidence bounds.
       
       output is the estimated p-value and the test statistic, if CI == False
       output is <estimated p-value, confidence bound on p-value, test statistic> if CI in {'lower','upper'}
       output is <estimated p-value, [lower confidence bound, upper confidence bound], test statistic> if CI == 'both'
       
       Dependencies: numpy, numpy.random, scipy.stats, binoUpperCL, binoLowerCL
       
    """
    z  = np.concatenate([x, y])   # pooled responses
    stats = dict( \
             mean = lambda u: np.mean(u[:len(x)])-np.mean(u[len(x):]),
             t = lambda u: scipy.stats.ttest_ind(u[:len(y)], u[len(y):], equal_var=True)[0] \
            )
    try:
        tst = stats[stat]
    except KeyError:
        raise ValueError("Unrecognized test statistic (stat): " + stat)    
    if side == 'greater_than':
        theStat = tst
    elif side == 'less_than':
        theStat = lambda u: -tst(u)
    elif side == 'both':
        theStat = lambda u: math.fabs(tst(u))
    else:
        raise ValueError("Unrecognized side choice: " + side)
    ts = theStat(z)
    hits = np.sum([ (theStat(np.random.permutation(z)) >= ts) for i in range(reps)])
    if CI == 'upper':
        return float(hits)/float(reps), binoUpperCL(reps, hits, cl = CL), ts
    elif CI == 'lower':
        return float(hits)/float(reps), binoLowerCL(reps, hits, cl = CL), ts
    elif CI == 'both':
        return float(hits)/float(reps),  \
                 (binoLowerCL(reps, hits, cl = 1-(1-CL)/2), binoUpperCL(reps, hits, cl = 1-(1-CL)/2)), \
                 ts
    else:
        return float(hits)/float(reps), ts

    z  = np.concatenate([x, y])   # pooled responses
    stats = dict( \
             mean = lambda u: np.mean(u[:len(x)])-np.mean(u[len(x):]),
             t = lambda u: scipy.stats.ttest_ind(u[:len(y)], u[len(y):], equal_var=True)[0] \
            )
    try:
        tst = stats[stat]
    except KeyError:
        raise ValueError("Unrecognized test statistic (stat): " + stat)    
    if side == 'greater_than':
        theStat = tst
    elif side == 'less_than':
        theStat = lambda u: -tst(u)
    elif side == 'both':
        theStat = lambda u: math.fabs(tst(u))
    else:
        raise ValueError("Unrecognized side choice: " + side)
    ts = theStat(z)
    hits = np.sum([ (theStat(np.random.permutation(z)) >= ts) for i in range(reps)])
    if CI == 'upper':
        return float(hits)/float(reps), binoUpperCL(reps, hits, cl = CL), ts
    elif CI == 'lower':
        return float(hits)/float(reps), binoLowerCL(reps, hits, cl = CL), ts
    elif CI == 'both':
        return float(hits)/float(reps),  \
                 (binoLowerCL(reps, hits, cl = 1-(1-CL)/2), binoUpperCL(reps, hits, cl = 1-(1-CL)/2)), \
                 ts
    else:
        return float(hits)/float(reps), ts

# <codecell>

def stratifiedPermutationTestMean(group, condition, response, groups, conditions, verbosity = 0):
    tst = 0.0
    '''
    Calculates difference in sample means between treatment conditions, within groups.
    '''
    c0 = condition == conditions[0]
    c1 = condition == conditions[1]
    tst_gg = []
    for g in groups:
        gg = group == g
        x = gg & c0
        y = gg & c1
        if (any(x) & any(y)):
            tst_gg.append(response[x].mean() - response[y].mean())
            tst += tst_gg[-1]
    if verbosity:
        print '\nTest statistic for each stratum:'
        print tst_gg
    return tst


def stratifiedPermutationTestPearsonr(group, condition, response, groups, conditions, verbosity = 0):
    tst = 0.0
    '''
    Calculates abs value of Pearson correlation between treatment condition and Y-Yhat (response), within groups.
    '''
    c0 = condition == conditions[0]
    c1 = condition == conditions[1]
    tst_gg = []
    for g in groups:
        gg = (group == g)
        tst_gg.append(abs(scipy.stats.pearsonr(response[gg],condition[gg])[0]))
        tst += tst_gg[-1]*(gg.mean())
    if verbosity:
        print '\nTest statistic for each stratum:'
        print tst_gg    
    return tst


def permuteWithinGroups(group, condition, groups):
    permuted = condition
    for g in groups:
        gg = group == g
        permuted[gg] = np.random.permutation(condition[gg])      
    return permuted


def stratifiedPermutationTest(group, condition, response, iterations, testStatistic=stratifiedPermutationTestMean, verbosity = 0):
    '''
    Stratified permutation test using the sum of the differences in means between two conditions in
    each group (stratum) as the test statistic.
    The test statistic is
        \sum_{g in groups} [
                            mean(response for cases in group g assigned to first condition) -
                            mean(response for cases in group g assigned to second condition)
                           ].
    There should be at least one group and no more than two conditions.
    Under the null hypothesis, all assignments to the two conditions that preserve the number of
    cases assigned to the conditions are equally likely.
    Groups in which all cases are assigned to the same condition are skipped; they do not contribute 
    to the p-value since all randomizations give the same contribution to the difference in means.
    
    Dependencies: numpy (as np)
    '''   
    groups = np.unique(group)
    conditions = np.unique(condition) 
    if len(conditions) > 2:
        raise ValueError('too many conditions:', conditions)
    elif len(conditions) < 2:
        return 1.0, 1.0, 1.0, np.nan, None
    else:
        tst = testStatistic(group, condition, response, groups, conditions, verbosity)
        dist = np.zeros(iterations)
        for i in range(iterations):
             dist[i] = testStatistic( group, 
                                      permuteWithinGroups(group, condition, groups),
                                      response, groups, conditions, verbosity=0
                                    )
            
    # define the conditions, then map count_nonzero over them
        conds = [dist <= tst, dist >= tst, abs(dist) >= abs(tst)]
        pLeft, pRight, pBoth = np.array(map(np.count_nonzero, conds))/float(iterations)
        return pLeft, pRight, pBoth, tst, dist

# <codecell>


# <codecell>


