# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# Last edited 7 August, 2014 by KO
# <p>
# Implements several types of propensity-score matching, balance diagnostics for group characteristics, average treatment effect estimates, and bootstraping to estimate standard errors of ATE estimates.
# <p>
# Propensity scoring methods:<br>
# <b>One-to-one:</b> Matches individuals in the treatment group to a single individual in the control group.  Variations include whether or not controls are matched with replacement, whether or not a caliper should be used, and the scale of the caliper.
# <br>
# <b>One-to-many:</b> Matches individuals in the treatment group to as many individuals in the control group as possible, subject to some criteria.  Variations include whether or not controls are matched with replacement, caliper methods as in one-to-one matching, and k nearest neighbor matching.

# <codecell>

import math
import numpy as np
import scipy
from scipy.stats import binom, hypergeom
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from ModelMatch import binByQuantiles
import statsmodels.api as sm

# <markdowncell>

# Implement one-to-one matching, caliper without replacement.  Variants of the method are examined in the following paper.  This is something to explore further.
# </br>
# Austin, P. C. (2014), A comparison of 12 algorithms for matching on the propensity score. Statist. Med., 33: 1057â1069. doi: 10.1002/sim.6004

# <codecell>

def computePropensityScore(predictors, groups):
    '''
    Compute propensity scores 
    
    Inputs:
    predictors = DataFrame containing covariates
    groups = Series containing treatment assignment. Indices should match those of predictors.
        Must be 2 groups

    
    dependencies: LogisticRegression from sklearn.linear_model
                  statsmodels as sm
    '''
    
    if len(groups.unique()) != 2:
        raise ValueError('wrong number of groups: expected 2')
        
    ####### Using LogisticRegression from sklearn.linear_model    
    #propensity = LogisticRegression()
    #propensity.fit(predictors, groups)
    #return propensity.predict_proba(predictors)[:,1]
    
    ####### Using sm.GLM
    predictors = sm.add_constant(predictors, prepend=False)
    glm_binom = sm.GLM(groups, predictors, family=sm.families.Binomial())
    res = glm_binom.fit()
    return res.fittedvalues

# <codecell>

def Match(groups, propensity, caliper = 0.05, caliper_method = "propensity", replace = False):
    ''' 
    Implements greedy one-to-one matching on propensity scores.
    
    Inputs:
    groups = Array-like object of treatment assignments.  Must be 2 groups
    propensity = Array-like object containing propensity scores for each observation. Propensity and groups should be in the same order (matching indices)
    caliper = a numeric value, specifies maximum distance (difference in propensity scores or SD of logit propensity) 
    caliper_method = a string: "propensity" (default) if caliper is a maximum difference in propensity scores,
            "logit" if caliper is a maximum SD of logit propensity, or "none" for no caliper
    replace = Logical for whether individuals from the larger group should be allowed to match multiple individuals in the smaller group.
        (default is False)
    
    Output:
    A series containing the individuals in the control group matched to the treatment group.
    Note that with caliper matching, not every treated individual may have a match.
    '''

    # Check inputs
    if any(propensity <=0) or any(propensity >=1):
        raise ValueError('Propensity scores must be between 0 and 1')
    elif not(0<=caliper<1):
        if caliper_method == "propensity" and caliper>1:
            raise ValueError('Caliper for "propensity" method must be between 0 and 1')
        elif caliper<0:
            raise ValueError('Caliper cannot be negative')
    elif len(groups)!= len(propensity):
        raise ValueError('groups and propensity scores must be same dimension')
    elif len(groups.unique()) != 2:
        raise ValueError('wrong number of groups: expected 2')
        
    
    # Transform the propensity scores and caliper when caliper_method is "logit" or "none"
    if caliper_method == "logit":
        propensity = log(propensity/(1-propensity))
        caliper = caliper*np.std(propensity)
    elif caliper_method == "none":
        caliper = 0
    
    # Code groups as 0 and 1
    groups = groups == groups.unique()[0]
    N = len(groups)
    N1 = groups[groups == 1].index; N2 = groups[groups == 0].index
    g1, g2 = propensity[groups == 1], propensity[groups == 0]
    # Check if treatment groups got flipped - the smaller should correspond to N1/g1
    if len(N1) > len(N2):
       N1, N2, g1, g2 = N2, N1, g2, g1
        
        
    # Randomly permute the smaller group to get order for matching
    morder = np.random.permutation(N1)
    matches = {}

    
    for m in morder:
        dist = abs(g1[m] - g2)
        if (dist.min() <= caliper) or not caliper:
            matches[m] = dist.argmin()    # Potential problem: check for ties
            if not replace:
                g2 = g2.drop(matches[m])
    return (matches)



def MatchMany(groups, propensity, method = "caliper", k = 1, caliper = 0.05, caliper_method = "propensity", replace = True):
    ''' 
    Implements greedy one-to-many matching on propensity scores.
    
    Inputs:
    groups = Array-like object of treatment assignments.  Must be 2 groups
    propensity = Array-like object containing propensity scores for each observation. Propensity and groups should be in the same order (matching indices)
    method = a string: "caliper" (default) to select all matches within a given range, "knn" for k nearest neighbors,
    k = an integer (default is 1). If method is "knn", this specifies the k in k nearest neighbors
    caliper = a numeric value, specifies maximum distance (difference in propensity scores or SD of logit propensity) 
    caliper_method = a string: "propensity" (default) if caliper is a maximum difference in propensity scores,
            "logit" if caliper is a maximum SD of logit propensity, or "none" for no caliper
    replace = Logical for whether individuals from the larger group should be allowed to match multiple individuals in the smaller group.
        (default is True)
    
    Output:
    A series containing the individuals in the control group matched to the treatment group.
    Note that with caliper matching, not every treated individual may have a match within calipers.
        In that case we match it to its single nearest neighbor.  The alternative is to throw out individuals with no matches, but then we'd no longer be estimating the ATT.
    '''

    # Check inputs
    if any(propensity <=0) or any(propensity >=1):
        raise ValueError('Propensity scores must be between 0 and 1')
    elif not(0<=caliper<1):
        if caliper_method == "propensity" and caliper>1:
            raise ValueError('Caliper for "propensity" method must be between 0 and 1')
        elif caliper<0:
            raise ValueError('Caliper cannot be negative')
    elif len(groups)!= len(propensity):
        raise ValueError('groups and propensity scores must be same dimension')
    elif len(groups.unique()) != 2:
        raise ValueError('wrong number of groups: expected 2')
        
    
    # Transform the propensity scores and caliper when caliper_method is "logit" or "none"
    if method == "caliper":
        if caliper_method == "logit":
            propensity = log(propensity/(1-propensity))
            caliper = caliper*np.std(propensity)
        elif caliper_method == "none":
            caliper = 0
    
    # Code groups as 0 and 1
    groups = groups == groups.unique()[0]
    N = len(groups)
    N1 = groups[groups == 1].index; N2 = groups[groups == 0].index
    g1, g2 = propensity[groups == 1], propensity[groups == 0]
    # Check if treatment groups got flipped - the smaller should correspond to N1/g1
    if len(N1) > len(N2):
       N1, N2, g1, g2 = N2, N1, g2, g1
        
        
    # Randomly permute the smaller group to get order for matching
    morder = np.random.permutation(N1)
    matches = {}
    
    for m in morder:
        dist = abs(g1[m] - g2)
        dist.sort()
        if method == "knn":
            caliper = dist.iloc[k-1]
        # PROBLEM: when there are ties in the knn. 
        # Need to randomly select among the observations tied for the farthest eacceptable distance
        keep = np.array(dist[dist<=caliper].index)
        if len(keep):
            matches[m] = keep
        else:
            matches[m] = [dist.argmin()]
        if not replace:
            g2 = g2.drop(matches[m])
    return (matches)


    
def whichMatched(matches, data, many = False, unique = False):
    ''' 
    Simple function to convert output of Matches to DataFrame of all matched observations
    Inputs:
    matches = output of Match
    data = DataFrame of covariates
    many = Boolean indicating if matching method is one-to-one or one-to-many
    unique = Boolean indicating if duplicated individuals (ie controls matched to more than one case) should be removed
    '''

    tr = matches.keys()
    if many:
        ctrl = [m for matchset in matches.values() for m in matchset]
    else:
        ctrl = matches.values()
    # need to remove duplicate rows, which may occur in matching with replacement
    temp = pd.concat([data.ix[tr], data.ix[ctrl]])
    if unique == True:
        return temp.groupby(temp.index).first()
    else:
        return temp

# <codecell>

#def averageTreatmentEffect(groups, response, matches):
#    '''
#    This works for one-to-one matching
#    
#    Inputs:
#    groups = Series containing treatment assignment. Must be 2 groups
#    response = Series containing response measurements. Indices should match those of groups.
#    matches = matched pairs returned from Match
#    '''
#    if len(groups.unique()) != 2:
#        raise ValueError('wrong number of groups: expected 2')
#    
#    groups = (groups == groups.unique()[0])            
#    response1, response0 = response[groups==1], response[groups==0]
#    response1_matched = response1[matches.notnull()]
#    response0_matched = response0[matches]
#    response0_matched = response0_matched[response0_matched.notnull()]
#    return response1.mean() - response0_matched.mean()

def averageTreatmentEffect(groups, response, matches):
    '''
    This works for one-to-one matching.  The data passed in should already have unmatched individuals (and duplicates?) removed.
    Weights argument will be added later for one-many matching
    
    Inputs:
    groups = Series containing treatment assignment. Must be 2 groups
    response = Series containing response measurements. Indices should match those of groups.
    matches = output of Match or MatchMany
    '''
    if len(groups.unique()) != 2:
        raise ValueError('wrong number of groups: expected 2')
    
    groups = (groups == groups.unique()[0])
    response = response.groupby(response.index).first()
    response1 = []; response0 = []
    for k in matches.keys():
        response1.append(response[k])
        response0.append( (response[matches[k]]).mean() )
    #response1, response0 = response[groups==1], response[groups==0]
    return np.mean( np.array(response1)-np.array(response0) )


def bootstrapATE(groups, response, propensity, B = 500, caliper = 0.05, caliper_method = "propensity", replace = False):
    '''
    Computes bootstrap standard error of the average treatment effect
    Sample observations with replacement, within each treatment group. Then match them and compute ATE
    Repeat B times and take standard deviation
    
    Inputs:
    groups = Series containing treatment assignment. Must be 2 groups
    response = Series containing response measurements
    propensity = Series containing propensity scores
    B = number of bootstrap replicates. Default is 500
    caliper, replace = arguments to pass to Match    
    '''
    if len(groups.unique()) != 2:
        raise ValueError('wrong number of groups: expected 2')

    data = pd.DataFrame({'groups':groups, 'response':response, 'propensity':propensity})
    boot_ate = np.empty(B)
    for i in range(B):
        bootdata = data.copy()
        for g in groups.unique():
            sample = np.random.choice(data.index[data.groups==g], sum(groups == g), replace = True)
            newdata =(data[data.groups==g]).ix[sample]
            newdata.index = bootdata.index[bootdata.groups == g]
            bootdata[bootdata.groups == g] = newdata
        pairs = Match(bootdata.groups, bootdata.propensity, caliper = caliper, caliper_method = caliper_method, replace = replace)
        boot_ate[i] = averageTreatmentEffect(bootdata.groups, bootdata.response, matches = pairs)
    return boot_ate.std()

def bootstrapManyATE(groups, response, propensity, B = 500, method = "caliper", k = 1, caliper = 0.05, caliper_method = "propensity", replace = True):
    '''
    Computes bootstrap standard error of the average treatment effect
    Sample observations with replacement, within each treatment group. Then match them and compute ATE
    Repeat B times and take standard deviation
    
    Inputs:
    groups = Series containing treatment assignment. Must be 2 groups
    response = Series containing response measurements
    propensity = Series containing propensity scores
    B = number of bootstrap replicates. Default is 500
    caliper, replace = arguments to pass to Match    
    '''
    if len(groups.unique()) != 2:
        raise ValueError('wrong number of groups: expected 2')

    data = pd.DataFrame({'groups':groups, 'response':response, 'propensity':propensity})
    boot_ate = np.empty(B)
    for i in range(B):
        bootdata = data.copy()
        for g in groups.unique():
            sample = np.random.choice(data.index[data.groups==g], sum(groups == g), replace = True)
            newdata =(data[data.groups==g]).ix[sample]
            newdata.index = bootdata.index[bootdata.groups == g]
            bootdata[bootdata.groups == g] = newdata
        pairs = MatchMany(bootdata.groups, bootdata.propensity, method = method, k = k, caliper = caliper, caliper_method = caliper_method, replace = replace)
        boot_ate[i] = averageTreatmentEffect(bootdata.groups, bootdata.response, matches = pairs)
    return boot_ate.std()


def regressAverageTreatmentEffect(groups, response, covariates, matches=None, verbosity = 0):
    '''
    This works for one-to-one matching.   The data passed in should already have unmatched individuals removed.
    Weights argument will be added later for one-many matching
    
    Inputs:
    groups = Series containing treatment assignment. Must be 2 groups
    response = Series containing response measurements. Indices should match those of groups.
    covariates = DataFrame containing the covariates to include in the linear regression
    matches = optional: if using one-many matching, should be the output of MatchMany.
            Use None for one-one matching.
    
    Dependencies: statsmodels.api as sm, pandas as pd
    '''
    if len(groups.unique()) != 2:
        raise ValueError('wrong number of groups: expected 2')
    
    weights = pd.Series(data = np.ones(len(groups)), index = groups.index)
    if matches:
        ctrl = [m for matchset in matches.values() for m in matchset]    
        matchcounts = pd.Series(ctrl).value_counts()
        for i in matchcounts.index:
            weights[i] = matchcounts[i]
        if verbosity:
            print weights.value_counts(), weights.shape
    X = pd.concat([groups, covariates], axis=1)
    X = sm.add_constant(X, prepend=False)
    linmodel = sm.WLS(response, X, weights = weights).fit()
    return linmodel.params[0], linmodel.bse[0]

# <codecell>

def Balance(groups, covariates):
    '''
    Computes absolute difference of means and standard error for covariates by group
    '''
    means = covariates.groupby(groups).mean()
    dist = abs(means.diff()).ix[1]
    std = covariates.groupby(groups).std()
    n = groups.value_counts()
    se = std.apply(lambda(s): np.sqrt(s[0]**2/n[0] + s[1]**2/n[1]))
    return dist, se

# <codecell>

def Stratify(groups, response, propensity, nbins = 5, verbosity = 0):
    ''' 
    Implements propensity score stratification on quantiles and computes average treatment effects
    
    Inputs:
    groups = Array-like object of treatment assignments.  Must be 2 groups
    response = Array-like object containing response measurements
    propensity = Array-like object containing propensity scores for each observation. Propensity and groups should be in the same order (matching indices)
    nbins = number of stratification groups of approximately equal size. 
            Default is to stratify on the propensity score quintiles (5)
    verbosity = flag for printing and debugging. 
        0 for no printed output, 1 for some verbosity, 2 for maximum verbosity

    dependencies: ModelMatch.py
    
    Output:

    '''    
    
    if any(propensity <=0) or any(propensity >=1):
        raise ValueError('Propensity scores must be between 0 and 1')
    elif len(groups)!= len(propensity):
        raise ValueError('groups and propensity scores must be same dimension')
    elif len(groups.unique()) != 2:
        raise ValueError('wrong number of groups: expected 2')
    
    groups = (groups == groups.unique()[0])
    bins = binByQuantiles(propensity, nbins = nbins, verbosity = verbosity)
    ate = np.empty(nbins)
    for b in arange(nbins):
        stratum = (bins == b)
        response0 = response[(stratum==1) & (groups==0)]
        response1 = response[(stratum==1) & (groups==1)]
        ate[b] = response1.mean() - response0.mean()
        # For diagnostics, print # of group0 and group1 in each stratum
    keep = where(isnan(ate) == False)[0]
    pooled_ate = np.average(ate[keep], weights = bins.value_counts().order(ascending=False)[keep])
    return ate, pooled_ate
        
        

# <codecell>


