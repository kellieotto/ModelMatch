



### Alt hypothesis: data from when the law was in effect comes from a different (mean) distribution than 
###   			from when the law was not in effect.  Do a matched pair permutation test


### Carry out matching on model predictions. Here we'll just do 1-1 matching

Matches <- function(treatment, prediction){
  ### function to match treatment to control based on their predicted outcomes
  ### Inputs: treatment vector, predicted outcomes vector
  ### Eventually need to input optional arguments for matching options to pass to Match
  ### Output: a list where each entry is a group of matched individuals with their treatment assignments
  pairs <- Match(Tr = treatment, X = prediction, replace = TRUE, M=1, ties = FALSE)
  pairs <- pairs$MatchLoopC[,1:2]
  strata <- lapply(unique(pairs[,1]), function(i){
                  matches <- pairs[pairs[,1]==i,2]
                  stratum <- c(i,matches)
                  tr      <- treatment[stratum]
                  cbind(stratum,tr)
                  })
  return(strata)
}



permute_within_groups <- function(groups){
  ### Permute treatment assignment within groups of matched individuals
  ### Input: groups = list from Matches (or from stratification function)
  permuted <- groups
  for(g in 1:length(groups)){
    permuted[[g]][,"tr"] <- sample(permuted[[g]][,"tr"])
  }
  return(permuted)
}

within_group_mean <- function(groups, prediction, response){
  ### Compute difference in mean residual between treated and control within each group
  ### Input: groups     = list from Matches (or from stratification function)
  ###        prediction = vector of predictions from the model
  ###        response   = vector of observed responses for each individual
  resid <- prediction-response
  sapply(groups, function(g){
    treated <- g[g[,2]==1,1]
    ctrl    <- g[g[,2]==0,1]
    return(mean(resid[treated]) - mean(resid[ctrl]))
  })
}

permu_test_mean <- function(groups, prediction, response, iters=1000){
  ### Carry out stratified permutation test for difference in mean residuals
  ### Input: groups     = list from Matches (or from stratification function)
  ###        prediction = vector of predictions from the model
  ###        response   = vector of observed responses for each individual
  ###        iters      = number of permutations (default 1000)
  truth <- sum(within_group_mean(groups, prediction, response))
  perm_dist <- replicate(iters, {perm_groups <- permute_within_groups(groups)
                                 sum(within_group_mean(perm_groups, prediction, response))})
  pval <- c("p_upper" = sum(perm_dist >= truth)/iters,
            "p_lower" = sum(perm_dist <= truth)/iters,	# 
            "twosided" = sum(abs(perm_dist) >= abs(truth))/iters)	# Alt hypoth: residual mean after law effected different from residual mean for pre-law
  return(list("diff_means"=truth, "perm_distribution"=perm_dist, "pvalue"=pval))
}
