



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

shift_treated <- function(groups, shift=0){
  ### Remove the null shift from the treatment group
  ### Input: groups     = list from Matches (or from stratification function)
  ###        shift      = shift (constant scalar) in treatment group response under the null (default 0)
  lapply(groups, function(x) {x[x[,"tr"]==1,1] = x[x[,"tr"]==1,1] - shift; return(x)})
}

permu_test_mean <- function(groups, prediction, response, iters=1000, shift = 0){
  ### Carry out stratified permutation test for difference in mean residuals
  ### Input: groups     = list from Matches (or from stratification function)
  ###        prediction = vector of predictions from the model
  ###        response   = vector of observed responses for each individual
  ###        iters      = number of permutations (default 1000)
  ###        shift      = shift (constant scalar) in treatment group response under the null (default 0)
  truth <- sum(within_group_mean(groups, prediction, response))
  groups_unshifted <- shift_treated(groups, shift = shift)
  perm_dist <- replicate(iters, {perm_groups <- permute_within_groups(groups)
                                 groups_reshifted <- shift_treated(groups, shift = -shift)
                                 sum(within_group_mean(perm_groups, prediction, response))})
  pval <- c("p_upper" = sum(perm_dist >= truth)/iters,
            "p_lower" = sum(perm_dist <= truth)/iters,	#
            "twosided" = sum(abs(perm_dist) >= abs(truth))/iters)	# Alt hypoth: residual mean after law effected different from residual mean for pre-law
  return(list("diff_means"=truth/length(groups), "perm_distribution"=perm_dist, "pvalue"=pval))
}


permu_CI_mean <- function(groups, prediction, response, side = "both", alpha=0.05, iters=1000, precision = 50, verbose = FALSE){
  ### Invert the stratified permutation test to get a 1-alpha confidence interval for the difference in mean residuals
  ### Input: groups     = list from Matches (or from stratification function)
  ###        prediction = vector of predictions from the model
  ###        response   = vector of observed responses for each individual
  ###        side       = Type of interval, either "both", "upper", or "lower". Default is "both"
  ###        alpha      = Significance level
  ###        iters      = number of permutations (default 1000)
  ###        precision  = rate at which to increase the confidence bounds under consideration.
  ###                     higher = slower run time but more precise. Default is 50 (relative to estimated diff in means)
  ###       verbose     = indicator whether to print p-values associated with each confidence interval endpoint considered.
  if(side == "both"){
    alpha <- alpha/2
    d1 <- permu_CI_mean(groups = groups, prediction = prediction, response = response, side = "lower", alpha = alpha, iters = iters, precision = precision, verbose = verbose)
    d2 <- permu_CI_mean(groups = groups, prediction = prediction, response = response, side = "upper", alpha = alpha, iters = iters, precision = precision, verbose = verbose)
    return(c(d1,d2))
  }
  pvalue_side <- which(c("lower", "upper", "both") == side)

  # initialize
  tr <- sapply(strata, function(x) x[x[,"tr"] == 1,"stratum"])
  response_alt <- response
  res <- permu_test_mean(groups, prediction, response, iters = iters)
  d_true <- res$diff_means; d <- res$diff_means; incr <- d/precision
  shift_pval <- rep(1,3)

  # Conduct permutation test for H0: shift = d until we reject H0
  while(shift_pval[pvalue_side] > alpha){
    d <- ifelse(side == "upper", d+incr*sign(incr), d+incr)
    response_alt[tr] <- response[tr] + d
    res <- permu_test_mean(groups, prediction, response_alt, iters = iters)
    shift_pval <- res$pvalue
    if(verbose){print(d); print(shift_pval)}
  }
  return(d)
}
