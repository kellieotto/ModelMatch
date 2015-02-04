###################################################### Modeling and matching ######################################################

#' Match on model predictions.
#'
#' At this point, Matches only supports one-to-one matching of treatment and control groups.
#' @param treatment Vector of treatment values
#' @param prediction Vector of predicted outcomes
#' @return a list. Each entry is a group of matched individuals with their treatments.
Matches <- function(treatment, prediction){
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



###################################################### Stratified permutation test for difference in mean residuals ######################################################

#' Permute within groups
#'
#' Permute treatment assignment within groups of matched individuals
#' @param groups List output from Matches, containing matched pairs or groups
#' @return a list of the same structure as the input, with treatment assignments permuted.
permute_within_groups <- function(groups){
  permuted <- groups
  for(g in 1:length(groups)){
    permuted[[g]][,"tr"] <- sample(permuted[[g]][,"tr"])
  }
  return(permuted)
}



#' Difference in means within groups
#'
#' Compute difference in mean residual between treated and control within each group
#' @inheritParams permu_test_mean
#' @return a vector of differences
within_group_mean <- function(groups, prediction, response, shift = 0){
  ### Compute difference in mean residual between treated and control within each group
  ### Input: groups     = list from Matches (or from stratification function)
  ###        prediction = vector of predictions from the model
  ###        response   = vector of observed responses for each individual
  ###        shift      = null shift (constant scalar) for treatment group (default 0)
  resid <- prediction-response
  sapply(groups, function(g){
    treated <- g[g[,2]==1,1]
    ctrl    <- g[g[,2]==0,1]
    return(-1*(mean(resid[treated] - shift) - mean(resid[ctrl])))
  })
}




#' Test for the difference in means within groups
#'
#' Carry out stratified permutation test for difference in mean residuals between the treatment and control groups.
#' @param groups List output from Matches, containing matched pairs or groups
#' @param prediction Vector of predicted outcomes
#' @param response Vector of responses
#' @param treatment Vector of treatments
#' @param iters The number of Monte Carlo iterations (default 1000)
#' @param shift Null shift (a scalar constant) for the treatment group (default 0)
#' @return a list containing the results: attributes diff_means (the estimated difference), perm_distribution (simulated permutation distribution), and pvalue (p-value for the test)
permu_test_mean <- function(groups, prediction, treatment, response, iters=1000, shift = 0){
  truth <- sum(within_group_mean(groups, prediction, response))
  response_shift <- response - shift*treatment
  perm_dist <- replicate(iters, {perm_groups <- permute_within_groups(groups)
                                 sum(within_group_mean(perm_groups, prediction, response_shift, shift = shift))})
  pval <- c("p_upper" = sum(perm_dist >= truth)/iters,
            "p_lower" = sum(perm_dist <= truth)/iters,	#
            "twosided" = sum(abs(perm_dist) >= abs(truth))/iters)	# Alt hypoth: residual mean after law effected different from residual mean for pre-law
  return(list("diff_means"=truth/length(groups), "perm_distribution"=perm_dist, "pvalue"=pval))
}




#' Confidence interval for the difference in means within groups
#'
#' Invert the stratified permutation test to get a \eqn{1-\alpha} confidence interval for the difference in mean residuals
#' @inheritParams permu_test_mean
#' @param side Type of interval, either "both", "upper", or "lower". Default is "both".
#' @param alpha Significance level
#' @param precision Rate at which to iteratively increment the confidence bounds, relative to the estimated difference in means. The value should be between 0 and 1. Smaller precision results in a slower run time but more precise confidence bounds (default is 0.02)
#' @param verbose Verbosity switch - print the p-value and confidence interval endpoint at each step? (Default FALSE)
#' @return a vector of differences
permu_CI_mean <- function(groups, prediction, response, treatment, side = "both", alpha=0.05, iters=1000, shift = 0, precision = 0.02, verbose = FALSE){
  if(side == "both"){
    alpha <- alpha/2
    d1 <- permu_CI_mean(groups = groups, prediction = prediction, treatment = treatment, response = response, side = "lower", alpha = alpha, iters = iters, precision = precision, verbose = verbose)
    d2 <- permu_CI_mean(groups = groups, prediction = prediction, treatment = treatment, response = response, side = "upper", alpha = alpha, iters = iters, precision = precision, verbose = verbose)
    return(c(d1,d2))
  }
  pvalue_side <- which(c("lower", "upper", "both") == side); print(pvalue_side)

  # initialize
  tr <- sapply(strata, function(x) x[x[,"tr"] == 1,"stratum"])
#  response_alt <- response
  res <- permu_test_mean(groups = groups, prediction = prediction, treatment = treatment, response = response, iters = iters)
  d_true <- res$diff_means; d <- res$diff_means; incr <- abs(d)*precision
  shift_pval <- rep(1,3)

  # Conduct permutation test for H0: shift = d until we reject H0
  while(shift_pval[pvalue_side] > alpha){
    d <- ifelse(side == "upper", d+incr, d-incr)
#    response_alt[tr] <- response[tr] + d
    res <- permu_test_mean(groups, prediction, treatment, response, shift = d, iters = iters)
    shift_pval <- res$pvalue
    if(verbose){print(d); print(shift_pval)}
  }
  return(d)
}




############################################# Pearson correlation - no stratification #################################################################################


#' Test the correlation between residuals and treatment
#'
#' Carry out a permutation test for the Pearson correlation between residuals and treatment.
#' @param prediction Vector of predicted outcomes
#' @param response Vector of responses
#' @param treatment Vector of treatments
#' @param iters The number of Monte Carlo iterations (default 1000)
#' @return a list containing the results: attributes estimate (the estimated difference), distr (simulated permutation distribution), and pvalue (p-value for the test)
permu_pearson <- function(prediction, response, treatment, iters = 1000){
  resid <- prediction-response
  pearson_r <- cor(resid, treatment)
  distr <- replicate(iters, {
    tr <- sample(treatment)
    cor(resid, tr)
  })
  pval <- c("p_upper" = sum(distr >= pearson_r)/iters,
            "p_lower" = sum(distr <= pearson_r)/iters,
            "twosided" = sum(abs(distr) >= abs(pearson_r))/iters)
  return(list("estimate" = pearson_r, "distr" = distr, "pvalue" = pval))
}



#' Confidence interval for the correlation between residuals and treatment
#'
#' Bootstrap to get a \eqn{1-\alpha} confidence interval for the Pearson correlation between residuals and treatment
#' @inheritParams permu_pearson_shift
#' @inheritParams permu_CI_pearson
#' @return a confidence interval (vector)
bootstrap_CI_pearson <- function(prediction, outcome, treatment, iters, alpha = 0.05, side = "both"){
  resid <- prediction-outcome
  distr <- replicate(iters, {
    resid_boot <- base::sample(resid, length(resid), replace=TRUE)
    cor(resid_boot, treatment)
  })
  if(side == "both"){
    alpha <- alpha/2
    quant <- c(floor(length(distr)*alpha), length(distr)-floor(length(distr)*alpha))
  }else{
    if(side == "upper"){
      quant <- length(distr) - floor(length(distr)*alpha)
    }else{quant <- floor(length(distr)*alpha)}
  }
  return(sort(distr)[quant])
}




#' Test the correlation between residuals and treatment
#'
#' Carry out a permutation test for the Pearson correlation between residuals and treatment.
#' @param prediction Vector of predicted outcomes
#' @param response Vector of responses
#' @param treatment Vector of treatments
#' @param iters The number of Monte Carlo iterations (default 1000)
#' @param shift Null shift away from 0 (a scalar constant) in correlation (default 0)
#' @return a list containing the results: attributes estimate (the estimated difference), distr (simulated permutation distribution), and pvalue (p-value for the test)
permu_pearson_shift <- function(prediction, response, treatment, iters = 1000, shift = 0){
  resid <- prediction-response
  pearson_r <- cor(resid, treatment)
  distr <- replicate(iters, {
    tr <- sample(treatment)
    cor(resid, tr) + shift
  })
  pval <- c("p_upper" = sum(distr >= pearson_r)/iters,
            "p_lower" = sum(distr <= pearson_r)/iters,
            "twosided" = sum(abs(distr) >= abs(pearson_r))/iters)
  return(list("estimate" = pearson_r, "distr" = distr, "pvalue"=pval))
}



#' Confidence interval for the correlation between residuals and treatment
#'
#' Invert the permutation test to get a \eqn{1-\alpha} confidence interval for the Pearson correlation between residuals and treatment
#' @inheritParams permu_pearson_shift
#' @param side Type of interval, either "both", "upper", or "lower". Default is "both".
#' @param alpha Significance level
#' @param verbose Verbosity switch - print the p-value and confidence interval endpoint at each step? (Default FALSE)
#' @return a confidence interval (vector)
permu_CI_pearson <- function(prediction, response, treatment, iters, alpha = 0.05, side = "both", verbosity = FALSE){
  resid <- prediction-response
  if(side == "both"){
    d1 <- permu_CI_pearson(prediction, response, treatment, iters, alpha = alpha/2, side = "lower", verbosity = verbosity)
    d2 <- permu_CI_pearson(prediction, response, treatment, iters, alpha = alpha/2, side = "upper", verbosity = verbosity)
    return(c(d1,d2))
  }
  which_p <- which(c("lower", "upper", "both") %in% side); if(verbosity){print(side)}
  pval <- rep(1,3)
  d <- ifelse(side == "upper", 0.005, -0.005); shift <- d
  while(pval[which_p] >= alpha){
    res <- permu_pearson_shift(prediction, response, treatment, iters, shift = shift)
    pval <- res$pvalue; shift <- shift+d
    if(verbosity == TRUE){cat(shift-d, "\n", pval, "\n")}
  }
  return(res$estimate +  shift)
}
