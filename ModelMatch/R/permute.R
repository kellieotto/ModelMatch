#library(Rcpp)
#sourceCpp("R/permute.cpp")
###################################################### Modeling and matching ######################################################

#' Match on model predictions.
#'
#' At this point, Matches only supports one-to-one matching of treatment and control groups.
#' @param treatment Vector of treatment values
#' @param prediction Vector of predicted outcomes
#' @return a list. Each entry is a group of matched individuals with their treatments.
Matches <- function(treatment, prediction){
  treated <- which(treatment == 1)
  treated <- sample(treated)  # Put them in a random order to reduce bias
  controls <- which(treatment == 0)
  matched_controls <- rep(NA, length(treated))
  ind <- 0
  for(i in treated){
    ind <- ind+1
    distance <- abs(prediction[i] - prediction[controls])
    matched_controls[ind] <- controls[which.min(distance)]
  }
  pairs <- cbind(treated, matched_controls)
  strata <- lapply(unique(pairs[,1]), function(i){
    matches <- pairs[pairs[,1]==i,2]
    names(matches) <- NULL
    stratum <- c(i,matches)
    data.frame("index" = stratum, "score" = prediction[stratum], "treatment" = treatment[stratum])
  })
  return(strata)
}

#' Stratify on model predictions.
#'
#' Note, we need sufficient distribution of the treatment values within each bin. User should check this.
#' @param treatment Vector of treatment values
#' @param prediction Vector of predicted outcomes
#' @param quantiles Vector of k-1 values that divide the k strata
#' @param strata Number of strata desired. Default is 5.
#' @return a list. Each entry is a group of matched individuals with their treatments.
Strata <- function(treatment, prediction, breaks = NULL, strata = 5){
  if(is.null(breaks)){
    if(is.null(strata)){
      stop("You must specify breaks or strata")
    }else{
      breaks <- quantile(prediction, probs = seq(0, 1, by = 1/strata))
    }
  }else{
    breaks <- sort(breaks)
    breaks <- c(min(c(prediction, breaks[1]-1)), breaks, max(c(prediction, breaks[length(breaks)]+1)))
  }
  group_labels <- cut(prediction, breaks, include.lowest = TRUE)
  dat <- data.frame("index" = 1:length(treatment), "score" = prediction, "treatment" = treatment)
  strata <- lapply(unique(group_labels), function(g) dat[group_labels == g,])
  return(strata)
}

###################################################### Stratified permutation test for difference in mean residuals ######################################################

#' Permute within groups
#'
#' Permute treatment assignment within groups of matched individuals
#' @param strata List output from Matches, containing matched pairs or groups
#' @return a list of the same structure as the input, with treatment assignments permuted.
permute_within_groups <- function(strata){
  permuted <- strata
  for(g in 1:length(strata)){
    permuted[[g]][,"treatment"] <- sample(strata[[g]][,"treatment"])
  }
  return(permuted)
}



#' Difference in means within groups
#'
#' Compute difference in mean residual between treated and control within each group
#' @inheritParams permu_test_mean
#' @return a vector of differences
within_group_mean <- function(strata, prediction, response, shift = 0){
  resid <- prediction-response
  sapply(strata, function(g){
    tt <- (g$treatment == 1)
    treated <- g$index[tt == 1]
    ctrl <- g$index[tt == 0]
    return(-1*(mean(resid[treated] - shift) - mean(resid[ctrl])))
  })
}




#' Test for the difference in means within groups
#'
#' Carry out stratified permutation test for difference in mean residuals between the treatment and control groups.
#' @param strata List output from Matches, containing matched pairs or groups
#' @param prediction Vector of predicted outcomes
#' @param response Vector of responses
#' @param treatment Vector of treatments
#' @param iters The number of Monte Carlo iterations (default 1000)
#' @param shift Null shift (a scalar constant) for the treatment group (default 0)
#' @return a list containing the results: attributes diff_means (the estimated difference), perm_distribution (simulated permutation distribution), and pvalue (p-value for the test)
permu_test_mean_slow <- function(strata, prediction, treatment, response, iters=1000, shift = 0){
  GG <- length(strata)
  if(length(prediction) == 1){prediction <- rep(prediction, length(treatment))}
  if(length(response) == 1){response <- rep(response, length(treatment))}

  truth <- sum(abs(within_group_mean(strata, prediction, response)))/GG
  response_shift <- response - shift*treatment
  perm_dist <- replicate(iters, {perm_groups <- permute_within_groups(strata)
                                 sum(abs(within_group_mean(perm_groups, prediction, response_shift, shift = shift)))/GG})
  pval <- c("p_upper" = sum(perm_dist >= truth)/iters,
            "p_lower" = sum(perm_dist <= truth)/iters)
  pval <- c(pval, "twosided" = min(1, 2*min(pval)))
  return(list("diff_means"=truth, "perm_distribution"=perm_dist, "pvalue"=pval))
}


#' Test for the difference in means within groups
#'
#' Carry out stratified permutation test for difference in mean residuals between the treatment and control groups.
#' @param strata List output from Matches, containing matched pairs or groups
#' @param prediction Vector of predicted outcomes
#' @param response Vector of responses
#' @param treatment Vector of treatments
#' @param iters The number of Monte Carlo iterations (default 1000)
#' @param shift Null shift (a scalar constant) for the treatment group (default 0)
#' @return a list containing the results: attributes diff_means (the estimated difference), perm_distribution (simulated permutation distribution), and pvalue (p-value for the test)
permu_test_mean <- function(strata, prediction, treatment, response, iters=1000, shift = 0){
  GG <- length(strata)
  if(length(prediction) == 1){prediction <- rep(prediction, length(treatment))}
  if(length(response) == 1){response <- rep(response, length(treatment))}

  truth <- sum(abs(within_group_mean_cpp(strata, prediction, response)))/GG
  response_shift <- response - shift*treatment
  perm_groups <- strata
  perm_dist <- replicate(iters, {
    perm_groups <- permute_within_groups_cpp(perm_groups)
    sum(abs(within_group_mean_cpp(perm_groups, prediction, response_shift, shift = shift)))
    })/GG
  pval <- c("p_upper" = sum(perm_dist >= truth)/iters,
            "p_lower" = sum(perm_dist <= truth)/iters)
  pval <- c(pval, "twosided" = min(1, 2*min(pval)))
  return(list("diff_means"=truth, "perm_distribution"=perm_dist, "pvalue"=pval))
}



#' Test for the difference in means within groups, refitting the model at each permutation
#'
#' Carry out stratified permutation test for difference in mean residuals between the treatment and control groups.
#' @inheritParams permu_test_mean
#' @param fit_function Function that takes arguments covariates and responses on which to fit the prediction model and returns that model object
#' example: fit_function = function(X, Y) lm(Y~., data.frame(X, Y))
#' @return a list containing the results: attributes diff_means (the estimated difference), perm_distribution (simulated permutation distribution), and pvalue (p-value for the test)
permu_test_mean_refit <- function(strata, fit_function, X, response, iters=1000, shift = 0){
  GG <- length(strata)
  if(length(response) == 1){response <- rep(response, length(treatment))}

  tmp <- do.call(rbind, strata)
  treatment <- tmp[order(tmp[,"index"]),"treatment"]
  prediction <- predict(fit_function(X[treatment == 0, ], response[treatment == 0]), X)
  truth <- sum(abs(within_group_mean_cpp(strata, prediction, response)))/GG
  response_shift <- response - shift*treatment
  perm_groups <- strata
  perm_dist <- replicate(iters, {
    perm_groups <- permute_within_groups_cpp(perm_groups)
    perm_groups_df <- do.call(rbind, perm_groups)
    treatment <- perm_groups_df[order(perm_groups_df[,"index"]),"treatment"]
    prediction <- predict(fit_function(X[treatment == 0, ], response[treatment == 0]), X)
    sum(abs(within_group_mean_cpp(perm_groups, prediction, response_shift, shift = shift)))
  })/GG
  pval <- c("p_upper" = sum(perm_dist >= truth)/iters,
            "p_lower" = sum(perm_dist <= truth)/iters)
  pval <- c(pval, "twosided" = min(1, 2*min(pval)))
  return(list("diff_means"=truth, "perm_distribution"=perm_dist, "pvalue"=pval))
}


#' Confidence interval for the difference in means within groups
#'
#' Invert the stratified permutation test to get a \eqn{1-\alpha} confidence interval for the difference in mean residuals
#' @inheritParams permu_test_mean
#' @param side Type of interval, either "both", "upper", or "lower". Default is "both".
#' @param alpha Significance level
#' @param precision Rate at which to iteratively increment the confidence bounds, relative to the estimated difference in means. The value should be between 0 and 1. Smaller precision results in a slower run time but more precise confidence bounds (default is 0.02)
#' @return a vector of differences
permu_CI_mean <- function(groups, prediction, response, treatment, side = "both", alpha=0.05, iters=1000){
  if(side == "both"){
    alpha <- alpha/2
    d1 <- permu_CI_mean(groups = groups, prediction = prediction, treatment = treatment, response = response, side = "lower", alpha = alpha, iters = iters)
    d2 <- permu_CI_mean(groups = groups, prediction = prediction, treatment = treatment, response = response, side = "upper", alpha = alpha, iters = iters)
    return(c(d1,d2))
  }
  pvalue_side <- which(c("lower", "upper", "both") == side)

  # initialize
  tr <- sapply(groups, function(x) x[x[,2] == 1,"stratum"])
  res <- permu_test_mean(groups = groups, prediction = prediction, treatment = treatment, response = response, iters = 1)
  d_prev <- res$diff_means; d_next <- res$diff_means; incr <- (max(prediction-response) - min(prediction-response))/3
  shift_pval <- rep(1,3)
  shift_permtest <- function(ss, pval = FALSE){
      # if pval = TRUE, return ONLY the relevant p-value
      res <- permu_test_mean(groups, prediction, treatment, response, shift = ss, iters = iters)
      if(pval == TRUE){return(res$pvalue[pvalue_side])}
      return(res)
    }

  # Conduct permutation test for H0: shift = d until we reject H0
  while(shift_pval[pvalue_side] > alpha){
    d_prev <- d_next
    d_next <- ifelse(side == "upper", d_prev+incr, d_prev-incr)
    res <- shift_permtest(d_next)
    shift_pval <- res$pvalue
  }

  # Bisection
  if(side == "upper"){
    l_int <- d_prev; u_int <- d_next
  }else{
    l_int <- d_next; u_int <- d_prev
  }
  d <- uniroot(function(x) {shift_permtest(x, pval=TRUE)-alpha}, lower = l_int, upper = u_int, tol = alpha/iters)
  return(d$root)
}




############################################# Pearson correlation - no stratification #################################################################################


#' Test the correlation between residuals and treatment
#'
#' Carry out a permutation test for the Pearson correlation between residuals and treatment.
#' @param prediction Vector of predicted outcomes
#' @param response Vector of responses
#' @param treatment Vector of treatments
#' @param iters The number of Monte Carlo iterations (default 1000)
#' #' @param shift Null shift away from 0 (a scalar constant) in correlation (default 0)
#' @return a list containing the results: attributes estimate (the estimated difference), distr (simulated permutation distribution), and pvalue (p-value for the test)
permu_pearson <- function(prediction, response, treatment, iters = 1000, shift = 0){
  resid <- prediction-response
  pearson_r <- cor(resid, treatment)
  distr <- replicate(iters, {
    tr <- sample(treatment)
    tmp <- cor(resid, tr)
    sign(tmp+shift)*min(abs(tmp)+shift, 1)
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



#
# #' Test the correlation between residuals and treatment
# #'
# #' Carry out a permutation test for the Pearson correlation between residuals and treatment.
# #' @param prediction Vector of predicted outcomes
# #' @param response Vector of responses
# #' @param treatment Vector of treatments
# #' @param iters The number of Monte Carlo iterations (default 1000)
# #' @param shift Null shift away from 0 (a scalar constant) in correlation (default 0)
# #' @return a list containing the results: attributes estimate (the estimated difference), distr (simulated permutation distribution), and pvalue (p-value for the test)
# permu_pearson_shift <- function(prediction, response, treatment, iters = 1000, shift = 0){
#   resid <- prediction-response
#   pearson_r <- cor(resid, treatment)
#   distr <- replicate(iters, {
#     tr <- sample(treatment)
#     cor(resid, tr) + shift
#   })
#   pval <- c("p_upper" = sum(distr >= pearson_r)/iters,
#             "p_lower" = sum(distr <= pearson_r)/iters,
#             "twosided" = sum(abs(distr) >= abs(pearson_r))/iters)
#   return(list("estimate" = pearson_r, "distr" = distr, "pvalue"=pval))
# }
#


#' Confidence interval for the correlation between residuals and treatment
#'
#' Invert the permutation test to get a \eqn{1-\alpha} confidence interval for the Pearson correlation between residuals and treatment
#' @inheritParams permu_pearson
#' @param side Type of interval, either "both", "upper", or "lower". Default is "both".
#' @param alpha Significance level
#' @param verbose Verbosity switch - print the p-value and confidence interval endpoint at each step? (Default FALSE)
#' @return a confidence interval (vector)
permu_CI_pearson <- function(prediction, response, treatment, iters = 1000, alpha = 0.05, side = "both", verbosity = FALSE){
  resid <- prediction-response

  if(side == "both"){
    d1 <- permu_CI_pearson(prediction, response, treatment, iters, alpha = alpha/2, side = "lower", verbosity = verbosity)
    d2 <- permu_CI_pearson(prediction, response, treatment, iters, alpha = alpha/2, side = "upper", verbosity = verbosity)
    return(c(d1,d2))
  }
  pvalue_side <- which(c("lower", "upper", "both") %in% side); if(verbosity){print(side)}
  res <- permu_pearson(prediction, response, treatment, iters=1)
  init_estimate <- res$estimate

  shift_permtest <- function(ss){
    # if pval = TRUE, return ONLY the relevant p-value
    res <- permu_pearson(prediction, treatment, response, iters = iters, shift = ss)
    return(res$pvalue[pvalue_side])
  }

  # Bisection
  if(side == "upper"){
    l_int <- init_estimate; u_int <- 1
  }else{
    l_int <- -1; u_int <- init_estimate
  }
  d <- uniroot(function(x) {shift_permtest(x)-alpha}, lower = l_int, upper = u_int, tol = alpha/iters)
  return(d$root)
}

