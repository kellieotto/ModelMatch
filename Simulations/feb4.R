library(ModelMatch)
library(ebal)
library(dplyr)
library(ggplot2)
library(reshape2)
library(xtable)

# Set parameters
gamma <- c(0.01, 0.1, 1, 10)
B <- 10
N <- 100

compute_diffmeans <- function(matches, Y){
  mean(sapply(matches, function(x){
    mean(Y[x$index[x$treatment==1]]) - mean(Y[x$index[x$treatment==0]])
  }))
}

compute_rmse <- function(x, truth){
  sqrt(mean((x-truth)^2))
}

simulate_estimation <- function(gamma, B, N, selection = "random", alpha = 0.5){
  # Run simulations
  # gamma     = the (constant additive) treatment effect
  # N         = number of individuals in the sample
  # B         = number of replications
  # selection = treatment assignment mechanism. Default is "random".
  #             Options: "random", "correlated", "misspecified pscore"
  # alpha     = covariance between treatment and X1 if selection == "correlated" or "misspecified pscore". Default 0.5

  beta0 <- 1
  beta1 <- 2
  beta2 <- 4
  # Set storage - things to report
  # Compare model-based matching, propensity score matching, t-test, and max entropy balancing

  # power: number of times pvalue < 0.05
  pvalue <- matrix(NA, nrow = B*length(gamma), ncol = 4)
  # estimate: mean(treated) - mean(control)
  estimate <- matrix(NA, nrow = B*length(gamma), ncol = 6)
  colnames(estimate) <- c("MM Pairs", "MM Strata", "Pscore Pairs", "Pscore Strata", "Ebal", "Unadjusted")
  estimate <- cbind(estimate, "Gamma" = rep(gamma, each = B))
  # covariance between T and X1, covariance between T and epsilon
  covx1 <- rep(NA, B*length(gamma))
  coveps <- rep(NA, B*length(gamma))
  rownum <- 0

  for(g in gamma){
    print(paste("Gamma = ", g))
    for(b in 1:B){
      rownum <- rownum + 1

      # Generate Xs and epsilon
      X1 <- rnorm(N); X2 <- rnorm(N, sd = 2)
      epsilon <- rnorm(N)

      # Treatment
      if(selection == "random"){
        tr <- 1*(rnorm(N) >= 0)
      }else{
        if(selection == "correlated"){
          # T = alpha*X_1 + delta, so cov(T, X_1) = alpha*var(X_1)
          tr <- 1*((alpha*X1 + rnorm(N)) >= 0)
        }else{
          if(selection == "misspecified pscore"){
            # T = alpha*X_1 + X_1*X_2 + delta, so cov(T, X_1) = alpha*var(X_1) still
            tr <- 1*((alpha*X1 + X1*X2 + rnorm(N)) >= 0)
          }
        }
      }

      # Generate Y
      Y <- beta0 + beta1*X1 + beta2*X2 + g*tr + epsilon
      dat <- data.frame(Y, X1, X2, tr)

      # Estimate Yhat using controls
      mm_model <- lm(Y~X1+X2, dat, subset = (tr==0))
      Yhat <- predict(mm_model, dat)

      # Estimate propensity score
      pscore_mod <- glm(tr~X1+X2, dat, family=binomial(link="probit"))
      pscore <- predict(pscore_mod, dat, type = "response")

      # Create matches/strata for estimation
      mm_matches <- Matches(treatment = tr, prediction = Yhat)
      pscore_matches <- Matches(treatment = tr, prediction = pscore)
      mm_strata <- Strata(treatment = tr, prediction = Yhat, strata = 5)
      pscore_strata <- Strata(treatment = tr, prediction = pscore, strata = 5)

      # Create entropy balanced groups
      eb.out <- ebalance(Treatment = tr, X = cbind(X1, X2))


      # Estimates
      estimate[rownum, "MM Pairs"] <- compute_diffmeans(mm_matches, Y-Yhat)
      estimate[rownum, "MM Strata"] <- compute_diffmeans(mm_strata, Y-Yhat)

      estimate[rownum, "Pscore Pairs"] <- compute_diffmeans(pscore_matches, Y)
      estimate[rownum, "Pscore Strata"] <- compute_diffmeans(pscore_strata, Y)

      estimate[rownum, "Ebal"] <- mean(Y[tr==1]) - weighted.mean(Y[tr==0], w=eb.out$w)
      estimate[rownum, "Unadjusted"] <- mean(Y[tr==1]) -  mean(Y[tr==0])
    }
  }
  return(as.data.frame(estimate))
}



res <- simulate_estimation(gamma, B = 10, N = 100, "random")
res_plot <- melt(res, id.vars = "Gamma")
ggplot(res_plot, aes(x = variable, y = value)) +
  geom_boxplot() +
  facet_wrap(~Gamma) +
  geom_hline(aes(yintercept = Gamma), linetype = "dashed") +
  xlab("Estimation Method") +
  ylab("Estimate") +
  ggtitle("Estimates of Varying Levels of Constant Additive Treatment Effects\n Random Treatment Assignment")
sapply(gamma, function(g){
  gamma_subset <- res[res$Gamma == g, -7]
  apply(gamma_subset, 2, compute_rmse, g)
})

res2 <- simulate_estimation(gamma, B = 10, N = 100, "correlated", alpha = 0.5)
res2_plot <- melt(res2, id.vars = "Gamma")
ggplot(res2_plot, aes(x = variable, y = value)) +
  geom_boxplot() +
  facet_wrap(~Gamma) +
  geom_hline(aes(yintercept = Gamma), linetype = "dashed") +
  xlab("Estimation Method") +
  ylab("Estimate") +
  ggtitle("Estimates of Varying Levels of Constant Additive Treatment Effects\n Treatment Assignment with Cov(T, X1) = 0.5")
sapply(gamma, function(g){
  gamma_subset <- res2[res2$Gamma == g, -7]
  apply(gamma_subset, 2, compute_rmse, g)
})

res3 <- simulate_estimation(gamma, B = 10, N = 100, "correlated", alpha = -0.9)
res3_plot <- melt(res3, id.vars = "Gamma")
ggplot(res3_plot, aes(x = variable, y = value)) +
  geom_boxplot() +
  facet_wrap(~Gamma) +
  geom_hline(aes(yintercept = Gamma), linetype = "dashed") +
  xlab("Estimation Method") +
  ylab("Estimate") +
  ggtitle("Estimates of Varying Levels of Constant Additive Treatment Effects\n Treatment Assignment with Cov(T, X1) = -0.9")
sapply(gamma, function(g){
  gamma_subset <- res3[res3$Gamma == g, -7]
  apply(gamma_subset, 2, compute_rmse, g)
})


res4 <- simulate_estimation(gamma, B = 10, N = 100, "misspecified pscore", alpha = 0.5)
res4_plot <- melt(res4, id.vars = "Gamma")
ggplot(res4_plot, aes(x = variable, y = value)) +
  geom_boxplot() +
  facet_wrap(~Gamma) +
  geom_hline(aes(yintercept = Gamma), linetype = "dashed") +
  xlab("Estimation Method") +
  ylab("Estimate") +
  ggtitle("Estimates of Varying Levels of Constant Additive Treatment Effects\n Treatment Assignment with Cov(T, X1) = 0.5\n Misspecified Propensity Score Model")
sapply(gamma, function(g){
  gamma_subset <- res4[res4$Gamma == g, -7]
  apply(gamma_subset, 2, compute_rmse, g)
})
