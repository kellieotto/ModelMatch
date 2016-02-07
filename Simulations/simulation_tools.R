library(ModelMatch)
library(ebal)
library(dplyr)
library(ggplot2)
library(reshape2)
library(VGAM)

### Helper functions to compute things

compute_diffmeans <- function(matches, Y){
  mean(sapply(matches, function(x){
    mean(Y[x$index[x$treatment==1]]) - mean(Y[x$index[x$treatment==0]])
  }))
}

compute_rmse <- function(x, truth){
  x <- x[!is.nan(x)]
  sqrt(mean((x-truth)^2, na.rm=TRUE))
}

### Main simulation function

simulate_estimation <- function(gamma, B, N, selection = "random", alpha = 0.5, errors = "normal"){
  # Run simulations
  # gamma     = the (constant additive) treatment effect
  # N         = number of individuals in the sample
  # B         = number of replications
  # selection = treatment assignment mechanism. Default is "random".
  #             Options: "random", "correlated", "misspecified pscore"
  # alpha     = covariance between treatment and X1 if selection == "correlated" or "misspecified pscore". Default 0.5
  # errors    = type of errors in the response model. Default is "normal".
  #             Options: "normal" (N(0,1)), "heteroskedastic" (N(0, abs(X2))), "heavy" (Laplacian(0,1))

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
  rownum <- 0

  for(g in gamma){
    print(paste("Gamma = ", g))
    for(b in 1:B){
      print(rownum)
      rownum <- rownum + 1

      # Generate Xs and epsilon
      X1 <- rnorm(N); X2 <- rnorm(N, sd = 2)

      if(errors == "normal"){
        epsilon <- rnorm(N)
      }else{
        if(errors == "heteroskedastic"){
          epsilon <- rnorm(N, sd = sqrt(abs(X2)))
        }else{
          if(errors == "heavy"){
            epsilon <- rlaplace(N)
          }else{
            stop("Invalid errors input")
          }
        }
      }
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
          }else{
            stop("Invalid selection input")
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
      print(cor(X1, tr))
      pscore_mod <- glm(tr~X1+X2, dat, family=binomial(link="probit"))
      pscore <- predict(pscore_mod, dat, type = "response")

      # Create matches/strata for estimation
      mm_matches <- Matches(treatment = tr, prediction = Yhat)
      pscore_matches <- Matches(treatment = tr, prediction = pscore)
      mm_strata <- Strata(treatment = tr, prediction = Yhat, strata = 5)
      pscore_strata <- Strata(treatment = tr, prediction = pscore, strata = 5)

      # Create entropy balanced groups
      eb.out <- try(ebalance(Treatment = tr, X = cbind(X1, X2)),
                    silent=FALSE)
      if(class(eb.out)=="try-error"){
        cat("EB did not converge in this simulation\n")
      }else{
        estimate[rownum, "Ebal"] <- mean(Y[tr==1]) - weighted.mean(Y[tr==0], w=eb.out$w)
      }

      # Estimates
      estimate[rownum, "MM Pairs"] <- compute_diffmeans(mm_matches, Y-Yhat)
      estimate[rownum, "MM Strata"] <- compute_diffmeans(mm_strata, Y-Yhat)

      estimate[rownum, "Pscore Pairs"] <- compute_diffmeans(pscore_matches, Y)
      estimate[rownum, "Pscore Strata"] <- compute_diffmeans(pscore_strata, Y)

      estimate[rownum, "Unadjusted"] <- mean(Y[tr==1]) -  mean(Y[tr==0])
    }
  }
  return(as.data.frame(estimate))
}

### Plots

plot_est_by_gamma <- function(estimates){
  # Input ``estimates'' should be the output of simulate_estimation
  res_plot <- melt(estimates, id.vars = "Gamma")
  p <- ggplot(res_plot, aes(x = variable, y = value)) +
    geom_boxplot() +
    facet_wrap(~Gamma) +
    geom_hline(aes(yintercept = Gamma), linetype = "dashed") +
    xlab("Estimation Method") +
    ylab("Estimate")
}
