library(ModelMatch)
library(ebal)
library(dplyr)
library(ggplot2)
library(reshape2)
library(VGAM)
library(grid)
library(gridExtra)
library(MASS)

### Helper functions to compute things

compute_diffmeans <- function(matches, Y){
  sizes <- sapply(matches, nrow)
  weighted.mean(
    x = sapply(matches, function(x){
      t1 <- mean(Y[x$index[x$treatment==1]])
      if(is.nan(t1)){t1 <- 0}
      t2 <- mean(Y[x$index[x$treatment==0]])
      if(is.nan(t2)){t2 <- 0}
      t1-t2
  }),
  w = sizes/sum(sizes))
}

compute_rmse <- function(x, truth){
  x <- x[!is.nan(x)]
  sqrt(mean((x-truth)^2, na.rm=TRUE))
}

compute_power <- function(pvalues){
  sapply((0:99)/100, function(p) mean(pvalues <= p, na.rm = TRUE))
}


### Main simulation functions

simulate_estimation <- function(gamma, B, N, selection = "random", nu = 0.5, errors = "normal", beta = c(1,2,4)){
  # Run simulations
  # gamma     = the (constant additive) treatment effect
  # N         = number of individuals in the sample
  # B         = number of replications
  # selection = treatment assignment mechanism. Default is "random".
  #             Options: "random", "correlated", "misspecified pscore"
  # nu        = covariance between treatment and X1 if selection == "correlated" or "misspecified pscore". Default 0.5
  # errors    = type of errors in the response model. Default is "normal".
  #             Options: "normal" (N(0,1)), "heteroskedastic" (N(0, abs(X2))), "heavy" (Laplacian(0,1))
  # beta      = vector of linear data-generating process coefficients. Defaults 1, 2, 4

  beta0 <- beta[1]
  beta1 <- beta[2]
  beta2 <- beta[3]
  # Set storage - things to report
  # Compare model-based matching, propensity score matching, t-test, and max entropy balancing

  # estimate: mean(treated) - mean(control)
  estimate <- matrix(NA, nrow = B*length(gamma), ncol = 7)
  colnames(estimate) <- c("Raw Diff in Means", "OLS with Controls", "MM Coarse Strata", "MM Fine Strata",
                          "Residualized, No Strata", "Pscore Matching", "Entropy Balancing")
  estimate <- cbind(estimate, "Gamma" = rep(gamma, each = B))
  rownum <- 0

  for(g in gamma){
    print(paste("Gamma = ", g))
    for(b in 1:B){
      print(rownum)
      rownum <- rownum + 1

      # Generate Xs and epsilon
      X1 <- rnorm(N); X2 <- rnorm(N)

      if(errors == "normal"){
        epsilon <- rnorm(N)
      }else{
        if(errors == "heteroskedastic"){
          epsilon <- rnorm(N, sd = sqrt(X1^2 + X2^2))
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
#        tr <- 1*(rnorm(N) >= 0)
        tr <- sample(rep(c(0,1), length.out = N))
      }else{
        if(selection == "correlated"){
          # T = nu*X_1 + delta, so cov(T, X_1) = nu*var(X_1) = nu
#          tr <- 1*((nu*X1 + rnorm(N)) >= 0)
          x1_lessthan <- which(X1 <= sort(X1)[50])
          x1_greaterthan <- which(X1 > sort(X1)[50])
          tr <- rep(0, N)
          tr[x1_lessthan] <- sample(rep(c(0, 0, 0, 0, 1), length.out = 50))
          tr[x1_lessthan]
          tr[x1_greaterthan] <- sample(rep(c(1, 1, 1, 1, 0), length.out = 50))
        }else{
          if(selection == "misspecified pscore"){
            # T = nu*X_1 + X_2 + X_2^3 + delta, so cov(T, X_1) = nu still
            tr <- 1*((nu*X1 + X2 + I(X_2^3) + rnorm(N)) >= 0)
          }else{
            stop("Invalid selection input")
          }
        }
      }

      # Generate Y
      Y <- beta0 + beta1*X1 + beta2*X2 + g*tr + epsilon
      dat <- data.frame(Y, X1, X2, tr)

      # Estimate Yhat using controls in a training set
      train <- dat[1:floor(0.3*nrow(dat)),]
      test <- dat[-(1:floor(0.3*nrow(dat))),]
      mm_model <- lm(Y~X1+X2, train, subset = (tr==0))
      Yhat <- predict(mm_model, test)

      # Estimate propensity score
      print(cor(X1, tr))
      pscore_mod <- glm(tr~X1+X2, dat, family=binomial(link="probit"))
      pscore <- predict(pscore_mod, dat, type = "response")

      # Create matches/strata for estimation
      mm_matches <- Strata(treatment = test$tr, prediction = Yhat, strata = 20)
      pscore_matches <- Matches(treatment = tr, prediction = pscore)
      mm_strata <- Strata(treatment = test$tr, prediction = Yhat, strata = 5)
      mm_nostrata <- Strata(treatment = test$tr, prediction = Yhat, strata = 1)

      # Create entropy balanced groups
      eb.out <- try(ebalance(Treatment = tr, X = cbind(X1, X2)),
                    silent=FALSE)
      if(class(eb.out)=="try-error"){
        cat("EB did not converge in this simulation\n")
      }else{
        estimate[rownum, "Entropy Balancing"] <- mean(Y[tr==1]) - weighted.mean(Y[tr==0], w=eb.out$w)
      }

      # Estimates
      estimate[rownum, "MM Fine Strata"] <- compute_diffmeans(mm_matches, test$Y-Yhat)
      estimate[rownum, "MM Coarse Strata"] <- compute_diffmeans(mm_strata, test$Y-Yhat)
      estimate[rownum, "Residualized, No Strata"] <- compute_diffmeans(mm_nostrata, test$Y-Yhat)

      estimate[rownum, "Pscore Matching"] <- compute_diffmeans(pscore_matches, Y)

      estimate[rownum, "Raw Diff in Means"] <- mean(Y[tr==1]) -  mean(Y[tr==0])
      estimate[rownum, "OLS with Controls"] <- lm(Y~X1 + X2 + tr, dat)$coeff["tr"]
    }
  }
  return(as.data.frame(estimate))
}


simulate_estimation_nonconstant <- function(gamma, B, N, selection = "random", nu = 0.5, errors = "normal", beta = c(1,2,4)){
  # Run simulations
  # gamma     = the (constant additive) treatment effect
  # N         = number of individuals in the sample
  # B         = number of replications
  # selection = treatment assignment mechanism. Default is "random".
  #             Options: "random", "correlated", "misspecified pscore"
  # nu        = covariance between treatment and X1 if selection == "correlated" or "misspecified pscore". Default 0.5
  # errors    = type of errors in the response model. Default is "normal".
  #             Options: "normal" (N(0,1)), "heteroskedastic" (N(0, abs(X2))), "heavy" (Laplacian(0,1))
  # beta      = vector of linear data-generating process coefficients. Defaults 1, 2, 4

  beta0 <- beta[1]
  beta1 <- beta[2]
  beta2 <- beta[3]
  # Set storage - things to report
  # Compare model-based matching, propensity score matching, t-test, and max entropy balancing

  # estimate: mean(treated) - mean(control)
  estimate <- matrix(NA, nrow = B*length(gamma), ncol = 7)
  colnames(estimate) <- c("Raw Diff in Means", "OLS with Controls", "MM Coarse Strata", "MM Fine Strata",
                          "Residualized, No Strata", "Pscore Matching", "Entropy Balancing")
  estimate <- cbind(estimate, "Gamma" = rep(gamma, each = B))
  rownum <- 0

  for(g in gamma){
    print(paste("Gamma = ", g))
    for(b in 1:B){
      print(rownum)
      rownum <- rownum + 1

      # Generate Xs and epsilon
      X1 <- rnorm(N); X2 <- rnorm(N)

      if(errors == "normal"){
        epsilon <- rnorm(N)
      }else{
        if(errors == "heteroskedastic"){
          epsilon <- rnorm(N, sd = sqrt(X1^2 + X2^2))
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
        #        tr <- 1*(rnorm(N) >= 0)
        tr <- sample(rep(c(0,1), length.out = N))
      }else{
        if(selection == "correlated"){
          # T = nu*X_1 + delta, so cov(T, X_1) = nu*var(X_1) = nu
          #          tr <- 1*((nu*X1 + rnorm(N)) >= 0)
          x1_lessthan <- which(X1 <= sort(X1)[50])
          x1_greaterthan <- which(X1 > sort(X1)[50])
          tr <- rep(0, N)
          tr[x1_lessthan] <- sample(rep(c(0, 0, 0, 0, 1), length.out = 50))
          tr[x1_lessthan]
          tr[x1_greaterthan] <- sample(rep(c(1, 1, 1, 1, 0), length.out = 50))
        }else{
          if(selection == "misspecified pscore"){
            # T = nu*X_1 + X_2 + X_2^3 + delta, so cov(T, X_1) = nu still
            tr <- 1*((nu*X1 + X2 + I(X_2^3) + rnorm(N)) >= 0)
          }else{
            stop("Invalid selection input")
          }
        }
      }

      # Generate Y
      Y <- beta0 + beta1*X1 + beta2*X2 + g*(1 + X1 + X2)*tr + epsilon
      dat <- data.frame(Y, X1, X2, tr)

      # Estimate Yhat using controls in a training set
      train <- dat[1:floor(0.3*nrow(dat)),]
      test <- dat[-(1:floor(0.3*nrow(dat))),]
      mm_model <- lm(Y~X1+X2, train, subset = (tr==0))
      Yhat <- predict(mm_model, test)

      # Estimate propensity score
      print(cor(X1, tr))
      pscore_mod <- glm(tr~X1+X2, dat, family=binomial(link="probit"))
      pscore <- predict(pscore_mod, dat, type = "response")

      # Create matches/strata for estimation
      mm_matches <- Strata(treatment = test$tr, prediction = Yhat, strata = 20)
      pscore_matches <- Matches(treatment = tr, prediction = pscore)
      mm_strata <- Strata(treatment = test$tr, prediction = Yhat, strata = 5)
      mm_nostrata <- Strata(treatment = test$tr, prediction = Yhat, strata = 1)

      # Create entropy balanced groups
      eb.out <- try(ebalance(Treatment = tr, X = cbind(X1, X2)),
                    silent=FALSE)
      if(class(eb.out)=="try-error"){
        cat("EB did not converge in this simulation\n")
      }else{
        estimate[rownum, "Entropy Balancing"] <- mean(Y[tr==1]) - weighted.mean(Y[tr==0], w=eb.out$w)
      }

      # Estimates
      estimate[rownum, "MM Fine Strata"] <- compute_diffmeans(mm_matches, test$Y-Yhat)
      estimate[rownum, "MM Coarse Strata"] <- compute_diffmeans(mm_strata, test$Y-Yhat)
      estimate[rownum, "Residualized, No Strata"] <- compute_diffmeans(mm_nostrata, test$Y-Yhat)

      estimate[rownum, "Pscore Matching"] <- compute_diffmeans(pscore_matches, Y)

      estimate[rownum, "Raw Diff in Means"] <- mean(Y[tr==1]) -  mean(Y[tr==0])
      estimate[rownum, "OLS with Controls"] <- lm(Y~X1 + X2 + tr + X1*tr + X2*tr, dat)$coeff["tr"]
    }
  }
  return(as.data.frame(estimate))
}


simulate_estimation_vary_mm <- function(gamma, B, N, selection = "random", nu = 0.5, errors = "normal", beta = c(1,2,4)){
  # Run simulations
  # gamma     = the (constant additive) treatment effect
  # N         = number of individuals in the sample
  # B         = number of replications
  # selection = treatment assignment mechanism. Default is "random".
  #             Options: "random", "correlated", "perfectly correlated", "misspecified pscore"
  # nu        = covariance between treatment and X1 if selection == "correlated" or "misspecified pscore". Default 0.5
  # errors    = type of errors in the response model. Default is "normal".
  #             Options: "normal" (N(0,1)), "heteroskedastic" (N(0, abs(X2))), "heavy" (Laplacian(0,1))
  # beta      = vector of linear data-generating process coefficients. Defaults 1, 2, 4

  beta0 <- beta[1]
  beta1 <- beta[2]
  beta2 <- beta[3]
  # Set storage - things to report
  # Compare model-based matching, propensity score matching, t-test, and max entropy balancing

  # estimate: mean(treated) - mean(control)
  estimate <- matrix(NA, nrow = B*length(gamma), ncol = 7)
  colnames(estimate) <- c("Pairs", "5 Strata", "Unstratified", "Pairs (ctrls)", "5 Strata (ctrls)", "Unstratified (ctrls)", "Unstratified (ctrls, regularized)")
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
          epsilon <- rnorm(N, sd = sqrt(X1^2 + X2^2))
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
          # T = nu*X_1 + delta, so cov(T, X_1) = nu*var(X_1)
          tr <- 1*((nu*X1 + rnorm(N)) >= 0)
        }else{
          if(selection == "perfectly correlated"){
            tr <- 1*(X1 >= 0)
          }else{
            if(selection == "misspecified pscore"){
            # T = nu*X_1 + X_1*X_2 + delta, so cov(T, X_1) = nu*var(X_1) still
              tr <- 1*((nu*X1 + X1*X2 + rnorm(N)) >= 0)
            }else{
              stop("Invalid selection input")
          }
          }
        }
      }

      # Generate Y
      Y <- beta0 + beta1*X1 + beta2*X2 + g*tr + epsilon
      dat <- data.frame(Y, X1, X2, tr)

      # Estimate Yhat using controls
      mm_model_ctrls <- lm(Y~X1+X2, dat, subset = (tr==0))
      Yhat_ctrls <- predict(mm_model_ctrls, dat)
      mm_model_ctrls <- lm.ridge(Y~X1+X2, data = dat, lambda = 2)
      Yhat_ctrls_reg <- scale(dat[,2:3],center = mm_model_ctrls$xm, scale = mm_model_ctrls$scales) %*% mm_model_ctrls$coef + mm_model_ctrls$ym

      # Estimate Yhat using all observations
      mm_model_all <- lm(Y~X1+X2, dat)
      Yhat_all <- predict(mm_model_all, dat)

      # Create matches/strata for estimation
      mm_matches_ctrls <- Matches(treatment = tr, prediction = Yhat_ctrls)
      mm_matches_all <- Matches(treatment = tr, prediction = Yhat_all)
      mm_strata_ctrls <- Strata(treatment = tr, prediction = Yhat_ctrls, strata = 5)
      mm_strata_all <- Strata(treatment = tr, prediction = Yhat_all, strata = 5)
      mm_unstr_ctrls <- Strata(treatment = tr, prediction = Yhat_ctrls, strata = 1)
      mm_unstr_all <- Strata(treatment = tr, prediction = Yhat_all, strata = 1)
      mm_unstr_reg <- Strata(treatment = tr, prediction = Yhat_ctrls_reg, strata = 1)

      # Estimates
      estimate[rownum, "Pairs"] <- compute_diffmeans(mm_matches_all, Y-Yhat_all)
      estimate[rownum, "5 Strata"] <- compute_diffmeans(mm_strata_all, Y-Yhat_all)

      estimate[rownum, "Pairs (ctrls)"] <- compute_diffmeans(mm_matches_ctrls, Y-Yhat_ctrls)
      estimate[rownum, "5 Strata (ctrls)"] <- compute_diffmeans(mm_strata_ctrls, Y-Yhat_ctrls)

      estimate[rownum, "Unstratified"] <- compute_diffmeans(mm_unstr_all, Y-Yhat_all)
      estimate[rownum, "Unstratified (ctrls)"] <- compute_diffmeans(mm_unstr_ctrls, Y-Yhat_ctrls)
      estimate[rownum, "Unstratified (ctrls, regularized)"] <- compute_diffmeans(mm_unstr_reg, Y-Yhat_ctrls_reg)

    }
  }
  return(as.data.frame(estimate))
}


simulate_tests <- function(gamma, N, B = 100, Bx = 100, selection = "random", nu = 0.5, errors = "normal"){
  # Run simulations
  # gamma     = the (constant additive) treatment effect
  # N         = number of individuals in the sample
  # B         = number of replications for each X draw
  # Bx        = number of random X draws
  # selection = treatment assignment mechanism. Default is "random".
  #             Options: "random", "correlated", "perfectly correlated", "misspecified pscore"
  # nu        = covariance between treatment and X1 if selection == "correlated" or "misspecified pscore". Default 0.5
  # errors    = type of errors in the response model. Default is "normal".
  #             Options: "normal" (N(0,1)), "heteroskedastic" (N(0, abs(X2))), "heavy" (Laplacian(0,1))

  beta0 <- 1
  beta1 <- 2
  beta2 <- 4
  # Set storage - things to report
  # Compare model-based matching, t-test after OLS, and randomization without controlling for covariates

  # power: number of times pvalue < 0.05
  pvalue <- matrix(NA, nrow = B*Bx*length(gamma), ncol = 5)
  colnames(pvalue) <- c("MM Fine Strata", "MM Coarse Strata","Residualized, No Strata","Wilcoxon", "OLS")
  pvalue <- cbind(pvalue, "Gamma" = rep(gamma, each = B))
  rownum <- 0

  for(g in gamma){
    print(paste("Gamma = ", g))
    for(bx in 1:Bx){
      X1 <- rnorm(N); X2 <- rnorm(N)

      for(b in 1:B){
        print(rownum)
        rownum <- rownum + 1

        # Generate epsilon

        if(errors == "normal"){
          epsilon <- rnorm(N)
        }else{
          if(errors == "heteroskedastic"){
            epsilon <- rnorm(N, sd = sqrt(X1^2 + X2^2))
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
          #        tr <- 1*(rnorm(N) >= 0)
          tr <- sample(rep(c(0,1), length.out = N))
        }else{
          if(selection == "correlated"){
            # T = nu*X_1 + delta, so cov(T, X_1) = nu*var(X_1) = nu
            #          tr <- 1*((nu*X1 + rnorm(N)) >= 0)
            x1_lessthan <- which(X1 <= sort(X1)[50])
            x1_greaterthan <- which(X1 > sort(X1)[50])
            tr <- rep(0, N)
            tr[x1_lessthan] <- sample(rep(c(0, 0, 0, 0, 1), length.out = 50))
            tr[x1_lessthan]
            tr[x1_greaterthan] <- sample(rep(c(1, 1, 1, 1, 0), length.out = 50))
          }else{
            if(selection == "perfectly correlated"){
              tr <- 1*(X1 >= 0)
            }else{
              if(selection == "misspecified pscore"){
                # T = nu*X_1 + X_1*X_2 + delta, so cov(T, X_1) = nu*var(X_1) still
                tr <- 1*((nu*X1 + X1*X2 + rnorm(N)) >= 0)
              }else{
                stop("Invalid selection input")
              }
            }
          }
        }

        # Generate Y
        Y <- beta0 + beta1*X1 + beta2*X2 + g*tr + epsilon
        dat <- data.frame(Y, X1, X2, tr)
        ntrainset <- floor(0.3*nrow(dat))
        train <- dat[1:ntrainset, ]
        test <- dat[(ntrainset+1):nrow(dat), ]

        # Estimate Yhat using controls
        mm_model <- lm(Y~X1+X2, train, subset = (tr==0))
        Yhat <- predict(mm_model, test)

        # Create matches/strata for estimation
        mm_matches <- Strata(treatment = test$tr, prediction = Yhat, strata = 20)
        mm_strata <- Strata(treatment = test$tr, prediction = Yhat, strata = 5)
        mm_nostrata <- Strata(treatment = test$tr, prediction = Yhat, strata = 1)

        # Tests
        strata_test <- permu_test_mean(strata = mm_strata, prediction = Yhat, treatment = test$tr,
                                       response = test$Y)
        fine_strata_test <- permu_test_mean(strata = mm_matches, prediction = Yhat, treatment = test$tr,
                                            response = test$Y)
        no_strata_test <- permu_test_mean(strata = mm_nostrata, prediction = Yhat, treatment = test$tr,
                                          response = test$Y)

        # Estimates
        pvalue[rownum, "MM Fine Strata"] <- fine_strata_test$pvalue["p_upper"]
        pvalue[rownum, "MM Coarse Strata"] <- strata_test$pvalue["p_upper"]
        pvalue[rownum, "Residualized, No Strata"] <- no_strata_test$pvalue["p_upper"]
        pvalue[rownum, "Wilcoxon"] <- wilcox.test(Y[tr == 1], Y[tr == 0])$p.value
        pvalue[rownum, "OLS"] <- summary(lm(Y~., dat))$coeff["tr", "Pr(>|t|)"]
      }
    }
  }
  return(as.data.frame(pvalue))
}


simulate_tests_nonconstant <- function(gamma, B, N, selection = "random", nu = 0.5, errors = "normal"){
  # Run simulations
  # gamma     = the (constant additive) treatment effect
  # N         = number of individuals in the sample
  # B         = number of replications for each X draw
  # selection = treatment assignment mechanism. Default is "random".
  #             Options: "random", "correlated", "perfectly correlated", "misspecified pscore"
  # nu        = covariance between treatment and X1 if selection == "correlated" or "misspecified pscore". Default 0.5
  # errors    = type of errors in the response model. Default is "normal".
  #             Options: "normal" (N(0,1)), "heteroskedastic" (N(0, abs(X2))), "heavy" (Laplacian(0,1))

  beta0 <- 1
  beta1 <- 2
  beta2 <- 4
  # Set storage - things to report
  # Compare model-based matching, t-test after OLS, and randomization without controlling for covariates

  # power: number of times pvalue < 0.05
  pvalue <- matrix(NA, nrow = B*Bx*length(gamma), ncol = 6)
  colnames(pvalue) <- c("MM Fine Strata", "MM Coarse Strata","Residualized, No Strata","Wilcoxon", "OLS", "OLS with Interactions")
  pvalue <- cbind(pvalue, "Gamma" = rep(gamma, each = B))
  rownum <- 0
  Bx <- 1

  for(g in gamma){
    for(bx in 1:Bx){
      X1 <- rnorm(N); X2 <- rnorm(N)

      for(b in 1:B){
        print(rownum)
        rownum <- rownum + 1

        # Generate epsilon

        if(errors == "normal"){
          epsilon <- rnorm(N)
        }else{
          if(errors == "heteroskedastic"){
            epsilon <- rnorm(N, sd = sqrt(X1^2 + X2^2))
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
          #        tr <- 1*(rnorm(N) >= 0)
          tr <- sample(rep(c(0,1), length.out = N))
        }else{
          if(selection == "correlated"){
            # T = nu*X_1 + delta, so cov(T, X_1) = nu*var(X_1) = nu
            #          tr <- 1*((nu*X1 + rnorm(N)) >= 0)
            x1_lessthan <- which(X1 <= sort(X1)[50])
            x1_greaterthan <- which(X1 > sort(X1)[50])
            tr <- rep(0, N)
            tr[x1_lessthan] <- sample(rep(c(0, 0, 0, 0, 1), length.out = 50))
            tr[x1_lessthan]
            tr[x1_greaterthan] <- sample(rep(c(1, 1, 1, 1, 0), length.out = 50))
          }else{
            if(selection == "perfectly correlated"){
              tr <- 1*(X1 >= 0)
            }else{
              if(selection == "misspecified pscore"){
                # T = nu*X_1 + X_1*X_2 + delta, so cov(T, X_1) = nu*var(X_1) still
                tr <- 1*((nu*X1 + X1*X2 + rnorm(N)) >= 0)
              }else{
                stop("Invalid selection input")
              }
            }
          }
        }

        # Generate Y
        Y <- beta0 + beta1*X1 + beta2*X2 + g*(1 + X1 + X2)*tr + epsilon
        dat <- data.frame(Y, X1, X2, tr)
        ntrainset <- floor(0.3*nrow(dat))
        train <- dat[1:ntrainset, ]
        test <- dat[(ntrainset+1):nrow(dat), ]

        # Estimate Yhat using controls
        mm_model <- lm(Y~X1+X2, train, subset = (tr==0))
        Yhat <- predict(mm_model, test)

        # Create matches/strata for estimation
        mm_matches <- Strata(treatment = test$tr, prediction = Yhat, strata = 20)
        mm_strata <- Strata(treatment = test$tr, prediction = Yhat, strata = 5)
        mm_nostrata <- Strata(treatment = test$tr, prediction = Yhat, strata = 1)

        # Tests
        strata_test <- permu_test_mean(strata = mm_strata, prediction = Yhat, treatment = test$tr,
                                       response = test$Y)
        fine_strata_test <- permu_test_mean(strata = mm_matches, prediction = Yhat, treatment = test$tr,
                                            response = test$Y)
        no_strata_test <- permu_test_mean(strata = mm_nostrata, prediction = Yhat, treatment = test$tr,
                                          response = test$Y)

        # Estimates
        pvalue[rownum, "MM Fine Strata"] <- fine_strata_test$pvalue["p_upper"]
        pvalue[rownum, "MM Coarse Strata"] <- strata_test$pvalue["p_upper"]
        pvalue[rownum, "Residualized, No Strata"] <- no_strata_test$pvalue["p_upper"]
        pvalue[rownum, "Wilcoxon"] <- wilcox.test(Y[tr == 1], Y[tr == 0])$p.value
        pvalue[rownum, "OLS"] <- summary(lm(Y~., dat))$coeff["tr", "Pr(>|t|)"]
        pvalue[rownum, "OLS with Interactions"] <- summary(lm(Y~X1 + X2 + tr + X1:tr + X2:tr, dat))$coeff["tr", "Pr(>|t|)"]
      }
    }
  }
return(as.data.frame(pvalue))
}


simulate_tests_which_residuals <- function(gamma, B, N, selection = "random", nu = 0.5, errors = "normal",
                                           refit_method = TRUE){
  # Run simulations, varying which model the residuals come from
  # gamma     = the (constant additive) treatment effect
  # N         = number of individuals in the sample
  # B         = number of replications
  # selection = treatment assignment mechanism. Default is "random".
  #             Options: "random", "correlated", "misspecified pscore"
  # nu        = covariance between treatment and X1 if selection == "correlated" or "misspecified pscore". Default 0.5
  # errors    = type of errors in the response model. Default is "normal".
  #             Options: "normal" (N(0,1)), "heteroskedastic" (N(0, abs(X2))), "heavy" (Laplacian(0,1))
  # refit_method = flag whether or not to run the method of refitting to controls at each permutation

  beta0 <- 1
  beta1 <- 2
  beta2 <- 4
  # Set storage - things to report
  # Compare model-based matching, t-test after OLS, and randomization without controlling for covariates

  # power: number of times pvalue < 0.05
  pvalue <- matrix(NA, nrow = B*length(gamma), ncol = 6)
  colnames(pvalue) <- c("MM fit to controls", "MM fit to controls, regularized",
                        "MM fit to treated","MM fit to all",
                        "MM with combined residuals",
                        "Match on Y, refit Yhat each permutation")
  pvalue <- cbind(pvalue, "Gamma" = rep(gamma, each = B))
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
          epsilon <- rnorm(N, sd = sqrt(X1^2 + X2^2))
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
          # T = nu*X_1 + delta, so cov(T, X_1) = nu*var(X_1)
          tr <- 1*((nu*X1 + rnorm(N)) >= 0)
        }else{
          if(selection == "misspecified pscore"){
            # T = nu*X_1 + X_1*X_2 + delta, so cov(T, X_1) = nu*var(X_1) still
            tr <- 1*((nu*X1 + X1*X2 + rnorm(N)) >= 0)
          }else{
            stop("Invalid selection input")
          }
        }
      }

      # Generate Y
      Y <- beta0 + beta1*X1 + beta2*X2 + g*tr + epsilon
      dat <- data.frame(Y, X1, X2, tr)

      # Estimate Yhat using controls
      mm_model_ctrls <- lm(Y~X1+X2, dat, subset = (tr==0))
      Yhat_ctrls <- predict(mm_model_ctrls, dat)
      mm_model_ctrls <- lm.ridge(Y~X1+X2, data = dat, lambda = 2)
      Yhat_ctrls_reg <- scale(dat[,2:3],center = mm_model_ctrls$xm, scale = mm_model_ctrls$scales) %*% mm_model_ctrls$coef + mm_model_ctrls$ym

      # Estimate Yhat using treated
      mm_model_tr <- lm(Y~X1+X2, dat, subset = (tr==1))
      Yhat_tr <- predict(mm_model_tr, dat)
      # Estimate Yhat on all data
      mm_model <- lm(Y~X1+X2, dat)
      Yhat <- predict(mm_model, dat)


      # Create matches/strata for estimation
      mm_strata_ctrls <- Strata(treatment = tr, prediction = Yhat_ctrls, strata = 5)
      mm_strata_ctrls_reg <- Strata(treatment = tr, prediction = Yhat_ctrls_reg, strata = 5)
      mm_strata_tr <- Strata(treatment = tr, prediction = Yhat_tr, strata = 5)
      mm_strata_all <- Strata(treatment = tr, prediction = Yhat, strata = 5)
      mm_strata_addresiduals <- Strata(treatment = tr, prediction = (Yhat_ctrls + Yhat_tr)/2,
                                       strata = 5)
      mm_strata_Y <- Strata(treatment = tr, prediction = Y, strata = 5)

      # Tests
      mm_ctrl_test <- permu_test_mean(strata = mm_strata_ctrls, prediction = Yhat_ctrls, treatment = tr,
                                          response = Y)
      mm_ctrl_reg_test <- permu_test_mean(strata = mm_strata_ctrls_reg, prediction = Yhat_ctrls_reg, treatment = tr,
                                          response = Y)
      mm_tr_test <- permu_test_mean(strata = mm_strata_tr, prediction = Yhat_tr, treatment = tr,
                                      response = Y)
      mm_all_test <- permu_test_mean(strata = mm_strata_all, prediction = Yhat, treatment = tr,
                                      response = Y)
      mm_addresiduals_test <- permu_test_mean(strata = mm_strata_addresiduals, prediction = (Yhat_ctrls + Yhat_tr)/2,
                                              treatment = tr, response = Y)

      if(refit_method == TRUE){
        mm_refit <- permu_test_mean_refit(strata = mm_strata_Y, response = Y, X = data.frame(X1, X2),
                                        fit_function = function(X, Y) lm(Y~., data.frame(X, Y)))
        pvalue[rownum, "Match on Y, refit Yhat each permutation"] <- mm_refit$pvalue["p_upper"]
      } else {
        pvalue[rownum, "Match on Y, refit Yhat each permutation"] <- NA
      }
      # Estimates
      pvalue[rownum, "MM fit to controls"] <- mm_ctrl_test$pvalue["p_upper"]
      pvalue[rownum, "MM fit to controls, regularized"] <- mm_ctrl_reg_test$pvalue["p_upper"]
      pvalue[rownum, "MM fit to treated"] <- mm_tr_test$pvalue["p_upper"]
      pvalue[rownum, "MM fit to all"] <- mm_all_test$pvalue["p_upper"]
      pvalue[rownum, "MM with combined residuals"] <- mm_addresiduals_test$pvalue["p_upper"]
    }
  }
  return(as.data.frame(pvalue))
}
### Plots

plot_est_by_gamma <- function(estimates, color_mm = FALSE, fill_by_fitmethod = FALSE){
  # Input ``estimates'' should be the output of simulate_estimation
  res_plot <- melt(estimates, id.vars = "Gamma")
  res_plot$color_mm = rep("black", nrow(res_plot))
  if(color_mm == TRUE){
    res_plot$color_mm[grepl("MM", res_plot$variable)] <- "#629e1f"
  }
  if(fill_by_fitmethod == TRUE){
    res_plot$color_mm <- "#DB7093"
    res_plot$color_mm[grepl("ctrls", res_plot$variable)] <- "#629e1f"
    res_plot$color_mm[grepl("regularized", res_plot$variable)] <- "black"
  }
  p <- ggplot(res_plot, aes(x = variable, y = value)) +
    geom_boxplot(aes(fill = color_mm), alpha = 0.25) +
    facet_wrap(~Gamma) +
    geom_hline(aes(yintercept = Gamma), linetype = "dashed") +
    xlab("Estimation Method") +
    ylab("Estimate") +
    scale_fill_manual(c("black", "#629e1f", "#DB7093"), values=c("black", "#629e1f", "#DB7093"), guide = FALSE)
}

plot_power_curves <- function(pvalues){
  # Input ``pvalues'' should be the output of simulate_tests
  gamma <- unique(pvalues[,"Gamma"])
  power_curves <- lapply(gamma, function(g){
    gamma_subset <- pvalues[pvalues$Gamma == g, -5]
    apply(gamma_subset, 2, compute_power)
  })
  power_curves <- do.call(rbind, power_curves)
  power_curves <- as.data.frame(cbind(power_curves,
                                      "alpha" = rep((0:99)/100, length(gamma)),
                                      "gamma" = rep(unique(as.character(pvalues$Gamma)), each = 100)
                                ))
  power_curves$gamma <- factor(power_curves$gamma, levels = unique(as.character(pvalues$Gamma)))
  power_curves[,"MM Coarse Strata"] <- as.numeric(as.character(power_curves[,"MM Coarse Strata"]))
  power_curves[,"MM Fine Strata"] <- as.numeric(as.character(power_curves[,"MM Fine Strata"]))
  power_curves[,"Wilcoxon"]      <- as.numeric(as.character(power_curves[,"Wilcoxon"]))
  power_curves[,"OLS"]           <- as.numeric(as.character(power_curves[,"OLS"]))
  power_curves[,"alpha"]         <- as.numeric(as.character(power_curves[,"alpha"]))
  power_curves_plot <- melt(power_curves, id.vars = c("alpha", "gamma"),
                            variable.name = "Method")
  ggplot(power_curves_plot, aes(x = alpha, y = value)) +
    geom_line(aes(color = Method)) +
    facet_wrap(~gamma) +
    xlab("Significance level") +
    ylab("Power")
  }

plot_power_curves_which_residuals <- function(pvalues){
  # Input ``pvalues'' should be the output of simulate_tests
  gamma <- unique(pvalues[,"Gamma"])
  refit_ind <- which(colnames(pvalues) == "Match on Y, refit Yhat each permutation")
  if(all(is.na(pvalues[, refit_ind]))){
    pvalues <- pvalues[,-refit_ind]
  }
  vv <- ncol(pvalues)
  power_curves <- lapply(gamma, function(g){
    gamma_subset <- pvalues[pvalues$Gamma == g, -vv]
    apply(gamma_subset, 2, compute_power)
  })
  power_curves <- do.call(rbind, power_curves)
  power_curves <- as.data.frame(cbind(power_curves,
                                      "alpha" = rep((0:99)/100, length(gamma)),
                                      "gamma" = rep(gamma, each = 100)
  ))
  power_curves[,"MM fit to controls"] <- as.numeric(as.character(power_curves[,"MM fit to controls"]))
  power_curves[,"MM fit to controls, regularized"] <- as.numeric(as.character(power_curves[,"MM fit to controls, regularized"]))
  power_curves[,"MM fit to treated"] <- as.numeric(as.character(power_curves[,"MM fit to treated"]))
  power_curves[,"MM fit to all"]      <- as.numeric(as.character(power_curves[,"MM fit to all"]))
  power_curves[,"MM with combined residuals"] <- as.numeric(as.character(power_curves[,"MM with combined residuals"]))
#  power_curves[,"MM with combined res and strat by ctrl"] <- as.numeric(as.character(power_curves[,"MM with combined res and strat by ctrl"]))
  if("Match on Y, refit Yhat each permutation" %in% colnames(power_curves)){
    power_curves[,"Match on Y, refit Yhat each permutation"] <- as.numeric(as.character(power_curves[,"Match on Y, refit Yhat each permutation"]))
  }
  power_curves[,"alpha"]         <- as.numeric(as.character(power_curves[,"alpha"]))
  power_curves_plot <- melt(power_curves, id.vars = c("alpha", "gamma"),
                            variable.name = "Method")
  ggplot(power_curves_plot, aes(x = alpha, y = value)) +
    geom_line(aes(color = Method)) +
    facet_wrap(~gamma) +
    xlab("Significance level") +
    ylab("Power") +
    theme(legend.position = "bottom")
}

plot_by_yhat_model <- function(gamma, N, selection = "random", nu = 0.5, errors = "normal"){
    # Run simulations
    # gamma     = the (constant additive) treatment effect
    # N         = number of individuals in the sample
    # B         = number of replications
    # selection = treatment assignment mechanism. Default is "random".
    #             Options: "random", "correlated", "misspecified pscore"
    # nu        = covariance between treatment and X1 if selection == "correlated" or "misspecified pscore". Default 0.5
    # errors    = type of errors in the response model. Default is "normal".
    #             Options: "normal" (N(0,1)), "heteroskedastic" (N(0, abs(X2))), "heavy" (Laplacian(0,1))

    beta0 <- 1
    beta1 <- 2
    beta2 <- 4
    # Set storage - things to report
    # Compare model-based matching, t-test after OLS, and randomization without controlling for covariates

    # Generate Xs and epsilon
    X1 <- rnorm(N); X2 <- rnorm(N, sd = 2)

    if(errors == "normal"){
      epsilon <- rnorm(N)
    }else{
      if(errors == "heteroskedastic"){
        epsilon <- rnorm(N, sd = sqrt(X1^2 + X2^2))
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
        # T = nu*X_1 + delta, so cov(T, X_1) = nu*var(X_1)
        tr <- 1*((nu*X1 + rnorm(N)) >= 0)
      }else{
        if(selection == "misspecified pscore"){
          # T = nu*X_1 + X_1*X_2 + delta, so cov(T, X_1) = nu*var(X_1) still
          tr <- 1*((nu*X1 + X1*X2 + rnorm(N)) >= 0)
        }else{
          stop("Invalid selection input")
        }
      }
    }

    # Generate Y
    Y <- beta0 + beta1*X1 + beta2*X2 + gamma*tr + epsilon
    dat <- data.frame(Y, X1, X2, tr)

    # Estimate Yhat using controls
    mm_model_ctrls <- lm(Y~X1+X2, dat, subset = (tr==0))
    Yhat_ctrls <- predict(mm_model_ctrls, dat)
    # Estimate Yhat on all data
    mm_model <- lm(Y~X1+X2, dat)
    Yhat <- predict(mm_model, dat)


    # Create matches/strata for estimation
    mm_strata <- Strata(treatment = tr, prediction = Yhat, strata = 5)
    mm_strata_ctrls <- Strata(treatment = tr, prediction = Yhat_ctrls, strata = 5)
    print(c(
      "All Obs"  = within_group_mean(mm_strata, Yhat, Y),
      "Controls" = within_group_mean(mm_strata_ctrls, Yhat_ctrls, Y)
    ))

    # Reshape for plotting
    stratum <- rep(1:length(mm_strata), sapply(mm_strata, nrow))
    mm_strata <- data.frame(stratum, "residual" = Yhat - Y, do.call(rbind, mm_strata))
    stratum <- rep(1:length(mm_strata_ctrls), sapply(mm_strata_ctrls, nrow))
    mm_strata_ctrls <- data.frame(stratum, "residual" = Yhat_ctrls - Y, do.call(rbind, mm_strata_ctrls))

    p1 <- ggplot(mm_strata, aes(residual)) + geom_histogram(aes(fill = factor(treatment)), binwidth = 0.5, position = "dodge") + facet_grid(~stratum) + ggtitle("Model fit to all observations")
    p2 <- ggplot(mm_strata_ctrls, aes(residual)) + geom_histogram(aes(fill = factor(treatment)), binwidth = 0.5, position = "dodge") + facet_grid(~stratum) + ggtitle("Model fit only to controls")

    grid.arrange(p1, p2)
}
