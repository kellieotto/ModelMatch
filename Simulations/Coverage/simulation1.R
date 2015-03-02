### Confidence interval simulation 1
### Based on JH Ebal simulations
### by Kellie Ottoboni; last edited 3/1/2015
### Simple design from Freedman (Weighting Regressions by P Scores)
### Y = a + bX + c1Z1 + c2Z2 + dU
### X = I(e + f1Z1 + f2Z2 + V > 0)
### U and V are N(0,1); U, V, (Z1, Z2) are mutually independent
### True treatment effect b is fixed at 1.
### Compare nominal (1-alpha)% coverage to true coverage of CIs


rm(list=ls())
library(Matching)
library(MASS)
library(ebal)
library(devtools)
install_github("kellieotto/ModelMatch/ModelMatch")
library(ModelMatch)
library(glmnet)
library(ipred)
library(randomForest)

## Setup
sims <- 100
alpha <- 0.05
B <- 500 # number of bootstrap replicates to obtain matching CIs

# (U, V, Z1, Z2) are multivariate normal with variances 1 and covars 0,0,.75,.75
a <- b <- c1 <- d <- 1
c2 <- 2
e <- 0.5
f1 <- 0.25
f2 <- 0.75
Sigma <- diag(rep(1,4))
Sigma[3,3] <- 2
Sigma[4,4] <- 1
Sigma[3,4] <- Sigma[4,3] <- 1
mu <- c(0, 0, 0.5, 1)



# pick target sample size
#nsim <- 1500
#nsim <- 600
nsim  <- 300
tr <- nsim/2; tr
co <- nsim - tr; co
est    <- c("RAW","MM1","MM2","MM3","MM4","PS1","PS2","PS3","PSW1","PSW2","PSW3","EB","GM")

# expand storage for three different outcomes
ci.store <- vector("list",length(est)); names(ci.store) <- est
for(name in est){ci.store[[name]] <- data.frame("lower" =rep(NA, sims), "upper"=rep(NA, sims))}
est.store <- matrix(NA, sims, length(est), dimnames = list(c(), est))

for(i in 1:sims){ # start simulations
  print(i)
  # draw Xs
  var <- mvrnorm(nsim, mu=mu, Sigma=Sigma, empirical = F)
  U <- var[,1]
  V <- var[,2]
  Z1 <- var[,3]
  Z2 <- var[,4]
  X <- as.numeric((e + f1*Z1 + f2*Z2 + V )> 0)
  Y <- a + b*X + c1*Z1 + c2*Z2 + d*U
  dat <- data.frame(Y, X, Z1, Z2)

  # true linear propensity scores
  PS.out  <- glm(X~ Z1 + Z2,family=binomial(link = "probit"))
  PS.true <- PS.out$fitted

  # propensity score specification formulae
  psxlist <- pslist <- list()
  # ps 1 - misspecified
  formula.PS1 <- formula(X ~ Z1)
  psxlist[[1]] <- model.matrix(formula.PS1)[,-1]
  # ps 2 - misspecified
  formula.PS2  <- formula(X ~ Z2)
  psxlist[[2]] <- model.matrix(formula.PS2)[,-1]
  # ps 3 - correctly specified
  formula.PS3  <- formula(X ~ Z1+Z2)
  psxlist[[3]] <- model.matrix(formula.PS3)[,-1]
  # PS estimation
  for(p in 1:3){
    pslist[[p]] <- glm(X~psxlist[[p]],family=binomial(link = "probit"))
  }


  # raw outcome
  ttest_ci <- confint(lm(Y~X), level = 1-alpha)[2,]
  ci.store[["RAW"]][i,] <- ttest_ci
  est.store[i,"RAW"] <- (b >= ttest_ci[1] & b<= ttest_ci[2])


  # Model matching predictions - predict on Z1 and Z2, NOT treatment var X
  mmlist <- list()
  # Model 1: simple linear regression on all covariates
  model.MM1 <- lm(Y~Z1 + Z2, data = dat)
  mmlist[[1]] <- model.MM1$fitted
  # Model 2: linear regression on all covariates with higher order polynomial terms
  model.MM2 <- lm(Y~Z1+I(Z1^2)+I(Z1^3)+Z2+I(Z2^2)+I(Z2^3), data = dat)
  mmlist[[2]] <- model.MM2$fitted
  # Model 3: bagged regression trees, default 25 bags
  model.MM3 <- bagging(Y~Z1 + Z2, data = dat, coob=TRUE)
  mmlist[[3]] <- predict(model.MM3, dat[,3:4])
  # Model 4: random forest
  model.MM4 <- randomForest(Y~Z1+Z2, data = dat)
  mmlist[[4]] <- predict(model.MM4, dat[,3:4])


  # model-based matching
  for(p in 1:4){
    pred <- mmlist[[p]]
    pairs <- Matches(X, pred)
    mm_res <- permu_CI_mean(groups = pairs, prediction = pred, response = Y, treatment = X, side = "both", alpha = alpha, iters=1000, precision = 0.25, verbose = F)
    ci.store[[paste("MM", p, sep="")]][i,] <- mm_res
    est.store[i,paste("MM", p, sep="")] <- (mm_res[[1]]<=1 & mm_res[[2]]>=1)
    }

#   ## GenMatching
#   dist_genmatch <- rep(NA, B)
#   for(b in 1:B){
#     dat_boot <- dat[sample(1:nrow(dat), replace=T),]
#     suppressWarnings(gout <- GenMatch(Tr=dat_boot$X, X=dat_boot[,3:4], BalanceMatrix=dat_boot[,3:4], estimand="ATT", M=1, weights=NULL,print.level=0))
#     out  <- Match(Tr=dat_boot$X,X=dat_boot[,3:4],Weight.matrix=gout,estimand="ATT",M=1,BiasAdjust=FALSE)
#     dist_genmatch[b] <- (sum(dat_boot$Y[out$index.treated]*out$weights)/sum(out$weights)) - (sum(dat_boot$Y[out$index.control]*out$weights)/sum(out$weights))
#   }
#   gm_ci <- quantile(dist_genmatch, na.rm=T, probs = c(alpha/2, 1-alpha/2))
#   est.store[i,"GM"] <- (gm_ci[1]<=1 & gm_ci[2]>=1)


  # Propensity scoring methods
  # Matching
  for(p in 1:3){
    dist_ps <- rep(NA, B)
    dist_psw <- rep(NA, B)
    if(sum(pslist[[p]]$fitted<(0+sqrt(.Machine$double.eps)))> 0 || sum(pslist[[p]]$fitted>(1-sqrt(.Machine$double.eps)))>0){ # error reporting if perfect sep
      cat("Perfect Seperation in glm occured in sim",i,"\n")
    } else {
      PS <- pslist[[p]]$linear
      for(bb in 1:B){
      boot <- sample(1:length(PS), length(PS), replace = TRUE)

      # Matching
      out <- Match(Tr=X[boot],X=PS[boot],estimand="ATT",M=1,BiasAdjust = FALSE)
      dist_ps[bb] <- (sum(Y[boot][out$index.treated]*out$weights)/sum(out$weights)) - (sum(Y[boot][out$index.control]*out$weights)/sum(out$weights))

      # Weighting
      PS.pr <- pslist[[p]]$fitted[boot]
      D <- PS.pr/(1-PS.pr)
      dist_psw[bb] <- mean(Y[boot][X==1])- (sum(Y[boot][X[boot]==0]*D[X[boot]==0])/sum(D[X[boot]==0]))
      }
    }
    ps_ci <- quantile(dist_ps, na.rm=T, probs = c(alpha/2, 1-alpha/2))
    ci.store[[paste("PS",p,sep="")]][i,] <- ps_ci
    est.store[i,paste("PS",p,sep="")] <- (ps_ci[1]<=1 & ps_ci[2]>=1)
    psw_ci <- quantile(dist_psw, na.rm=T, probs = c(alpha/2, 1-alpha/2))
    ci.store[[paste("PSW",p,sep="")]][i,] <- psw_ci
    est.store[i,paste("PSW",p,sep="")] <-(psw_ci[1]<=1 & psw_ci[2]>=1)
    }



  est.summ <- round(apply(est.store, 2, mean, na.rm=T), 5)
  print(est.summ)
}
ci.length.summ <- sapply(ci.store, function(x) x[,2] - x[,1])

save(est.store, est.summ, ci.store, ci.length.summ, file = "simulation1.Rdata")


