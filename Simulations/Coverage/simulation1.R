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

# (U, V, Z1, Z2) are multivariate normal with variances 1 and covars 0,0,.75,.75
a <- b <- c1 <- d <- 1
c2 <- 2
e <- 0.5
f1 <- 0.25
f2 <- 0.75
sigma <- 1
rho <- 1
Sigma <- diag(rep(1,4))
Sigma[3,3] <- Sigma[4,4] <- sigma
Sigma[3,4] <- Sigma[4,3] <- rho
mu <- c(0, 0, 0.5, 1)



# pick target sample size
#nsim <- 1500
#nsim <- 600
nsim  <- 300
tr <- nsim/2; tr
co <- nsim - tr; co
est    <- c("RAW","MM1","MM1ATT","MM2","MM2ATT","MM3","MM3ATT","MM4","MM4ATT","PS1","PS2","PS3","PSW1","PSW2","PSW3","EB","GM")

# expand storage for three different outcomes
est.store <- matrix(NA, sims,length(est), dimnames = list(c(), est))

for(i in 1:sims){ # start simulations
  print(i)
  # draw Xs
  var <- mvrnorm(nsim, mu=mu, Sigma=Sigma, empirical = F)
  U <- var[,1]
  V <- var[,2]
  Z1 <- var[,3]
  Z2 <- var[,4]
  X <- (e + f1*Z1 + f2*Z2 + V > 0)
  Y <- a + b*X + c1*Z1 + c2*Z2 + d*U
  
#   # sample to target sizes
#   sim.dat <- data.frame(W,X1,X2,X3,X4,X5,X6)
#   tr.keep <- sample(which(sim.dat$W==1),tr,replace=F)
#   co.keep <- sample(which(sim.dat$W==0),co,replace=F)
#   sim.dat <- sim.dat[c(tr.keep,co.keep),]
#   W <- sim.dat[,"W"]
#   X <- as.matrix(sim.dat[,names(sim.dat)[-1]])
  
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
  est.store[i,"RAW"] <- (b >= ttest_ci[1] & b<= ttest_ci[2])
  
  
  # Model matching predictions
  mmlist <- list()
  # Model 1: simple linear regression on all covariates
  model.MM1 <- lm(Y~X + Z1 + Z2)
  mmlist[[1]] <- model.MM1$fitted
  # Model 2: linear regression on all covariates with higher order polynomial terms
  model.MM2 <- lm(Y~X+I(X^2)+I(X^3)+Z1+I(Z1^2)+I(Z1^3)+Z2+I(Z2^2)+I(Z2^3))    
  mmlist[[2]] <- model.MM2$fitted
  # Model 3: bagged regression trees, default 25 bags
  model.MM3 <- bagging(Y~., data = data.frame(Y,X,Z1,Z2), coob=TRUE)
  mmlist[[3]] <- predict(model.MM3, data.frame(X,Z1,Z2))
  # Model 4: random forest
  model.MM4 <- randomForest(Y~., data = data.frame(Y,X,Z1,Z2))
  mmlist[[4]] <- predict(model.MM4, X)
  
  
  # model-based matching
  for(p in 1:4){
    pred <- mmlist[[p]]
    pairs <- Matches(X, pred)
    pairs.df <- data.frame(do.call(rbind,pairs))
    w0 <- 1/table(pairs.df[pairs.df$tr == 0,1]); Y0 <- Y[unique(pairs.df[pairs.df$tr==0,1])]
      mm_res <- permu_CI_mean(groups = pairs, prediction = pred, response = Y, treatment = X, side = "both", alpha = alpha, iters=1000, verbose = TRUE)
      est.store[[k]][i,paste("MM", p, sep="")] <- mm_res[[1]] #mm_res[[1]]/tr # edit: I changed MM function to divide by number treated
      est.store[[k]][i,paste("MM", p, "ATT", sep="")] <- mean(Y[W==1]) - sum(Y0*w0)/sum(w0)
    }}
  
  # entropy balancing
  for(k in 1:3){
    Y <- Y.mat[,k]
    out.eb <- try(ebalance(
      Treatment=W,
      X=X,
      print.level = -1
    ),silent=FALSE)
    if(class(out.eb)=="try-error"){
      cat("EB did not converge in sims",i,"\n")
    } else {
      d <- out.eb$w
      est.store[[k]][i,"EB"]  <- mean(Y[W==1]) - (sum(Y[W==0]*d)/sum(d))
    }}
  
  ## GenMatching
  for(k in 1:3){
    Y <- Y.mat[,k]
    suppressWarnings(gout <- GenMatch(Tr=W, X=X, BalanceMatrix=X, estimand="ATT", M=1, weights=NULL,
                                      print.level=0))
    out  <- Match(Tr=W,X=X,Weight.matrix=gout,estimand="ATT",M=1,BiasAdjust=FALSE)
    est.store[[k]][i,"GM"] <- (sum(Y[out$index.treated]*out$weights)/sum(out$weights)) - (sum(Y[out$index.control]*out$weights)/sum(out$weights))
  }
  
  # Propensity scoring methods
  # Matching
  for(k in 1:3){
    Y <- Y.mat[,k]
    for(p in 1:3){
      if(sum(pslist[[p]]$fitted<(0+sqrt(.Machine$double.eps)))> 0 || sum(pslist[[p]]$fitted>(1-sqrt(.Machine$double.eps)))>0){ # error reporting if perfect sep
        cat("Perfect Seperation in glm occured in sim",i,"\n")
      } else {
        # ps matching
        PS <- pslist[[p]]$linear
        out <- Match(Tr=W,X=PS,estimand="ATT",M=1,BiasAdjust = FALSE)
        est.store[[k]][i,paste("PS",p,sep="")] <- (sum(Y[out$index.treated]*out$weights)/sum(out$weights)) - (sum(Y[out$index.control]*out$weights)/sum(out$weights))
      }
      
      # Weighting
      PS.pr <- pslist[[p]]$fitted
      d <- PS.pr/(1-PS.pr)
      est.store[[k]][i,paste("PSW",p,sep="")] <- mean(Y[W==1])- (sum(Y[W==0]*d[W==0])/sum(d[W==0]))
    }
    # store housekeeping vars
    para.store[[k]][i,"FracCtoT"] <- sum(W==0)/sum(W==1)
    para.store[[k]][i,"corY.PS"] <- cor(Y,PS.true)
    para.store[[k]][i,"corW.PS"] <- cor(W,PS.true)
    para.store[[k]][i,"corW.MM1"] <- cor(W, mmlist[[k]][[1]])
    para.store[[k]][i,"corW.MM2"] <- cor(W, mmlist[[k]][[2]])
    para.store[[k]][i,"corW.MM3"] <- cor(W, mmlist[[k]][[3]])
    para.store[[k]][i,"corPS.PS1"] <- cor(pslist[[1]]$fitted,PS.true)
    para.store[[k]][i,"corPS.PS2"] <- cor(pslist[[2]]$fitted,PS.true)
    para.store[[k]][i,"corPS.PS3"] <- cor(pslist[[3]]$fitted,PS.true)
    para.store[[k]][i,"N"]  <- nrow(sim.dat)
    para.store[[k]][i,"tr"] <- sum(W==1)
    para.store[[k]][i,"co"] <- sum(W==0)
  }
  
  est.summ <- lapply(est.store, function(x)round(apply(x, 2, mean, na.rm=T), 5))
  para.summ <- lapply(para.store, function(x)round(apply(x, 2, mean, na.rm=T), 5))
  mse <- lapply(est.store, function(x) apply(x, 2, function(y) mean((y-0)^2, na.rm=T)))
  bias <- lapply(est.store, function(x) apply(x, 2, function(y) mean(y-0, na.rm=T)))
}

save(est.store, para.store, est.summ, para.summ, mse, bias, file = "simulation1.Rdata")


