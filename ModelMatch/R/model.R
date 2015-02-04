# library(rpart)
# library(ipred)
# library(randomForest)
# library(e1071)
# library(gam)
# library(gbm)
# library(nnet)
### NEED TO UPDATE - don't use!
choose_model_old <- function(dat){
  ### dat = dataframe of inputs, including the outcome and all predictors. e0 is the outcome
  mod_lm <- lm(e0~., data = dat); mod_lm$fitted
  mod_poly <- lm(dat$e0~poly(as.matrix(dat[,-1], degree = 2, raw = TRUE))); predict(mod_poly)
  mod_tree <- rpart(e0~., data=dat); predict(mod_tree)
  mod_bag <- bagging(e0~., data=dat); predict(mod_bag)
  mod_rf <- randomForest(e0~., data=dat); predict(mod_rf)
  mod_svm <- svm(e0~., data = dat)
  mod_gam <- gam(e0~., data = dat)
  mod_boost <- gbm(e0~., data = dat, bag.fraction = 1, train.fraction = 1)
  mod_nnet <- nnet(e0/6 ~ ., data = dat, size = 2)
  models <- list(mod_lm, mod_poly, mod_tree, mod_bag, mod_rf, mod_svm, mod_gam, mod_nnet, mod_boost)
  mse <- sapply(models[1:length(models)-1], function(x) mean((predict(x) - dat$e0)^2))
  mse <- c(mse, mean((predict(mod_nnet)*6 - dat$e0)^2))
  mse <- c(mse, mean((predict(mod_boost, dat, n.trees = 100) - dat$e0)^2))
  best <- which.min(mse)[1]
  return(models[[best]])
}


