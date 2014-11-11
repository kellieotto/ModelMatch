dat_ctrl <- read.table("../Data/cps_controls.txt", header=F)
dat_tr   <- read.table("../Data/nswre74_treated.txt", header=F)

dat <- rbind(dat_ctrl, dat_tr)
colnames(dat) <-  c('Treated', 'Age', 'Education', 'Black', 'Hispanic', 'Married',
                              'Nodegree', 'RE74', 'RE75', 'RE78')
mod <- lm(RE78~., data = dat[,-1])


## test:
strata <- Matches(treatment=dat$Treated, prediction=mod$fitted)
perm_strata <- permute_within_groups(strata)
resid <- mod$residuals
res <- permu_test_mean(strata, prediction = mod$fitted, response = dat$RE78)
hist(res$perm_distribution); abline(v=res$diff_means)