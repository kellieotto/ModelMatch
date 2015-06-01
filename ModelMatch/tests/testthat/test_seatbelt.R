

data(Seatbelts)
seat <- as.data.frame(Seatbelts)
mod <- loess(DriversKilled~PetrolPrice, data = seat)
## test:
strata <- Matches(treatment=seat$law, prediction=mod$fitted)
perm_strata <- permute_within_groups(strata)
resid <- mod$residuals


### test 1: null hypothesis shift = 0
res <- permu_test_mean(strata, prediction = mod$fitted, treatment = seat$law, response = seat$DriversKilled)
hist(res$perm_distribution/length(strata)); abline(v=res$diff_means, col = "red")
### test 2: null hypothesis shift = -23
res <- permu_test_mean(strata, prediction = mod$fitted, treatment = seat$law, response = seat$DriversKilled, shift = -23)
hist(res$perm_distribution/length(strata)); abline(v=res$diff_means, col = "red")


d <- permu_CI_mean(groups=strata, prediction=mod$fitted, treatment = seat$law, response=seat$DriversKilled, side = "both", iters=1000, verbose = TRUE)
abline(v = d[1], lty = 2); abline(v = d[2], lty = 2)

context("Strata")

easypairs <- Matches(c(rep(0,5), rep(1,5)),1:10)
test_that("Matched pairs", {
  expect_equal((sapply(easypairs, function(x)x[2])), rep(5,5))
})

