

data(Seatbelts)
seat <- as.data.frame(Seatbelts)
mod <- loess(DriversKilled~PetrolPrice, data = seat)
## test:
strata <- Matches(treatment=seat$law, prediction=mod$fitted)
perm_strata <- permute_within_groups(strata)
resid <- mod$residuals
res <- permu_test_mean(strata, prediction = mod$fitted, response = seat$DriversKilled)
hist(res$perm_distribution/length(strata)); abline(v=res$diff_means, col = "red")
d <- permu_CI_mean(groups=strata, prediction=mod$fitted, response=seat$DriversKilled, side = "lower",iters=100, precision = 25)
abline(v = d[1], lty = 2); abline(v = d[2], lty = 2)

context("Strata")

easypairs <- Matches(c(rep(0,5), rep(1,5)),1:10)
test_that("Matched pairs", {
  expect_equal((sapply(easypairs, function(x)x[2])), rep(5,5))
})

