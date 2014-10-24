

data(Seatbelts)
seat <- as.data.frame(Seatbelts)
mod <- loess(DriversKilled~PetrolPrice, data = seat)
## test:
strata <- Matches(treatment=seat$law, prediction=mod$fitted)
perm_strata <- permute_within_groups(strata)
resid <- mod$residuals
permu_test_mean(strata, prediction = mod$fitted, response = seat$DriversKilled)


context("Strata")

easypairs <- Matches(c(rep(0,5), rep(1,5)),1:10)
test_that("Matched pairs", {
  expect_equal((sapply(easypairs, function(x)x[2])), rep(5,5))
})

