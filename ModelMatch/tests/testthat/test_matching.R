context("Strata")

test_that("Matched pairs", {
  easypairs <- Matches(c(rep(0,5), rep(1,5)), 1:10)
  expect_equal((sapply(easypairs, function(x) x[2,1])), rep(5,5))
})

test_that("Stratified groups", {
  easypairs <- Strata(c(rep(0,5), rep(1,5)), 1:10)
  expect_equal((sapply(easypairs, function(x) x[2,1])), 2*(1:5))
})
