context("Matching functions")

test_that("Matched pairs", {
  easypairs <- Matches(c(rep(0,5), rep(1,5)), 1:10)
  expect_equal((sapply(easypairs, function(x) x[2,1])), rep(5,5))

  treatment <- c(1, 1, 0, 0, 0)
  X <- c(1, 2, 1, 3, 7)
  matches <- cbind(c(1,2), c(3,3))
  set.seed(2)
  res <- Matches(treatment, X)
  index <- sapply(res, function(x) x[,"index"])
  expect_equal(matches, t(index))
})

test_that("Stratified groups", {
  easypairs <- Strata(c(rep(0,5), rep(1,5)), 1:10)
  expect_equal((sapply(easypairs, function(x) x[2,1])), 2*(1:5))

  treatment <- c(1, 1, 0, 0, 0)
  X <- c(1, 2, 1, 3, 7)
  matches <- cbind(c(1,2), c(3,3))
  set.seed(2)
  res <- Strata(treatment, X, 2)
  index <- sapply(res, function(x) x[,"index"])
  expect_equal(index, list(1:3, 4:5))
})
