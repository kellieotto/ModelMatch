context("Permutation test helper functions")

test_that("Rcpp and slow versions of within_group_mean match"{
  matches <- list(data.frame("index"=c(1,2), "score"=c(3,4), "treatment"=c(1,0)), data.frame("index"=c(3,4), "score"=c(3,4), "treatment"=c(1,0)))
  Yhat <- rep(1,4)
  Y <- c(2,1,2,1)
  res1 <- within_group_mean(matches, Yhat, Y)
  res2 <- within_group_mean_cpp(matches, Yhat, Y)
  expect_equal(res1, res2)
})

test_that("permute within groups"{
  matches <- list(data.frame("index"=c(1,2), "score"=c(3,4), "treatment"=c(1,0)), data.frame("index"=c(3,4), "score"=c(3,4), "treatment"=c(1,0)))
  Yhat <- rep(1,4)
  Y <- c(2,1,2,1)

  set.seed(1)
  matches1 <- permute_within_groups(matches)
  perm <- sapply(matches1, function(x) x[,"treatment"])
  expect_equal(perm, matrix(c(1, 0, 0, 1), nrow = 2))
  set.seed(1)
  matches2 <- permute_within_groups_cpp(matches)
  perm <- sapply(matches2, function(x) x[,"treatment"])
  expect_equal(perm, matrix(c(0, 1, 0, 1), nrow = 2))
})
