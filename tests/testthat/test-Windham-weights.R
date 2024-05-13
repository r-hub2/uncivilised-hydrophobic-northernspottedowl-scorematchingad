test_that("Windham_weights() fails properly", {
  Y <- matrix(rnorm(30), ncol = 3)
  expect_error(Windham_weights(ldenfun = NULL, Y = Y, theta = rep(0, 10), cW = rep(0.1, 10)))
  
  expect_error({Windham_weights(ldenfun = function(Y, theta){c(Inf, rep(1, nrow(Y)-1))}, Y = Y, theta = rep(0, 10), cW = rep(0.1, 10))}, "[Nn]on-finite")
  
  expect_error({Windham_weights(ldenfun = function(Y, theta){c(NA, rep(1, nrow(Y)-1))}, Y = Y, theta = rep(0, 10), cW = rep(0.1, 10))}, "[Nn]on-finite")
})

