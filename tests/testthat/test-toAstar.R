test_that("fromAstar() has same result for compositional vectors", {
  set.seed(3645)
  p <- 3
  Astar <- rWishart(1, 6, diag(3))[,,1]
  ALbL <- ppi_fromAstar(Astar)

  # two orthogonal vectors from rep(1, p)/p using rows of the Helmert submatrix
  u <- c(1, -1, 0)/5 + rep(1, p)/p
  uL <- u[-p]
  expect_equal(t(u) %*% Astar %*% u,
               t(uL) %*% ALbL$AL %*% uL + t(ALbL$bL) %*% uL + ALbL$const)


  u <- c(1, 1, -2)/10 + rep(1, p)/p
  uL <- u[-p]
  expect_equal(t(u) %*% Astar %*% u,
               t(uL) %*% ALbL$AL %*% uL + t(ALbL$bL) %*% uL + ALbL$const)

})

test_that("toAstar() reverses fromAstar()", {
  set.seed(3645)
  Astar <- rWishart(1, 6, diag(3))[,,1]
  ALbL <- ppi_fromAstar(Astar)
  newAstar <- ppi_toAstar(ALbL$AL, ALbL$bL)
  expect_equal(newAstar + sum(Astar)/(3*3), Astar)
})

