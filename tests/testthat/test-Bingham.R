test_that("taped Bingham log-likelihood gives correct values", {
  p <- 4
  set.seed(345)
  theta <- runif(p-1 + (p - 1) * p/2)

  A <- Bingham_theta2Amat(theta)
  expect_equal(Bingham_Amat2theta(A), theta)

  set.seed(123)
  sample <- simdd::rBingham(2, A)

  pman <- manifoldtransform("sph", "identity", "sph")
  lltape <- tapell(ll = "Bingham", ytape = sample[2, ],
                    usertheta = rep(NA, length(theta)),
                    tranobj = pman$tran)

  u <- t(sample[1, , drop = FALSE])

  expect_equal(pForward0(lltape$ptr, u, theta), t(u) %*% A %*% u,
               ignore_attr = TRUE)

  expect_equal(pJacobian(lltape$ptr, u, theta), 2 * A %*% u,
               ignore_attr = TRUE)

  # test deriv wrt theta
  llBingham <- function(theta){
    A <- Bingham_theta2Amat(theta)
    return(t(u) %*% A %*% u)
  }
  Rgradt <- numericDeriv(quote(llBingham(theta)), c("theta"))
  lltape_t <- tapeSwap(lltape)
  expect_equal(pJacobian(lltape_t$ptr, theta, u), attr(Rgradt, "gradient"),
               tolerance = 1E-5, ignore_attr = TRUE)
})

test_that("Bingham() with smfull optimiser works", {
  p <- 4
  set.seed(345)
  theta <- runif(p-1 + (p - 1) * p/2)
  A <- Bingham_theta2Amat(theta)

  set.seed(123)
  sample <- simdd::rBingham(100, A)
  est <- Bingham(sample, method = "smfull")
  expect_absdiff_lte_v(est$est$paramvec, theta, 3 * est$SE$paramvec)
})

test_that("Bingham() w Mardia optimiser works", {
  p <- 5
  set.seed(345)
  theta <- runif(p-1 + (p - 1) * p/2)
  A <- Bingham_theta2Amat(theta)

  set.seed(123)
  sample <- simdd::rBingham(100, A)
  est <- Bingham(sample, method = "Mardia")
  A_es <- eigen(A)
  expect_absdiff_lte_v(est$est$evals[-p], A_es$values[-p], 3 * est$SE$evals[-p])
})

test_that("Bingham() chooses methods correctly", {
  p <- 4
  set.seed(345)
  theta <- runif(p-1 + (p - 1) * p/2)
  A <- Bingham_theta2Amat(theta)

  set.seed(123)
  sample <- simdd::rBingham(10, A)

  expect_equal(Bingham_full(sample), Bingham(sample, method = "smfull"))
  expect_equal(Bingham_Mardia(sample), Bingham(sample, method = "Mardia"))
  expect_equal(Bingham_Mardia(sample), Bingham(sample, method = "hybrid"))
})

test_that("Bingham() works with highly skewed trace", {
  skip_on_cran() #slow
  p <- 4
  set.seed(111)
  theta <- runif(p-1 + (p - 1) * p/2)
  A <- Bingham_theta2Amat(theta)
  diag(A) <- c(30, 1, 0.1, NA) * diag(A)
  diag(A)[p] <- -sum(diag(A), na.rm = TRUE)
  theta <- Bingham_Amat2theta(A)

  set.seed(123)
  sample <- simdd::rBingham(1000, A)

  est <- Bingham(sample, method = "smfull")
  expect_absdiff_lte_v(est$est$paramvec, theta, 3 * est$SE$paramvec)

  estM <- Bingham(sample, method = "Mardia")
  A_es <- eigen(A)
  expect_absdiff_lte_v(estM$est$evals[-p], A_es$values[-p], 3 * estM$SE$evals[-p])
})

test_that("Bingham() with various fixed elements works", {
  skip_on_cran() #slow
  p <- 4
  set.seed(345)
  A <- rsymmetricmatrix(p, -10, 10)
  A[p,p] <- -sum(diag(A)[-p])
  theta <- Bingham_Amat2theta(A)
  sample <- simdd::rBingham(1000, A)

  #a fixed off diagonal element
  inA <- matrix(NA, nrow = p, ncol = p)
  inA[p, 1] <- inA[1, p] <- A[1, p]
  est <- Bingham(sample, A = inA, method = "smfull")
  expect_lte_v(abs(est$est$A - A)[-(p*p)], 3 * est$SE$A[-(p*p)])
  expect_equal(est$est$A[!is.na(inA)], A[!is.na(inA)])
  expect_equal(est$SE$A[!is.na(inA)], 0 * A[!is.na(inA)])
  expect_error(Bingham(sample, A = inA, method = "Mardia"))
 
  #a fixed diagonal element
  inA <- matrix(NA, nrow = p, ncol = p)
  inA[2, 2] <- A[2, 2]
  est <- Bingham(sample, A = inA, method = "smfull")
  expect_lte_v(abs(est$est$A - A)[-(p*p)], 3 * est$SE$A[-(p*p)])
  expect_equal(est$est$A[!is.na(inA)], A[!is.na(inA)])
  expect_equal(est$SE$A[!is.na(inA)], 0 * A[!is.na(inA)])
  expect_error(Bingham(sample, A = inA, method = "Mardia"))

  #fixed final diagonal element should error
  inA <- matrix(NA, nrow = p, ncol = p)
  inA[p, p] <- A[p, p]
  expect_error(Bingham(sample, A = inA, method = "smfull"))
})

test_that("Bingham methods match for simulated weights, ignoring SE, which shouldn't match", {
  p <- 4
  set.seed(345)
  theta <- runif(p-1 + (p - 1) * p/2)
  A <- Bingham_theta2Amat(theta)

  set.seed(123)
  Y <- simdd::rBingham(100, A)
  
  #simulate weights
  set.seed(1342)
  vw <- virtualweights(Y)

  sim1 <- Bingham(vw$newY, method = "Mardia")
  dir1 <-  Bingham(Y, method = "Mardia", w = vw$w)
  expect_equal(sim1$est, dir1$est)

  sim2 <- Bingham(vw$newY, method = "smfull")
  dir2 <-  Bingham(Y, method = "smfull", w = vw$w)
  expect_equal(sim2$est, dir2$est)
})
