test_that("Fisher-Bingham likelihood runs and matches R code", {
  set.seed(321)
  p <- 3
  A <- matrix(NA, ncol = p, nrow = p)
  A[upper.tri(A)] <- runif(sum(upper.tri(A)))
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  diag(A) <- c(runif(p-1), NA)
  diag(A)[2] <- 0
  A[p,p] <- -sum(diag(A)[1:(p-1)]) #to satisfy the trace = 0 constraint
  m <- runif(p)
  m <- m / sqrt(sum(m^2))
  k <- 2
  theta <- FB_mats2theta(k, m, A)
  expect_equal(FB_theta2mats(theta),
    list(k = k,
         m = m,
         km = k*m,
         A = A))


  set.seed(123)
  sample <- matrix(runif(2 * p, -10, 10), nrow = 2)
  sample <- sample / sqrt(rowSums(sample^2))
  stopifnot(all(abs(sqrt(rowSums(sample^2)) - 1) < 1E-5))

  pman <- manifoldtransform("sph", "identity", "sph")
  lltape <- tapell("FB", sample[1,], rep(NA, length(theta)), pman$tran, function(n){rep(1, n)})$ptr

  expect_equal(pForward0(lltape, sample[1, ], theta), log(qdFB(sample[1, ], k, m, A)),
               ignore_attr = TRUE)## very important to check a tape
  #deriv wrt u
  u <- sample[1, ]
  expect_equal(pJacobian(lltape, sample[1, ], theta),
               attr(numericDeriv(quote(log(qdFB(u, k, m, A))), "u"), "gradient"),
               ignore_attr = TRUE, tolerance = 1E-5)
  #deriv wrt theta
  lltape_theta <- swapDynamic(lltape, theta, sample[1, ])
  qdFB_theta <- function(u, theta){
    mats <- FB_theta2mats(theta)
    qdFB(u, mats$k, mats$m, mats$A)
  }
  expect_equal(pJacobian(lltape_theta, theta, u),
    attr(numericDeriv(quote(log(qdFB_theta(u, theta))), "theta"), "gradient"),
    ignore_attr = TRUE, tolerance = 1E-5)
})

test_that("FB() fits for p = 3", {
  skip_on_cran() #test with various fixed elements is sufficient
  p <- 3
  set.seed(111)
  theta <- runif(p-1 + (p - 1) * p/2 + p, -10, 10)
  thetamats <- FB_theta2mats(theta)

  #simulate
  set.seed(123451)
  sample <- simdd::rFisherBingham(10000, mu = thetamats$k * thetamats$m, Aplus = thetamats$A)

  #Fit
  estobj <- FB(sample)
  expect_absdiff_lte_v(estobj$est$paramvec, theta, 3 * estobj$SE$paramvec)
  expect_lt_v(estobj$SE$paramvec, rep(5))
})

test_that("FB() fits with various fixed elements", {
  p <- 3
  set.seed(111)
  theta <- runif(p-1 + (p - 1) * p/2 + p, -10, 10)
  thetamats <- FB_theta2mats(theta)

  #simulate
  set.seed(12345)
  sample <- simdd::rFisherBingham(10000, mu = thetamats$k * thetamats$m, Aplus = thetamats$A)

  #a fixed A element
  inA <- matrix(NA, nrow = p, ncol = p)
  inA[p, 1] <- inA[1, p] <- thetamats$A[1, p]
  estobj <- FB(sample, A = inA)
  expect_lte_v(abs(estobj$est$A - thetamats$A)[-(p*p)], 3 * estobj$SE$A[-(p*p)])
  expect_equal(estobj$est$A[!is.na(inA)], thetamats$A[!is.na(inA)])
  expect_equal(estobj$SE$A[!is.na(inA)], 0 * thetamats$A[!is.na(inA)])

  #fixed final diagonal element should error
  inA <- matrix(NA, nrow = p, ncol = p)
  inA[p, p] <- thetamats$A[p, p]
  expect_error(FB(sample, A = inA))

  # a fixed Fisher element
  inkm <- rep(NA, p)
  inkm[p] <- thetamats$m[p] * thetamats$k
  estobj <- FB(sample, km = inkm)
  expect_lte_v(abs(estobj$est$km - thetamats$km), 3 * estobj$SE$km + 1E-10)
  expect_equal(estobj$est$km[!is.na(inkm)], thetamats$km[!is.na(inkm)])
  expect_equal(estobj$SE$km[!is.na(inkm)], 0 * thetamats$km[!is.na(inkm)])
})

test_that("smdbjgrad at true parameters is poor for FB()", {
  skip("Only for research.")
  p <- 3
  set.seed(111)
  theta <- runif(p-1 + (p - 1) * p/2 + p, -10, 10)
  thetamats <- FB_theta2mats(theta)

  #simulate
  set.seed(12345)
  Y <- simdd::rFisherBingham(1E6, mu = thetamats$k * thetamats$m, Aplus = thetamats$A)

  # evaluate gradient
  tapes <- buildsmdtape("sph", "identity", "sph", 
                        ll = "FB",
                        Y[1, ], 
                        usertheta = NA * theta)
  smvals <- smvalues_wsum(tapes$smdtape, Y, theta)
  expect_gt(sqrt(sum((smvals$grad/nrow(sample))^2)), 0.001)
})

test_that("FB() with many fixed elements leads to small smdbjgrad", {
  skip("Fixing many of the elements doesn't improve smdbjgrad")
  p <- 3
  set.seed(111)
  theta <- runif(p-1 + (p - 1) * p/2 + p, -10, 10)
  thetamats <- FB_theta2mats(theta)

  #simulate
  set.seed(12345)
  sample <- simdd::rFisherBingham(10000, mu = thetamats$k * thetamats$m, Aplus = thetamats$A)

  #many fixed elements
  intheta <- theta
  intheta[8] <- NA
  tapes <- buildsmdtape("sph", "identity", "sph", 
                        ll = "FB",
                        sample[1, ], 
                        usertheta = intheta)
  smvals <- smvalues_wsum(tapes$smdtape, sample, theta[is.na(intheta)])
  expect_gt(abs(smvals$grad/nrow(sample)), 0.001)
})

test_that("FB() fits for p = 3, small sample, but with terrible accuracy", {
  skip_on_cran() #test with various fixed elements is sufficient
  p <- 3
  set.seed(111)
  theta <- runif(p-1 + (p - 1) * p/2 + p, -10, 10)
  thetamats <- FB_theta2mats(theta)

  #simulate
  set.seed(123)
  sample <- simdd::rFisherBingham(1000, mu = thetamats$k * thetamats$m, Aplus = thetamats$A)

  #Fit
  estobj <- FB(sample)
  expect_absdiff_lte_v(estobj$est$paramvec, theta, 2 * estobj$SE$paramvec)
  expect_lt_v(estobj$SE$paramvec, rep(20))
})

