test_that("von-Mises Fisher likelihood runs and fits", {
  set.seed(123)
  theta <-  3 * c(1, -1) / sqrt(2)
  if (!requireNamespace("movMF")){skip("Need movMF package")}
  sample <- movMF::rmovMF(1000, theta)
  thetatape <- c(1, -1)

  p <- length(theta)
  intheta <- c(NA, rep(0, p - 1))
  tapes <- buildsmdtape("sph", "identity", "sph", "vMF",
                        rep(1, p)/sqrt(p), rep(NA, p),
                        bdryw = "ones",
                        verbose = FALSE)
  expect_equal(pForward0(tapes$lltape$ptr, sample[1, ], theta), sum(sample[1, ]  * theta)) ## very important to check a tape
  out <- cppad_closed(tapes$smdtape, Y = sample)
  expect_absdiff_lte_v(out$est, theta, 3 * out$SE)
})

test_that("vMF_Mardia() function works for data centred off the north pole", {
  set.seed(123)
  k <- 3
  m <- c(1, -1) / sqrt(2)
  km <-  k * m
  if (!requireNamespace("movMF")){skip("Need movMF package")}
  sample <- movMF::rmovMF(1000, km)
  out <- vMF(Y = sample, method = "Mardia")
  expect_equal(out$est$m , m, tolerance = 1E-1) #moment estimate part
  expect_lt_v(abs(out$est$k - k), 3 * out$SE$k)
})

test_that("vMF_Full() function works", {
  set.seed(123)
  k <- 3
  m <- c(1, -1) / sqrt(2)
  km <-  k * m
  if (!requireNamespace("movMF")){skip("Need movMF package")}
  sample <- movMF::rmovMF(100, km) #faithful to seed
  out <- vMF(Y = sample, method = "smfull")
  expect_lt_v(abs(out$est$paramvec - km), 3 * out$SE$paramvec)

  #with a fixed component
  inkm <- km
  inkm[2] <- NA
  out <- vMF(Y = sample, paramvec = inkm, method = "smfull")
  expect_lte_v(abs(out$est$paramvec - km), 3 * out$SE$paramvec)
})

test_that("vMF() fitting works on dimension 5", {
  set.seed(123)
  p <- 5
  k <- 3
  m <- runif(p, min = -10, 10)
  m <- m / sqrt(sum(m^2))
  km <-  k * m
  set.seed(12412)
  if (!requireNamespace("movMF")){skip("Need movMF package")}
  sample <- movMF::rmovMF(1000, m * k)
  #full method
  out <- vMF(Y = sample, method = "smfull")
  expect_lt_v(abs(out$est$paramvec - km), 3 * out$SE$paramvec)
  #full with a fixed components
  inkm <- km
  inkm[2] <- NA
  inkm[p] <- NA
  out <- vMF(Y = sample, paramvec = inkm, method = "smfull")
  expect_lte_v(abs(out$est$paramvec - km), 3 * out$SE$paramvec)
  expect_equal(out$SE$paramvec[!is.na(inkm)], rep(0, sum(!is.na(inkm))))

  #Mardia method
  out <- vMF(Y = sample, method = "Mardia")
  expect_equal(out$est$m , m, tolerance = 1E-1) #moment estimate part
  expect_lt_v(abs(out$est$k - k), 3 * out$SE$k)
})

test_that("vMF matches for simulated weights, ignoring SE, which shouldn't match", {
  set.seed(123)
  p <- 3
  k <- 3
  m <- runif(p, min = -10, 10)
  m <- m / sqrt(sum(m^2))
  km <-  k * m
  set.seed(1231)
  if (!requireNamespace("movMF")){skip("Need movMF package")}
  Y <- movMF::rmovMF(10, km)
  #simulate weights
  set.seed(1342)
  vw <- virtualweights(Y)

  set.seed(321)
  sim1 <- vMF(vw$newY, method = "Mardia")
  set.seed(321)
  dir1 <-  vMF(Y, method = "Mardia", w = vw$w)
  expect_equal(sim1$est[c("k", "m")], dir1$est[c("k", "m")])

  sim2 <- vMF(vw$newY, method = "smfull")
  dir2 <-  vMF(Y, method = "smfull", w = vw$w)
  expect_equal(sim2$est[c("k", "m")], dir2$est[c("k", "m")])

  sim3_m <- vMF_m(vw$newY)
  dir3_m <- vMF_m(Y, w = vw$w)
  expect_equal(sim3_m, dir3_m)
  expect_error(expect_equal(dir3_m, vMF_m(Y)), "not equal to")
})


test_that("vMF() robust fitting works on dimension 5 with direction outliers", {
  set.seed(123)
  p <- 5
  k <- 3
  m <- rep(1, p) #uniform direction poorness due to outliers is evenly distributed in each element
  m <- m / sqrt(sum(m^2))
  km <-  k * m
  set.seed(121)
  if (!requireNamespace("movMF")){skip("Need movMF package")}
  sample <- movMF::rmovMF(100, km)
  # add outliers
  set.seed(2151)
  outliers <- movMF::rmovMF(10, -km)
  sample_o <- rbind(sample, outliers)
  #full method, not robust
  out1 <- vMF(Y = sample_o, method = "smfull")
  #full method, robust, expect to be closer to true value (due to the outliers)
  out2 <- vMF_robust(Y = sample_o, cW = rep(0.1, 5), method = "smfull")
  expect_true(all(abs(out2$est$paramvec - km) < abs(out1$est$paramvec - km)))

  #check that fixed components remain fixed for full method
  inkm <- km
  inkm[2] <- NA
  inkm[p] <- NA
  out <- vMF_robust(Y = sample_o, cW = 0.1 * is.na(inkm), paramvec = inkm, method = "smfull")
  expect_equal(out$est$paramvec[!is.na(inkm)], km[!is.na(inkm)])

  #Mardia method
  out1 <- vMF(Y = sample_o, method = "Mardia")
  #full method, robust, expect to be closer to true value (due to the outliers)
  out2 <- vMF_robust(Y = sample_o, cW = rep(0.1, length(km)), method = "Mardia")
  expect_true(all(abs(out2$paramvec - km) < abs(out1$est$paramvec - km)))
  #robust estimate of kappa, the outlying points look like low concentration so I'd expect robustness of kappa to do better than non-robust method, and potentially about the same as the full robust method.
  out3 <- vMF_kappa_robust(Y = sample_o, cW = 0.1)
  expect_lt(abs(out3$est$k - k), abs(out1$est$k - k))
})

test_that("robust vMF() with concentration outliers: ok with full robust, better with hybrid, p = 5", {
  set.seed(123)
  p <- 5
  k <- 10 #probably more realistic than 3 - higher concentration means outliers are possible - low concentration means distribution looks uniform
  m <- rep(1, p) #uniform direction poorness due to outliers is evenly distributed in each element
  m <- m / sqrt(sum(m^2))
  km <-  k * m
  set.seed(121)
  if (!requireNamespace("movMF")){skip("Need movMF package")}
  sample <- movMF::rmovMF(100, km)
  # add outliers in concentration only
  set.seed(2151)
  outliers <- movMF::rmovMF(5, 0.1 * km)
  sample_o <- rbind(sample, outliers)

  #Mardia method
  out1 <- vMF(Y = sample_o, method = "Mardia")
  expect_equal(out1$est$m, vMF_m(sample_o)) #because m is estimated as the mean direction
  #full method, robust, expect to be closer to true value overall, but have mixed results
  out2 <- vMF_robust(Y = sample_o, cW = rep(1E-2, ncol(sample_o)), method = "Mardia")
  #expect partially robust Mardia method to be better at k
  out3 <- vMF_kappa_robust(Y = sample_o, cW = 0.01)

  #mixed results for full robust - choosing mean direction
  expect_lt(abs(out2$est$k - k), abs(out1$est$k - k))
  expect_false(all(abs(out2$est$paramvec - km) < abs(out1$est$paramvec - km)))

  # for hybrid, results more solid
  expect_true(all(abs(out3$k * out3$m - km) < abs(out1$est$paramvec - km)))
  expect_lt(abs(out3$est$k - k), abs(out1$est$k - k))
})


test_that("controls of FixedPoint() are passed", {
  set.seed(123)
  p <- 3
  k <- 3
  m <- runif(p, min = -10, 10)
  m <- m / sqrt(sum(m^2))
  km <-  k * m
  set.seed(123)
  if (!requireNamespace("movMF")){skip("Need movMF package")}
  Y <- movMF::rmovMF(10, m * k)

  out_default <- vMF_robust(Y, cW = rep(0.1, ncol(Y)), method = "Mardia") #use this packages defaults this is pretty fussy!

  suppressWarnings(out1 <- vMF_robust(Y, cW = rep(0.1, ncol(Y)), method = "Mardia",
             fpcontrol = list(MaxIter = 2))) #Fixed point iterations of only 2
  expect_equal(out1$info$fpevals, 2)
  expect_error(expect_equal(out_default, out1))
})

test_that("vMF(), vMF_stdY() and vMF_m() differs when weights differ", {
  p <- 5
  k <- 3
  m <- runif(p, min = -10, 10)
  m <- m / sqrt(sum(m^2))
  if (!requireNamespace("movMF")){skip("Need movMF package")}
  Y <- movMF::rmovMF(100, k * m)

  m <- vMF_m(Y)
  set.seed(13423)
  m_w <- vMF_m(Y, w = runif(nrow(Y)))
  expect_gt(sqrt(sum((m_w - m)^2)), 1E-5)


  expect_equal(vMF(Y, method = "Mardia")$est$m, m)
  set.seed(13423)
  expect_equal(vMF(Y, method = "Mardia", w = runif(nrow(Y)))$est$m, m_w)

  expect_equal(colMeans(vMF_stdY(Y))[-1], rep(0, p-1))
  set.seed(13423)
  expect_gt(sqrt(sum(colMeans(vMF_stdY(Y, w = runif(nrow(Y))))[-1]^2)), 1E-5)
})
