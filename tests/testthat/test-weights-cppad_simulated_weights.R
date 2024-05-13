set.seed(1234)
m <- rppi_egmodel(1000, maxden = 4)
#simulate weights
set.seed(134)
vw <- virtualweights(m$sample)
acut = 0.1

test_that("cppad_closed() w = rep(1, nrow(Y)) is near the result as if w omitted", {
  p <- ncol(m$sample)
  tapes <- buildsmdtape("sim","sqrt", "sph", "ppi",
                        ytape = rep(1/p, m$p),
                        usertheta = rep(NA, length(m$theta)),
                        bdryw = "minsq", acut = acut)
  out_constant <- cppad_closed(tapes$smdtape, m$sample, w = rep(1, nrow(m$sample)))
  out_ommit <- cppad_closed(tapes$smdtape, m$sample)

  expect_equal(out_ommit$est, out_constant$est)
  expect_equal(out_ommit$value, out_constant$value)
})

test_that("evaltape_wsum() matches for simulated weights and constant weights", {
  intheta <- ppi_paramvec(m$p)
  tapes <- buildsmdtape("sim","sqrt", "sph", "ppi",
                        m$sample[1, ], intheta,
                        bdryw = "minsq",
                        acut = acut)
  smd_u <- tapeSwap(tapes$smdtape)
  smd_sim <- evaltape_wsum(smd_u, vw$newY, m$theta)
  smd_hardcoded <- evaltape_wsum(smd_u, m$sample, m$theta, w=vw$w)
  expect_equal(smd_sim, smd_hardcoded)

  # compare results to manual calculation
  smd_sim_manual <- sum(vapply(1:nrow(vw$newY), function(i){pForward0(tapes$smdtape$ptr, m$theta, vw$newY[i, ])}, FUN.VALUE = 1.3))
  smd_dir_manual_v <- vapply(1:nrow(m$sample), function(i){pForward0(tapes$smdtape$ptr, m$theta, m$sample[i, ])}, FUN.VALUE = 1.3)
  smd_dir_manual <- sum(smd_dir_manual_v*vw$w)
  expect_equal(smd_sim_manual, smd_sim)
  expect_equal(smd_sim_manual, smd_dir_manual)
})

test_that("evaltape_wsum() matches for simulated weights and constant weights with boundary data", {
  intheta <- ppi_paramvec(m$p)
  tapes <- buildsmdtape("sim","sqrt", "sph", "ppi",
                        m$sample[1, ], intheta,
                        bdryw = "minsq",
                        acut = acut)
  smd_u <- tapeSwap(tapes$smdtape)
  Y <- m$sample
  isbdry <- simplex_isboundary(Y, 1E-2)
  Yapproxcentres <- Y 
  Yapproxcentres[!isbdry, ] <- NA
  Yapproxcentres[isbdry, ] <- simplex_boundaryshift(Y[isbdry, , drop = FALSE], shiftsize = 1E-3)
  origYcentres <- Yapproxcentres
  
  Y <- vw$newY 
  isbdry <- simplex_isboundary(Y, 1E-2)
  Yapproxcentres <- Y 
  Yapproxcentres[!isbdry, ] <- NA
  Yapproxcentres[isbdry, ] <- simplex_boundaryshift(Y[isbdry, , drop = FALSE], shiftsize = 1E-3)
  simYcentres <- Yapproxcentres

  smd_sim <- evaltape_wsum(smd_u, vw$newY, pmat = m$theta, xcentres = simYcentres)
  smd_hardcoded <- evaltape_wsum(smd_u, m$sample, pmat = m$theta, xcentres = origYcentres, w=vw$w)
  expect_equal(smd_sim, smd_hardcoded)
})

test_that("cppad_search() for ppi with minsq matches itself", {
  tapes <- buildsmdtape("sim","sqrt", "sph", "ppi",
               ytape = rep(1/m$p, m$p),
               usertheta = ppi_paramvec(m$p),
               bdryw = "minsq",
               acut = acut)

  suppressWarnings({out_sim <- cppad_search(tapes$smdtape, m$theta *0 + 1, vw$newY, control = list(tol = 1E-12, maxit = 10))})
  suppressWarnings({out_dir <- cppad_search(tapes$smdtape, m$theta *0 + 1, m$sample, control = list(tol = 1E-12, maxit = 10), w = vw$w)})
  expect_equal(out_sim[!(names(out_sim) %in% c("counts", "SE"))], 
     out_dir[!(names(out_sim) %in% c("counts", "SE"))],
     tolerance = 2E-3)

  expect_equal(out_dir[c("est", "value", "sqgradsize")], out_sim[c("est", "value", "sqgradsize")], tolerance = 2E-3)
})


