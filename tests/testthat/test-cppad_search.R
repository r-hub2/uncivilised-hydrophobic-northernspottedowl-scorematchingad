# test of iterative solver for ppi (compare to cppad_closet()):
## microbiomfit with outliers

test_that("cppad_search gives similar result to cppad_closed", {
  set.seed(354)
  m <- rppi_egmodel(100, maxden = 4)
  tapes <- buildsmdtape("sim","sqrt", "sph", "ppi",
                        ytape = rep(1/m$p, m$p),
                        usertheta = rep(NA, length(m$theta)),
                        bdryw = "minsq", acut = 0.1)
  estsearch <- cppad_search(tapes$smdtape, m$theta *0 + 1, m$sample, control = list(tol = 1E-13))
  estclosed <- cppad_closed(tapes$smdtape, m$sample)
  expect_equal(estsearch$est, estclosed$est, tolerance = 1E-3, ignore_attr = "names")
  expect_equal(estsearch$SE, estclosed$SE, tolerance = 1E-3)
})


