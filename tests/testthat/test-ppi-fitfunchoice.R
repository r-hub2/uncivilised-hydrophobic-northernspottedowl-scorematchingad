set.seed(314)
m <- rppi_egmodel(10)

test_that("Correctly chooses Dirichlet", {
  out <- ppi(m$sample, ppi_paramvec(p=3, AL = 0, bL = 0), trans = "sqrt", method = "hardcoded",
             acut = 0.1,
             bdryw = "minsq")
  expect_ppi_str(out, m$p)
  expect_equal(out$info$method, "dir_sqrt_minimah")

  out <- ppi(m$sample, ppi_paramvec(p=3, AL = 0, bL = 0), trans = "sqrt", method = "hardcoded",
             acut = 0.1,
             bdryw = "prodsq")
  expect_ppi_str(out, m$p)
  expect_equal(out$info$method, "dir_sqrt_prodh")
})

test_that("Correctly chooses estimatorlog_ratio", {
  out <- ppi(m$sample, ppi_paramvec(p=3, bL = 0, betap = tail(m$beta, 1)), trans = "alr", method = "hardcoded")
  expect_ppi_str(out, m$p)
  expect_equal(out$info$method, "ppi_alr_gengamma")
})

test_that("Correctly chooses sphere estimators with fixed beta", {
  out <- ppi(m$sample, ppi_paramvec(bL = 0, beta = m$beta), trans = "sqrt", method = "hardcoded",
             acut = 0.1,
             bdryw = "minsq")
  expect_ppi_str(out, m$p)
  expect_equal(out$info$method, "estimator1_zerob")
  out <- ppi(m$sample, ppi_paramvec(bL = 0, beta = m$beta), trans = "sqrt", method = "hardcoded",
             acut = 0.1,
             bdryw = "prodsq")
  expect_ppi_str(out, m$p)
  expect_equal(out$info$method, "ppi_sqrt_prodh_zerob")

  out <- ppi(m$sample, ppi_paramvec(beta = m$beta), trans = "sqrt", method = "hardcoded",
             acut = 0.1,
             bdryw = "minsq")
  expect_ppi_str(out, m$p)
  expect_equal(out$info$method, "estimator1_incb")
  out <- ppi(m$sample, ppi_paramvec(beta = m$beta), trans = "sqrt", method = "hardcoded",
             acut = 0.1,
             bdryw = "prodsq")
  expect_ppi_str(out, m$p)
  expect_equal(out$info$method, "ppi_sqrt_prodh")
})

test_that("Correctly chooses sphere estimators for unfixed beta", {
  suppressWarnings({out <- ppi(m$sample, ppi_paramvec(p=3, betap = tail(m$beta, 1)), trans = "sqrt", method = "hardcoded",
             acut = 0.1,
             bdryw = "minsq")})
  expect_ppi_str(out, m$p)
  expect_equal(out$info$method, "ppi_sqrt_minimah_betap")
  out <- ppi(m$sample, trans = "sqrt", method = "hardcoded",
             acut = 0.1,
             bdryw = "minsq")
  expect_ppi_str(out, m$p)
  expect_equal(out$info$method, "ppi_sqrt_minimah_full")
})

test_that("Correctly chooses cppad", {
  # full hardcoded with prodsq doesn't exist
  expect_warning(out <- ppi(m$sample, ppi_paramvec(p=3, betap = tail(m$beta, 1)), trans = "sqrt", method = "hardcoded",
             acut = 0.1,
             bdryw = "prodsq"))
  expect_ppi_str(out, m$p)
  expect_equal(out$info$method, "closed")

  expect_warning({out <- ppi(m$sample, trans = "sqrt", method = "hardcoded",
             acut = 0.1,
             control = list(tol = 1E-10),
             bdryw = "prodsq")}, "hard-coded")
  expect_ppi_str(out, m$p)

  out <- ppi(m$sample, ppi_paramvec(beta = m$beta), trans = "sqrt", method = "closed",
             acut = 0.1,
             bdryw = "prodsq")
  expect_ppi_str(out, m$p)
  expect_equal(out$info$method, "closed")
})
