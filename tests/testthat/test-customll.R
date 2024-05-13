test_that("customll_test() returns TRUE or FALSE", {
  suppressWarnings({out <- customll_test()})
  expect_true(out %in% c(TRUE, FALSE))
})

# warning: for interactive testing load_all() with install doesn't include scorematchingad.h properly
skip_if_not(suppressWarnings(customll_test()), "Taping customll not available on this system")
test_that("customll can compile a ll function that is evaluated by evalll and taping later", {
  dirll <- customll("
  a1type dirichlet(const veca1 &u, const veca1 &beta) {
          size_t d  = u.size();
          a1type y(0.);  // initialize summation
          for(size_t i = 0; i < d; i++)
          {   y   += beta[i] * log(u[i]);
          }
          return y;
  }")
  expect_equal(evalll(dirll, rep(1/3, 3), rep(-0.5, 3)), 3 * (-0.5 * log(1/3)))
  
  psimplex <- manifoldtransform("sim", "identity", "sim")
  lltape <- tapell(dirll, c(0.1, 0.4, 0.5), rep(NA, 3), psimplex$tran, verbose = FALSE)

  expect_equal(pForward0(lltape$ptr, rep(1/3, 3), rep(-0.5, 3)), 3 * (-0.5 * log(1/3)))
  expect_equal(pJacobian(lltape$ptr, rep(1/3, 3), rep(-0.5, 3)), rep(-0.5 * 3, 3))
})

test_that("customll errors correctly with wrong signature", {
  expect_error({dirll <- customll("
a1type dirichlet(const vecd &u, const veca1 &beta) {
        size_t d  = u.size();
        a1type y(0.);  // initialize summation
        for(size_t i = 0; i < d; i++)
        {   y   += beta[i] * log(u[i]);
        }
        return y;
}")})
})

test_that("a customll gets all the way to the correct score matching estimate", {
  ll <- customll("
  a1type customppi(const veca1 &u, const veca1 &theta) {
          a1type y;
          y = ll::ll_ppi(u, theta);
          return y;
  }")
  set.seed(13411)
  mod <- rppi_egmodel(100)
  tran <- manifoldtransform("sim", "sqrt", "sph")$tran
  Y <- mod$sample

  tapes_custom <- buildsmdtape(
     start = "sim",
     tran = "alr",
     end = "Euc",
     ll = ll,
     ytape = c(0.2, 0.3, 0.5),
     usertheta = ppi_paramvec(p = 3, betap = tail(mod$beta, 1)),
     bdryw = "ones",
     verbose = FALSE)
  tapes <- buildsmdtape(
     start = "sim",
     tran = "alr",
     end = "Euc",
     ll = "ppi",
     ytape = c(0.2, 0.3, 0.5),
     usertheta = ppi_paramvec(p = 3, betap = tail(mod$beta, 1)),
     bdryw = "ones",
     verbose = FALSE)

  est_custom <- cppad_closed(tapes_custom$smdtape, Y) 
  est <- cppad_closed(tapes$smdtape, Y) 
  expect_equal(est_custom$est, est$est)
  expect_equal(est_custom$covar, est$covar)
})
