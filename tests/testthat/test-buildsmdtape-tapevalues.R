test_that("Dirichlet with smd values and derivatives independent of tape", {
  u1 = matrix(c(0.001, 0.011, 1 - 0.01 - 0.011), nrow = 1)
  u2 = matrix(c(1,1,1), nrow = 1)
  theta1 = rep(-0.5, 3)
  theta2 = rep(0, 3)
  fixedtheta = rep(FALSE, 3)
  acut = 0.1
  ueval <- matrix(c(0.4, 0.011, 1 - 0.4 - 0.011), nrow = 1)
  thetaeval <- c(-0.1, -0.5, 2)

  compare <- function(uval, thetaeval, tapes1, tapes2){
  expect_equal(pForward0(tapes1$lltape$ptr, ueval, thetaeval), pForward0(tapes2$lltape$ptr, ueval, thetaeval))
  expect_equal(pJacobian(tapes1$lltape$ptr, ueval, thetaeval), pJacobian(tapes2$lltape$ptr, ueval, thetaeval))
  expect_equal(pHessian(tapes1$lltape$ptr, ueval, thetaeval), pHessian(tapes2$lltape$ptr, ueval, thetaeval))

  expect_equal(pForward0(tapes1$smdtape$ptr, thetaeval, ueval), pForward0(tapes2$smdtape$ptr, thetaeval, ueval))
  expect_equal(pJacobian(tapes1$smdtape$ptr, thetaeval, ueval), pJacobian(tapes2$smdtape$ptr, thetaeval, ueval))
  expect_equal(pHessian(tapes1$smdtape$ptr, thetaeval, ueval), pHessian(tapes2$smdtape$ptr, thetaeval, ueval))
  return(NULL)
  }

  #Sphere with minsq
  tapes1 <- buildsmdtape("sim", "sqrt", "sph",
               "dirichlet", u1, rep(NA, 3), thetatape_creator = function(n){theta1},
               bdryw = "minsq", acut = acut)
  tapes2 <- buildsmdtape("sim", "sqrt", "sph",
               "dirichlet", u2, rep(NA, 3), thetatape_creator = function(n){theta2},
               bdryw = "minsq", acut = acut)

  compare(ueval, thetaeval, tapes1, tapes2)

  #Sphere, prodsq
  tapes1 <- buildsmdtape("sim", "sqrt", "sph",
               "dirichlet", u1, rep(NA, 3), thetatape_creator = function(n){theta1},
               bdryw = "prodsq", acut = acut)
  tapes2 <- buildsmdtape("sim", "sqrt", "sph",
               "dirichlet", u2, rep(NA, 3), thetatape_creator = function(n){theta2},
               bdryw = "prodsq", acut = acut)
  compare(ueval, thetaeval, tapes1, tapes2)

  #Simplex
  tapes1 <- buildsmdtape("sim", "identity", "sim",
               "dirichlet", u1, rep(NA, 3), thetatape_creator = function(n){theta1},
               bdryw = "minsq", acut = acut)
  tapes2 <- buildsmdtape("sim", "identity", "sim",
               "dirichlet", u2, rep(NA, 3), thetatape_creator = function(n){theta2},
               bdryw = "minsq", acut = acut)
  compare(ueval, thetaeval, tapes1, tapes2)

  #Simplex, prodsq
  tapes1 <- buildsmdtape("sim", "identity", "sim",
               "dirichlet", u1, rep(NA, 3), thetatape_creator = function(n){theta1},
               bdryw = "prodsq", acut = acut)
  tapes2 <- buildsmdtape("sim", "identity", "sim",
               "dirichlet", u2, rep(NA, 3), thetatape_creator = function(n){theta2},
               bdryw = "prodsq", acut = acut)
  compare(ueval, thetaeval, tapes1, tapes2)
})

