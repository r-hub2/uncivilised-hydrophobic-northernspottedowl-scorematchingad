#interior point
set.seed(124)
u1 <- matrix(runif(5, 0, 1), nrow = 1)
u1 <- u1 / sum(u1)
#boundary point
set.seed(12465)
u2 <- matrix(c(0, runif(4, 0, 1)), nrow = 1)
u2 <- u2 / sum(u2)
utabl <- rbind(u1, u2)

test_that("simplex_isboundary() works on a single measurement", {
  #interior point
  out <- simplex_isboundary(u1, bdrythreshold = 1E-15)
  expect_equal(out, FALSE, ignore_attr = TRUE)

  #boundary point
  out <- simplex_isboundary(u2, bdrythreshold = 1E-15)
  expect_equal(out, TRUE, ignore_attr = TRUE)

  #a table
  out <- simplex_isboundary(utabl, bdrythreshold = 1E-15)
  expect_equal(out, c(FALSE, TRUE), ignore_attr = TRUE)
})

test_that("simplex_boundaryshift gives results equal to shiftsize", {
  centres <- simplex_boundaryshift(Y = utabl, shiftsize = 1E-5)
  expect_equal(sqrt(rowSums((centres - utabl)^2)), rep(1E-5, nrow(utabl)))
})

test_that("simplex_boundarysplit() combines simplex_boundaryshift and simplex_isboundary correctly", {
  splitl <- simplex_boundarysplit(utabl, bdrythreshold = 1E-10, shiftsize = 1E-10)
  expect_equal(splitl$interior, u1)
  expect_equal(splitl$uboundary, u2)
  expect_equal(sqrt(rowSums((splitl$uboundary - splitl$boundaryapprox)^2)), rep(1E-10, nrow(splitl$uboundary)))
})
