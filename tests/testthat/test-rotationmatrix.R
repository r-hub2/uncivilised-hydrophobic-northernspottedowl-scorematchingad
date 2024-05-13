test_that("Rotation matrix generation on high dimensions", {
  u <- c(1, 2, 3, 4, 5, 6)
  b <- c(3,4,5,1,2,6)
  Rmat <- rotationmatrix(u, b)
  expect_equal(u, Rmat %*% b, ignore_attr = TRUE)
})


test_that("Handedness of coord system is preserved", {
  vector.cross <- function(a, b) { #code from https://stackoverflow.com/questions/15162741/what-is-rs-crossproduct-function
    if(length(a)!=3 || length(b)!=3){
      stop("This cross product function is only defined for 3D vectors.");
    }
    i1 <- c(2,3,1)
    i2 <- c(3,1,2)
    return (a[i1]*b[i2] - a[i2]*b[i1])
  }

  e1 <- c(1, 0, 0)
  e2 <- c(0, 1, 0)
  e3 <- c(0, 0, 1)

  stopifnot(isTRUE(all.equal(vector.cross(e1, e2), e3)))

  #the vector cross product gives direction in the right-hand rule
  #so a transformation R of e1, e2, e3 that preserves sizes, angles and handedness should be s.t.
  # Re1 x Re2 = Re3
  u <- c(1, 2, 3)
  Rmat <- rotationmatrix(e1, u)

  expect_equal(vector.cross(Rmat %*% e1, Rmat %*% e2), Rmat %*% e3,
               ignore_attr = TRUE)

  # also rotating e2 to e1 doesn't just switch the order of the dimensions
  Rmat <- rotationmatrix(e1, e2)
  expect_equal(vector.cross(e3, Rmat %*% e2), e2,
               ignore_attr = TRUE)
})

test_that("rotation is the simplest for e1, e2, e3", {
  # rotating v1 to v2 can be acheived in two ways.
  # Suppose the angle between v1 and v2 is theta.
  # A rotation of theta in the axis perpendicular to v1 and v2
  # OR
  # rotation of pi around v2 axis, then rotation of 2*pi - theta in the axis perpendicular to v1 and v2
  # the former seems SIMPLER and better
  # for z to x that means y remains unchanged, x goes to -z
  e1 <- c(1, 0, 0)
  e2 <- c(0, 1, 0)
  e3 <- c(0, 0, 1)
  Rmat <- rotationmatrix(e1, e3)
  expect_equal(Rmat %*% e3, e1, ignore_attr = TRUE)
  expect_equal(Rmat %*% e2, e2, ignore_attr = TRUE)
  expect_equal(Rmat %*% e1, -e3, ignore_attr = TRUE)
})

test_that("Rotation matrix rotates opposite direction to one", {
  u <- c(0, 2, 0, 0, 0, 0)
  uend <- c(0, -1, 0, 0, 0, 0)
  Rmat <- rotationmatrix(uend, u)

  expect_equal(uend * 2, Rmat %*% u, ignore_attr = TRUE) 
})
