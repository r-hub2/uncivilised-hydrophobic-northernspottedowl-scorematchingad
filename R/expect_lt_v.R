# Custom expecations for vector comparisons in testthat
# @examples
# object = c(1, 3, 5)
# expected = 3
# expect_lt_v(c(1,3,5), 3)

expect_lt_v <- function(object, expected){
  expect_op_v(object, expected, operation = function(x, y) x < y, comparelang = "NOT less than")
}

expect_lte_v <- function(object, expected){
  expect_op_v(object, expected, operation = function(x, y) x <= y, comparelang = "larger than")
}

expect_absdiff_lte_v <- function(object1, object2, expected){
  expect_lte_v(abs(object1 - object2), expected)
}

expect_op_v <- function(object, expected, operation = function(x, y) x < y, comparelang = "larger than"){
  act <- testthat::quasi_label(rlang::enquo(object), arg = "object")

  # check
  stopifnot(length(expected) > 0) #important because passing vectors
  stopifnot(length(object) > 0) #important because passing vectors
  if (length(expected) != length(object)){
    if (length(object) %% length(expected) == 0){
      expected <- rep(expected, length(object) / length(expected))
    } else {
      testthat::fail(sprintf("Length mismatch: %s has length %i but 'expected' has length %i", act$lab, length(object), length(expected)))
    }
  }
  comparison <- operation(object, expected)
  ok = all(comparison)
  if (ok){
    testthat::succeed()
    return(invisible(act$val))
  }
  failid <- which(!comparison)
  faildesc <- paste0(sprintf("%s was %s expected at index c(%s).", act$lab, comparelang, paste(failid, collapse = ", ")))
  faildetail <- paste(failid, ": ", object[failid], "!<",expected[failid])
  testthat::fail(paste(c(faildesc, faildetail), collapse = "\n"))
}
