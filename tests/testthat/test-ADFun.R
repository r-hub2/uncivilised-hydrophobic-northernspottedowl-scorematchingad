test_that("ADFun object is returned by tapell, and its values can be accessed", {
maninfo <- manifoldtransform("sim", "sqrt", "sph")
ppitape <- tapell(ll = "ppi",
                  ytape = c(0.2, 0.3, 0.5),
                  usertheta = ppi_paramvec(p = 3), 
                  tranobj = maninfo$tran)

expect_type(ppitape$xtape, "double") #vectors are neither S3 objects nor S4 objects
expect_type(ppitape$dyntape, "double")
expect_type(ppitape$name, "character")
expect_type(ppitape$ptr, "externalptr")

expect_error(ppitape$xtape <- rep(0, 3))
expect_error(ppitape$dyntape <- rep(0, 3))
expect_error(ppitape$name <- rep(0, 3))
expect_error(ppitape$ptre <- rep(0, 3))

})
