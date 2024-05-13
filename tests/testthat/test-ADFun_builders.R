test_that("Jacobian, Hessian, GradOffset, Swap and LogJacDet all produce new tapes", {
maninfo <- manifoldtransform("sim", "sqrt", "sph")
ppitape <- tapell(ll = "ppi",
                  ytape = c(0.2, 0.3, 0.5),
                  usertheta = ppi_paramvec(p = 3), 
                  tranobj = maninfo$tran) 
expect_s3_class(tapeJacobian(ppitape), "ADFun")
expect_s3_class(tapeHessian(ppitape), "ADFun")
expect_s3_class(tapeGradOffset(ppitape), "ADFun")
expect_condition(tapeLogJacDet(ppitape), class = "Rcpp::exception", regexp = "equal")
expect_s3_class(tapeLogJacDet(tapeJacobian(ppitape)), "ADFun")
expect_s3_class(tapeSwap(ppitape), "ADFun")
})

