# test guttman.R functions

test_that("Guttman's mu2 works",
          {expect_equal(print(mu2(x = 1:12, y = c(1,9,10,2,3,4,5,8,7,6,11,12)), 12),
             0.657060518732) } )

test_that("Guttman's mu2df returns a matrix",
          {expect_equal(class(mu2df(df = data.frame(a = 1:5, b = c(2, 3, 1, 5, 4), c = 5:1))), c("matrix", "array"))})

