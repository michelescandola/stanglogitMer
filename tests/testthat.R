library(testthat)
library(stanglogitMer)

test_check("stanglogitMer")

test_that("glogit and invglogit", {
  tmp = glogit(0:10,-10,10,0.5,0.2)
  expect_equal(tmp, c(-0.4995837,  1.9737532,  4.2189901,  6.0436778,  7.3978305,  8.3365461,  8.9569287,  9.3540907,  9.6031939,  9.7574313,  9.8521692))
  expect_equal(invglogit(tmp,-10,10,0.5,0.2),0:10)
  expect_warning(invglogit(tmp,-10,3,0.5,0.2))
  expect_error(glogit("1"))
  expect_error(invglogit("1"))
})

test_that("stanglogitMer", {
  set.seed(5)
  tmp = data.frame(y=c(glogit(-10:10,-10,10,0.5,0.2)+rnorm(21,0,0.01),glogit(-10:10,-10,10,0.5,0.2)+rnorm(21,0,0.01),
                       glogit(-10:10,-10,10,0.7,0.3)+rnorm(21,0,0.01),glogit(-10:10,-10,10,0.7,0.3)+rnorm(21,0,0.01),
                       glogit(-10:10,-10,10,0.1,0.2)+rnorm(21,0,0.01),glogit(-10:10,-10,10,0.1,0.2)+rnorm(21,0,0.01),
                       glogit(-10:10,-10,10,0.1,0.3)+rnorm(21,0,0.01),glogit(-10:10,-10,10,0.1,0.3)+rnorm(21,0,0.01)),
                       covariate = c(-10:10,-10:10,-10:10,-10:10,-10:10,-10:10,-10:10,-10:10),
                       grouping = factor(rep(c(1,2),each=42)),
                       condition = factor(rep(c(1,2),each=21))
                     )

  tmp$y = tmp$y + (as.numeric(tmp$grouping)*0.5)*as.numeric(tmp$condition)

  out = stanglogitMer(tmp$y,~condition,~condition,~condition,~grouping,covariate=tmp$covariate,data=tmp,cores = 4, chains = 4,control=list(adapt_delta=0.99))

  expect_equal(out[[1]],tmp)
})
