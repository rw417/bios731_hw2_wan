# load libraries
library(stats)

get_simdata <- function(n, beta0, beta1) {
  x <- rnorm(n, 0, 1)
  logit_p <- beta0 + beta1 * x
  p <- 1 / (1 + exp(-logit_p))

  y <- rbinom(n, 1, p)
  
  # Put x and y into a matrix
  simdata <- cbind(y, x)
  return(simdata)
}
