# Use optim() to find the estimates
optim_function <- function(sim_data, init = c(0,0)){
  # Grab data
  y = sim_data[,1]
  x = cbind(rep(1, nrow(sim_data)), sim_data[,2])
  
  log_likelihood <- function(beta, x, y) {
    p <- 1 / (1 + exp(-x %*% beta)) 
    return(-sum(y * log(p) + (1 - y) * log(1 - p)))
  }
  
  # Initial guess for beta
  beta_init <- init

  # Optimize using BFGS method
  result <- optim(par = beta_init, fn = log_likelihood, x = x, y = y, method = "BFGS")
  
  return(list(
    solution = result$par,
    iter = result$counts["gradient"]
    )
  )
}
