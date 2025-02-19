# Implement MM algorithm for the logistic regression

mm = function(sim_data, init, tol=1e-8, max_iter=100, verbose=FALSE){
  # Grab data
  y = sim_data[,1]
  x = cbind(rep(1, nrow(sim_data)), sim_data[,2])
  x0 = x[,1]
  x1 = x[,2]
  
  # Initialize the solution
  iter = 1
  beta0_history = beta1_history = tol_value = rep(NA, length.out = max_iter)
  
  beta0_history[iter] = init[1]
  beta1_history[iter] = init[2]
  tol_value[iter] = 10000
  
  # Define functions for uniroot
  beta0_f <- function(beta0_new, beta, x, y){
    beta0_old = beta[1]
    gradient = 0
    for(i in 1:length(y)){
      gradient = gradient - exp(sum(x[i,] %*% beta)) * x0[i] * exp(-2*x0[i]*beta0_old) /  (1+exp(sum(x[i,] %*% beta))) * exp(2*x0[i]*beta0_new) + y[i] * x0[i]
    }
    return(gradient)
  }
  
  beta1_f <- function(beta1_new, beta, x, y){
    beta1_old = beta[2]
    gradient = 0
    for(i in 1:length(y)){
      gradient = gradient - exp(x[i,] %*% beta) * x1[i] * exp(-2*x1[i]*beta1_old) /  (1+exp(x[i,] %*% beta)) * exp(2*x1[i]*beta1_new) + y[i] * x1[i]
    }
    return(gradient)
  }
  
  while (iter < max_iter) {
    beta = c(beta0_history[iter], beta1_history[iter])

    # Compute new beta
    beta0_uniroot = uniroot(beta0_f, c(-10, 10), beta=beta, x=x, y=y)
    beta0_new = beta0_uniroot$root
    gradient0 = beta0_uniroot$f.root
    
    beta1_uniroot = uniroot(beta1_f, c(-10, 10), beta=beta, x=x, y=y)
    beta1_new = beta1_uniroot$root
    gradient1 = beta1_uniroot$f.root

    beta0_history[iter+1] = beta0_new
    beta1_history[iter+1] = beta1_new

    tol_value[iter+1] = sqrt(gradient0^2 + gradient1^2)
    
    # Check stopping criterion
    if(tol_value[iter+1] < tol){
      # message("Converged in", iter, "iterations.\n")
      break
    }

    if(iter %% 100 == 0 & verbose==TRUE){
      message("Iteration: ", iter, " | Tol: ", tol_value[iter+1])
    }
    
    # Update iter
    iter = iter + 1
  }
  
  beta_new = c(beta0_new, beta1_new)
  
  return(
    list(
      solution = beta_new,
      iter=iter,
      history = cbind(beta0_history, beta1_history),
      tol_value = tol_value
    )
  )
}




## Scratch
# beta0_new = newtonRaphson(function(beta0_new) beta0_f(beta0_new, beta, x, y), x0=beta[1])$root
# beta1_new = newtonRaphson(function(beta1_new) beta1_f(beta1_new, beta, x, y), x0=beta[2])$root
