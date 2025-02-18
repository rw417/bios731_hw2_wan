# Implement MM algorithm for the logistic regression

mm = function(sim_data, init, tol=1e-12, max_iter=100){
  # Grab data
  y = sim_data[,1]
  x = cbind(rep(1, nrow(sim_data)), sim_data[,2])
  
  # Initialize beta values
  beta0=init[1]
  beta1=init[2]
  
  beta0_history = beta1_history = rep(NA, length.out = max_iter)
  
  # Initialize the solution
  iter = 1
  beta0_history[iter] = beta0
  beta1_history[iter] = beta1
  
  while (iter < max_iter-1) {
    beta = c(beta0_history[iter], beta1_history[iter])
    
    # Compute beta0
    
    
    # Check stopping criterion
    if(sqrt(sum(gradient^2)) < tol){
      message("Converged in", iter, "iterations.\n")
      break
    }
    
    # Update iter
    iter = iter + 1
  }
  
  return(
    list(
      solution = beta_new,
      history = cbind(beta0_history, beta1_history)
    )
  )
}