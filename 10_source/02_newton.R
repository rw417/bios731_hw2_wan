# Implement Newton's method to perform gradient descent
newton = function(sim_data,init,tol=1e-8, max_iter=100, verbose=FALSE){
  # Grab data
  y = sim_data[,1]
  x = cbind(rep(1, nrow(sim_data)), sim_data[,2])
  
  # Initialize beta values
  beta0=init[1]
  beta1=init[2]
  
  beta0_history = beta1_history = tol_value = rep(NA, length.out = max_iter)
  
  # Initialize the solution
  iter = 1
  beta0_history[iter] = beta0
  beta1_history[iter] = beta1
  tol_value[iter] = 10000
  
  while (iter < max_iter) {
    beta = c(beta0_history[iter], beta1_history[iter])
    
    # Compute the gradient
    gradient = rep(0, length(beta))
    for(j in 1:length(beta)){
      gradient[j] = sum((y - exp(x %*% beta)/(1+exp(x %*% beta))) * x[,j])
    }

    # Compute the hessian
    hessian = matrix(0, nrow = length(beta), ncol = length(beta))
    for(j in 1:length(beta)){
      for(k in 1:length(beta)){
        hessian[j,k] = -1 * sum(exp(x %*% beta)/(1+exp(x %*% beta))^2 * x[,j]*x[,k])
      }
    }
    
    # Solve Newton
    beta_new = beta - solve(hessian) %*% gradient
    beta0_history[iter+1] = beta_new[1]
    beta1_history[iter+1] = beta_new[2]
    tol_value[iter+1] = sqrt(sum(gradient^2))
    
    # message("Iteration ", iter, " complete \n")
    # message("Gradient: ", gradient[1], " ", gradient[2], "\n")
    # message("Beta: ", beta[1], " ", beta[2], "\n")
    # message("Beta_new: ", beta_new[1], " ", beta_new[2], "\n")

    # Check stopping criterion
    if(tol_value[iter+1] < tol){
      # message("Converged in", iter, "iterations.\n")
      break
    }
    
    # Update iter
    iter = iter + 1
  }
  
  return(
    list(
      solution = as.numeric(beta_new),
      iter = iter,
      history = cbind(beta0_history, beta1_history),
      tol_value = tol_value
      )
    )
}
  
