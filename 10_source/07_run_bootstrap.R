####################
# Run bootstrap on simulated data
run_bootstrap <- function(sim_data, nboot, seed = NULL, verbose=FALSE) {
  if (!is.null(seed)) set.seed(seed)
  
  # Variable to record the rows of data to sample
  n_bs <- nrow(sim_data)
  
  # Create a list of 5 matrices to store results
  bs_results <- lapply(1:4, function(i) matrix(NA, nrow=nboot, ncol=4))
  names(bs_results) <- c("newton", "mm", "glm", "optim")

  ####################
  # Start of Bootstrap
  for (b in 1:nboot) {
    ####################
    # sample with replacement from sim_data
    bs_idx <- sample(1:n_bs, size = n_bs, replace = TRUE)
    bs_sample <- sim_data[bs_idx, ]  # Direct indexing

    ####################
    # apply estimation methods to the data
    # Newton
    time_start <- proc.time()[3] # for timing
    newton_results <- newton(bs_sample, init=c(0.5,0.5), tol=1e-08, max_iter=100)
    time_newton <- proc.time()[3] - time_start
    bs_results$newton[b,] <- c(newton_results$solution, newton_results$iter, time_newton)
    
    # MM
    time_start <- proc.time()[3] # for timing
    mm_results <- mm(bs_sample, init=c(0.5,0.5), tol=1e-08, max_iter=100)
    time_mm <- proc.time()[3] - time_start
    bs_results$mm[b,] <- c(mm_results$solution, mm_results$iter, time_mm)
    
    # GLM
    time_start <- proc.time()[3] # for timing
    glm_results <- glm_function(bs_sample)
    time_glm <- proc.time()[3] - time_start
    bs_results$glm[b,] <- c(glm_results$solution, glm_results$iter, time_glm)
    
    # Optim
    time_start <- proc.time()[3] # for timing
    optim_results <- optim_function(bs_sample, init=c(0.5,0.5))
    time_optim <- proc.time()[3] - time_start
    bs_results$optim[b,] <- c(optim_results$solution, optim_results$iter, time_optim)
    
    if(b %% 10 == 0 & verbose==TRUE) {
      message(paste("Completed bootstrap iteration: ", b))
      flush.console()
    }
  }
  return(bs_results)
}
