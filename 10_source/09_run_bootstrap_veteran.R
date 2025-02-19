####################
# Run bootstrap on the veteran dataset
run_bootstrap_veteran <- function(nboot, seed = NULL, verbose=FALSE) {
  if (!is.null(seed)) set.seed(seed)
  
  # Load data
  library(survival)
  
  # Variable to record the rows of data to sample
  n_bs <- nrow(veteran)
  
  # Create a list of 5 matrices to store results
  bs_results <- lapply(1:1, function(i) matrix(NA, nrow=nboot, ncol=3))
  names(bs_results) <- c("em")
  
  ####################
  # Start of Bootstrap
  for (b in 1:nboot) {
    ####################
    # sample with replacement from veteran data
    bs_idx <- sample(1:n_bs, size = n_bs, replace = TRUE)
    bs_sample <- veteran[bs_idx, ]  # Direct indexing
    
    ####################
    # apply estimation methods to the data
    # em
    time_start <- proc.time()[3] # for timing
    em_results <- em_exp(dta=bs_sample, survival="time", censor="status", lambda_init=0.5, tol=1e-8, max_iter=100)
    time_em <- proc.time()[3] - time_start
    bs_results$em[b,] <- c(em_results$solution, em_results$iter, time_em)
    
    # # AFT
    # time_start <- proc.time()[3] # for timing
    # aft_results <- survreg(Surv(time, status) ~ 1, data=veteran, dist="exponential")
    # time_aft <- proc.time()[3] - time_start
    # bs_results$aft[b,] <- c(aft_results$solution, aft_results$iter, time_aft)

    if(b %% 10 == 0 & verbose==TRUE) {
      message(paste("Completed bootstrap iteration: ", b))
      flush.console()
    }
  }
  return(bs_results)
}
