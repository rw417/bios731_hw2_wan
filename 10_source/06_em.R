# EM Algorithm for censored exponential data
em_exp = function(dta, survival = "time", censor = "status", lambda_init, tol = 1e-12, max_iter = 200){
  lambda = lambda_init
  delta = dta[[censor]]
  y = dta[[survival]]
  
  n = dim(dta)[1]
  iter = 1
  tol_criteria = 10000
  
  # define vectors to store elements of interest
  observed_ll  = rep(NA, length = max_iter)
  
  # message(paste0("delta: ", delta))
  # message(paste0("y: ", y))

  while(iter <= max_iter  & tol_criteria > tol){
    
    ###############################################################
    ## E-step
    ###############################################################
    
    z = y + 1/lambda
    
    ###############################################################
    ## M-step
    ###############################################################
    
    lambda = n/sum(delta * y + (1 - delta)*z)
    
    ###############################################################

    ## define log likelihood at current iteration
    observed_ll[iter] = n*log(lambda) - lambda * (sum(delta * y + (1-delta)*z))
    
    ## check convergence criteria and increase iteration
    if(iter > 1){
      tol_criteria = abs(observed_ll[iter] - observed_ll[iter-1])
    }
    # message(paste0("iteration: ", iter, "; lambda:", lambda, "; ll: ", round(observed_ll[iter], 4)), "; tol: ", round(tol_criteria, 4))
    iter = iter + 1
  }
  
  ### return parameters of interest
  return(
    list(solution = log(1/lambda),
       iterations = iter,
       observed_ll = observed_ll[1:iter]
     )
  )
}


