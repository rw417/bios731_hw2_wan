glm_function <- function(sim_data){
  if(class(sim_data)[1] != "data.frame"){
    sim_data = data.frame(sim_data)
  }
  glm_model <- glm(y~x, data = sim_data, family = binomial)
  return(list(
    solution = as.numeric(glm_model$coef),
    iter = glm_model$iter
    )
  )
}
