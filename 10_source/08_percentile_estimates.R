####################
# Calculate the Wald standard error and percentile-based CI for all bootstrap runs
percentile_estimates <- function(result, true_value, alpha = 0.05) {
  # Estimate standard error of beta hat from bootstrap results
  se_b <- sd(result, na.rm=TRUE)
  
  # Estimate CI and coverage
  ci_lower <- quantile(result, alpha / 2, na.rm=TRUE)
  ci_upper <- quantile(result, 1 - alpha / 2, na.rm=TRUE)
  coverage <- (true_value > ci_lower & true_value < ci_upper)
  
  bs_estimates <- c(mean(result), se_b, ci_lower, ci_upper, coverage)
  names(bs_estimates) <- c("mean", "se_b", "ci_lower", "ci_upper", "coverage")

  return(bs_estimates)
}
