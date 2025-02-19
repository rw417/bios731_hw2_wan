####################################################################
# Robert Wan
# January 2025
#
# This file produces simulations for linear regression under different data
# generating scenarios
####################################################################

library(stats)
library(here)

###############################################################
## define or source functions used in code below
###############################################################

# # Get all R script files in the directory
# script_files <- list.files(here::here("10_source"), pattern = "\\.R$", full.names = TRUE)
# 
# # Source all R scripts
# lapply(script_files, source)

source(here::here("10_source", "01_simulate_data.R"))
source(here::here("10_source", "02_newton.R"))
source(here::here("10_source", "03_mm.R"))
source(here::here("10_source", "04_glm.R"))
source(here::here("10_source", "05_optim.R"))
source(here::here("10_source", "06_em.R"))
source(here::here("10_source", "07_run_bootstrap.R"))

###############################################################
## set simulation design elements
###############################################################

n <- 200
nboot <- 200

###############################################################
## Run simulations
###############################################################
# Since we are only running 1 simulation, we only need to run 
# the bootstrap for this study.
# No extra simulation needed.

set.seed(0219)

# Simulate data
sim_data <- get_simdata(n=n, beta0=1, beta1=0.3)

# Run bootstrap
bs_results <- run_bootstrap(sim_data, nboot, verbose=TRUE)

# Save Results
saveRDS(bs_results, file = here::here("30_results", "bs_results.RDS"))

###############################################################
## EM on the veterans dataset
###############################################################
