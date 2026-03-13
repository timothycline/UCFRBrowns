#******************************************************************************
# Brown trout simulation
#******************************************************************************
# Inputs ----------------------------------------------------------------
n_tagged <- 500

survival <- 0.5
migration <- 0.8
p_efish <- 0.4
p_scan <- 0.7


# Simulate population -----------------------------------------------------
# Biological process
n_surv <- rbinom(n = 1, size = n_tagged, prob = survival)
n_stay <- rbinom(n = 1, size = n_surv, prob = (1 - migration))
n_mig <- n_surv - n_stay

# Observation process
n_recapped <- rbinom(n = 1, size = n_stay, prob = p_efish)
n_scanned <- rbinom(n = 1, size = n_mig, prob = p_scan)

# detection probability 
n_mig <- rbinom(n = n_mig, size = 2, prob = p_scan)


# Run JAGS model ----------------------------------------------------------
# Data list
jags_data <- list(n_recapped = n_recapped,
                  n_scanned = n_scanned,
                  n_tagged = n_tagged,
                  p_efish_mean_logit = p_efish_mean_logit,
                  p_efish_sd_logit = p_efish_sd_logit,
                  p_scan_mean_logit = p_scan_mean_logit,
                  p_scan_sd_logit = p_scan_sd_logit)

# Initials

jags_inits <- function(){
  list(
    surv = runif(1, 0.1, 0.9),
    mig = runif(1, 0.1, 0.9)
  )
}


# parameters to track
jags_parm <- c("surv", "mig")

# run the model
jags_out <- jagsUI::jags(
  jags_data,
  jags_inits,
  jags_parm,
  "./survival_basic.txt",
  n.chains = 3,
  n.iter = 10000,
  n.burnin = 5000,
  n.thin = 2
)

jags_out
  

