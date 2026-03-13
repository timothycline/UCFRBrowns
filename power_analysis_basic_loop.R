#******************************************************************************
# Brown trout simulation
#
# Scenarios: 1) mid survival, high migration
#            2) low survival, low migration - hypothesis
#
#******************************************************************************

library(tidyverse)
sim_data <- read.csv(file = "./simulation_input.csv")

out_list <- list()

n_reps <- 200
n_scenarios <- nrow(sim_data)



# Start simulation --------------------------------------------------------
for( s in 1:n_scenarios){
  # Inputs 
  n_tagged <- sim_data$n_tagged[s]
  
  survival <- sim_data$survival[s]
  migration <- sim_data$migration[s]
  p_efish <- sim_data$p_efish[s]
  p_scan <- sim_data$p_scan[s]
  
  # Detection probabilities ------------------------------------------------
  p_efish_mu <- p_efish
  p_efish_sd_logit <- 0.25
  
  p_efish_mean_logit <- qlogis(p_efish_mu)
  
  p_scan_mu <- p_scan
  p_scan_sd_logit <- 0.25
  
  p_scan_mean_logit <- qlogis(p_scan_mu)
  
  out_tibble <- 
    tibble(
      scenario = rep(s, n_reps),
      surv_mu = rep(NA, n_reps),
      surv_sd = rep(NA, n_reps),
      mig_mu = rep(NA, n_reps),
      mig_sd = rep(NA, n_reps)
    )
  
  for( i in 1:n_reps){
    
    # Simulate population -----------------------------------------------------
    # Biological process
    n_surv <- round(n_tagged * survival)
    n_stay <- round(n_surv * (1 - migration))
    n_mig <- n_surv - n_stay
    
    # Observation process
    n_recapped <- rbinom(n = 1, size = n_stay, prob = p_efish)
    n_scanned <- rbinom(n = 1, size = n_mig, prob = p_scan)
    
    
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
      n.thin = 2, 
      parallel = T
    )
    
    print(i)
    
    out_tibble$surv_mu[i] <- jags_out$mean$surv
    out_tibble$surv_sd[i] <- jags_out$sd$surv
    out_tibble$mig_mu[i] <- jags_out$mean$mig
    out_tibble$mig_sd[i] <- jags_out$sd$mig
  
  } # End iteration loop
  
  out_list[[s]] <- out_tibble
  
} # End scenario loop

all_data <- bind_rows(out_list)

write_rds(all_data, "./simulation_output.rds")

