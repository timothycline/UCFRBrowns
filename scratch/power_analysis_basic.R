#******************************************************************************
# Brown trout simulation
#
# Scenarios: 1) mid survival, high migration
#            2) low survival, low migration - hypothesis
#
#******************************************************************************
# Inputs ----------------------------------------------------------------
n_tagged <- 500 # 100, 500, 1000 - hyp - 500


survival <- 0.3 # 0.1, 0.3 - hyp - 0.3
migration <- 0.05 # 0.05, 0.5, hyp - 0.05
p_efish <- 0.3 # 0.3, 0.45, hyp - 0.3
p_scan <- 0.8 # 0.8


# Detection probabilities -------------------------------------------------
p_efish_mu <- p_efish
p_efish_sd_logit <- 0.3

p_efish_mean_logit <- qlogis(p_efish_mu)

p_scan_mu <- p_scan
p_scan_sd_logit <- 0.3

p_scan_mean_logit <- qlogis(p_scan_mu)


par(mfrow = c(2, 1), mar = c(4.5, 4.5, 1, 1))

hist(plogis(rnorm(100000, mean = p_efish_mean_logit, sd = p_efish_sd_logit)),
     xlab = "E-fish detection prob.", main = "",
     xlim = c(0, 1), breaks = seq(0, 1, 0.01))

hist(plogis(rnorm(100000, mean = p_scan_mean_logit, sd = p_scan_sd_logit)),
     xlab = "PIT array detection prob.", main = "",
     xlim = c(0, 1), breaks = seq(0, 1, 0.01))

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
  n.thin = 2
)

jags_out
  

