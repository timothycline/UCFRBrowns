#******************************************************************************
# Brown trout simulation
#******************************************************************************

n_tagged <- 1000
n_states <- 3

survival <- 0.4
migration <- 0.4
p_efish <- 0.4
p_scan <- 0.9


# State matrix
state_matrix <- matrix(NA, nrow = n_states, ncol = n_states)

state_matrix[1,1] <- survival * (1 - migration)
state_matrix[1,2] <- survival * migration
state_matrix[1,3] <- 1 - survival
state_matrix[2,1] <- 0 # Assume emmigration is permenant
state_matrix[2,2] <- 1
state_matrix[2,3] <- 0
state_matrix[3,1] <- 0 # Dead fish must stay dead
state_matrix[3,2] <- 0
state_matrix[3,3] <- 1


# Detection matrix
obs_matrix <- matrix(NA, nrow = n_states, ncol = n_states)

obs_matrix[1,1] <- p_efish
obs_matrix[1,2] <- 0
obs_matrix[1,3] <- 1 - p_efish
obs_matrix[2,1] <- 0 
obs_matrix[2,2] <- p_scan
obs_matrix[2,3] <- 1 - p_scan
obs_matrix[3,1] <- 0 
obs_matrix[3,2] <- 0
obs_matrix[3,3] <- 1

# Transition states
z <- matrix(NA, nrow = n_tagged, ncol = 2)
z[,1] <- 1

for( i in 1:n_tagged){
  z[i,2] <- which(rmultinom(1, 1, state_matrix[z[i,1], ]) == 1)
}

# Observation process
ch <- matrix(NA, nrow = n_tagged, ncol = 2)
ch[,1] <- 1

for( i in 1:n_tagged){
  ch[i,2] <- which(rmultinom(1, 1, obs_matrix[z[i,2], ]) == 1)
}

# 


# Alternative model -------------------------------------------------------
n_tagged <- 1000

survival <- 0.4
migration <- 0.4
p_efish <- 0.4
p_scan <- 0.9

# Biological process
n_surv <- rbinom(n = 1, size = n_tagged, prob = survival)
n_stay <- rbinom(n = 1, size = n_surv, prob = (1 - migration))
n_mig <- n_surv - n_stay

# Observation process
n_recapped <- rbinom(n = 1, size = n_stay, prob = p_efish)
n_scanned <- rbinom(n = 1, size = n_mig, prob = p_efish)



# Prepare JAGS ------------------------------------------------------------

# Chains
ni <- 50000
nb <- 25000
nt <- 1
nc <- 3

# Data list
jags.data <- list( n_recapped = n_recapped,
                   n_scanned = n_scanned,
                   n_tagged = n_tagged)

# Initials

jags.inits <- function(){
  list(
    surv = runif(1, 0.1, 0.9),
    mig = runif(1, 0.1, 0.9)
  )
}


# parameters to track
jags.parm <- c("surv", "mig")

# run the model
jags_out <- jagsUI::jags(
  jags.data,
  jags.inits,
  jags.parm,
  "./jags/Nmixture_basic.txt",
  n.chains=nc,
  n.iter=ni,
  n.burnin=nb,
  n.thin=nt
)

jags_out
  













