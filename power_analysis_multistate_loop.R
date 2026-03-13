#******************************************************************************
# Multi-state brown trout power analysis
#
# Study design:
#   - 5 stationary PIT antennas divide the study reach into 6 spatial zones
#   - Fish tagged across multiple reaches (mixed tagging locations)
#   - Two-way movement (upstream and downstream)
#   - Seasonal electrofishing surveys (5 occasions per year)
#   - Antennas running continuously; detection summarized per occasion
#
# Model: Hidden Markov Model (multi-state capture-recapture)
#   States 1-6: alive in reach 1 (most upstream) through reach 6
#   State 7:    dead
#
# Key parameters estimated:
#   phi[r]       reach-specific survival per occasion
#   psi_down[r]  probability of downstream movement per occasion
#   psi_up[r]    probability of upstream movement per occasion
#   p_scan[k]    antenna k detection probability
#   p_efish      electrofishing detection probability
#
#******************************************************************************

library(tidyverse)
library(jagsUI)

set.seed(123)

# ==============================================================================
# True parameter values for simulation (biologically plausible baseline)
# ==============================================================================

n_reaches <- 6
n_ant     <- 5
n_occ     <- 5   # seasonal occasions (e.g., spring, early summer, late summer,
                  # fall, following spring)

# Reach-specific survival per occasion
# 0.85/occasion ≈ 0.44 cumulative survival over 5 occasions
phi_true <- rep(0.85, n_reaches)

# Movement: probability of moving at all per occasion
# Slightly lower in the terminal downstream reach
psi_move_true <- c(0.25, 0.25, 0.25, 0.25, 0.25, 0.15)

# Of movement, proportion going downstream (interior reaches 2-5)
prop_down_true <- c(NA, 0.65, 0.65, 0.65, 0.65, NA)

# Derived movement probabilities (enforce boundary conditions)
psi_down_true <- numeric(n_reaches)
psi_up_true   <- numeric(n_reaches)

psi_down_true[1] <- psi_move_true[1]   # reach 1: all movement downstream
psi_up_true[1]   <- 0

for(r in 2:5) {
  psi_down_true[r] <- psi_move_true[r] * prop_down_true[r]
  psi_up_true[r]   <- psi_move_true[r] * (1 - prop_down_true[r])
}

psi_down_true[6] <- 0                  # reach 6: all movement upstream
psi_up_true[6]   <- psi_move_true[6]

# Detection probabilities
p_scan_true  <- rep(0.85, n_ant)
p_efish_true <- 0.30

# Initial reach distribution for tagged fish
# Reflects mixed tagging locations (mostly upper reaches)
z_init_probs <- c(0.40, 0.30, 0.20, 0.07, 0.03, 0.00)


# ==============================================================================
# Simulation function
# ==============================================================================

simulate_fish <- function(n_tagged, phi, psi_down, psi_up, p_scan, p_efish,
                          n_occ, z_init_probs) {

  n_reaches  <- length(phi)
  n_ant      <- length(p_scan)
  dead_state <- n_reaches + 1

  # Build transition matrix
  tr <- matrix(0, nrow = dead_state, ncol = dead_state)

  for(r in 1:n_reaches) {
    pd   <- psi_down[r]
    pu   <- psi_up[r]
    stay <- 1 - pd - pu

    if(r > 1)         tr[r, r - 1]    <- phi[r] * pu    # move upstream
                      tr[r, r]         <- phi[r] * stay  # stay
    if(r < n_reaches) tr[r, r + 1]    <- phi[r] * pd    # move downstream
                      tr[r, dead_state] <- 1 - phi[r]   # die
  }
  tr[dead_state, dead_state] <- 1  # absorbing dead state

  # Simulate initial states from tagging location distribution
  z_init <- sample(1:n_reaches, size = n_tagged, replace = TRUE,
                   prob = z_init_probs)

  # Simulate hidden state trajectories
  z <- matrix(NA_integer_, nrow = n_tagged, ncol = n_occ)
  z[, 1] <- z_init

  for(t in 2:n_occ) {
    for(i in 1:n_tagged) {
      z[i, t] <- which(rmultinom(1, 1, tr[z[i, t-1], ]) == 1)
    }
  }

  # Simulate observations (occasions 2 through n_occ)
  y_efish <- matrix(0L, nrow = n_tagged, ncol = n_occ)
  y_ant   <- array(0L,  dim = c(n_tagged, n_occ, n_ant))

  for(t in 2:n_occ) {
    for(i in 1:n_tagged) {
      s <- z[i, t]

      # Electrofishing detects fish in reaches 1 through n_reaches-1
      if(s >= 1 && s <= (n_reaches - 1)) {
        y_efish[i, t] <- rbinom(1, 1, p_efish)
      }

      # Antenna k detects fish currently in reach k+1
      for(k in 1:n_ant) {
        if(s == (k + 1)) {
          y_ant[i, t, k] <- rbinom(1, 1, p_scan[k])
        }
      }
    }
  }

  list(z = z, z_init = z_init, y_efish = y_efish, y_ant = y_ant)
}


# ==============================================================================
# Scenario grid
# ==============================================================================
#
# Baseline: 500 fish, all other parameters at true values above
# Vary one factor at a time to isolate sensitivity
#
# Scenario | n_tagged | Note
# ---------|----------|------------------------------
#    1     |    100   | Small sample
#    2     |    250   | Moderate sample
#    3     |    500   | Baseline
#    4     |   1000   | Large sample
#    5     |    500   | Low survival (phi = 0.70)
#    6     |    500   | Low antenna detection (p_scan = 0.60)
#    7     |    500   | Low efish detection (p_efish = 0.15)

scenarios <- tibble(
  scenario   = 1:7,
  n_tagged   = c(100, 250, 500, 1000, 500, 500, 500),
  phi_scalar = c(0.85, 0.85, 0.85, 0.85, 0.70, 0.85, 0.85),
  p_scan_val = c(0.85, 0.85, 0.85, 0.85, 0.85, 0.60, 0.85),
  p_efish_val= c(0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.15),
  label      = c("100 tagged", "250 tagged", "500 tagged (baseline)",
                 "1000 tagged", "Low survival (phi=0.70)",
                 "Low antenna p (0.60)", "Low efish p (0.15)")
)

n_reps <- 50   # replications per scenario


# ==============================================================================
# Power analysis loop
# ==============================================================================

out_list <- vector("list", nrow(scenarios))

for(s in 1:nrow(scenarios)) {

  n_tagged    <- scenarios$n_tagged[s]
  phi_s       <- rep(scenarios$phi_scalar[s], n_reaches)
  p_scan_s    <- rep(scenarios$p_scan_val[s], n_ant)
  p_efish_s   <- scenarios$p_efish_val[s]

  cat("\n=== Scenario", s, ":", scenarios$label[s], "===\n")

  out_tibble <- tibble(
    scenario      = s,
    n_tagged      = n_tagged,
    rep           = 1:n_reps,
    phi_mean      = vector("list", n_reps),   # vector of length 6 per rep
    phi_sd        = vector("list", n_reps),
    psi_down_mean = vector("list", n_reps),
    psi_down_sd   = vector("list", n_reps),
    psi_up_mean   = vector("list", n_reps),
    psi_up_sd     = vector("list", n_reps),
    converged     = logical(n_reps)
  )

  for(i in 1:n_reps) {

    cat("  rep", i, "\n")

    # --- Simulate data ---
    sim <- simulate_fish(
      n_tagged     = n_tagged,
      phi          = phi_s,
      psi_down     = psi_down_true,
      psi_up       = psi_up_true,
      p_scan       = p_scan_s,
      p_efish      = p_efish_s,
      n_occ        = n_occ,
      z_init_probs = z_init_probs
    )

    # --- JAGS data ---
    jags_data <- list(
      n_fish  = n_tagged,
      n_occ   = n_occ,
      z_init  = as.integer(sim$z_init),
      y_efish = sim$y_efish,
      y_ant   = sim$y_ant
    )

    # --- Initial values ---
    # Use true latent states for efficient MCMC convergence
    # (acceptable practice in power analysis simulations)
    z_inits       <- sim$z
    z_inits[, 1]  <- NA   # first occasion is a logical node; do not initialize

    jags_inits <- function() {
      list(
        z             = z_inits,
        phi           = runif(6, 0.60, 0.95),
        psi_move_raw  = runif(6, 0.10, 0.40),
        # prop_down_raw is only defined for reaches 2-5 (indices 2:5)
        # Provide NA at positions 1 and 6 so jagsUI skips those nodes
        prop_down_raw = c(NA, runif(4, 0.40, 0.80), NA),
        p_scan        = runif(5, 0.60, 0.95),
        p_efish       = runif(1, 0.10, 0.50)
      )
    }

    # --- Parameters to monitor ---
    jags_parm <- c("phi", "psi_down", "psi_up", "p_scan", "p_efish")

    # --- Run JAGS ---
    jags_out <- tryCatch(
      jagsUI::jags(
        data       = jags_data,
        inits      = jags_inits,
        parameters = jags_parm,
        model.file = "./survival_multistate.txt",
        n.chains   = 3,
        n.iter     = 12000,
        n.burnin   = 6000,
        n.thin     = 2,
        parallel   = TRUE,
        verbose    = FALSE
      ),
      error = function(e) {
        cat("    JAGS error:", conditionMessage(e), "\n")
        NULL
      }
    )

    if(!is.null(jags_out)) {
      out_tibble$phi_mean[[i]]      <- jags_out$mean$phi
      out_tibble$phi_sd[[i]]        <- jags_out$sd$phi
      out_tibble$psi_down_mean[[i]] <- jags_out$mean$psi_down
      out_tibble$psi_down_sd[[i]]   <- jags_out$sd$psi_down
      out_tibble$psi_up_mean[[i]]   <- jags_out$mean$psi_up
      out_tibble$psi_up_sd[[i]]     <- jags_out$sd$psi_up
      # Flag convergence: Rhat < 1.1 for all monitored parameters
      out_tibble$converged[i] <- all(unlist(jags_out$Rhat) < 1.1, na.rm = TRUE)
    }

  } # end rep loop

  out_list[[s]] <- out_tibble

} # end scenario loop

all_data <- bind_rows(out_list)
write_rds(all_data, "./simulation_output_multistate.rds")
cat("\nDone. Output saved to simulation_output_multistate.rds\n")
