#******************************************************************************
# Power analysis: Upper Clark Fork juvenile brown trout survival study
#
# System: Warm Springs -> Gembeck Rd -> Sager Ln (~11 river miles each)
#         3 tributaries (2 above Gembeck, 1 between Gembeck and Sager)
#         Fish tagged primarily in tributaries
#
# Key scientific questions:
#   1. Is mainstem survival very low?
#   2. Are fish emigrating from tribs permanently leaving the upper river?
#
# Empirical detection estimates (pilot work, 2023-2025):
#   p_trib   : 0.93  (trib mouth antennas; directional)
#   p_gem    : 0.76  (Gembeck Rd mainstem antenna)
#   p_ws     : 0.76  (Warm Springs mainstem antenna)
#   p_sag    : 0.76  (Sager Ln mainstem antenna)
#   p_efish  : 0.15-0.35 (barge/boat electrofishing)
#
# Survival expectation: 10-30% annual
#   -> per-occasion phi = annual_phi^(1/n_occ)
#
# Feasible tagging: 4,000-6,000 via spring electrofishing + screw traps
#
#******************************************************************************

library(tidyverse)
library(jagsUI)
library(here)
library(furrr)

# Parallel setup:
# Each JAGS rep uses 3 cores (parallel chains). With 10 logical cores available,
# run 3 scenarios simultaneously (3 scenarios × 3 chains = 9 cores, 1 free).
n_jags_chains <- 3
n_workers     <- floor(parallel::detectCores() / n_jags_chains)
plan(multisession, workers = n_workers)

# here() always resolves paths relative to the project root (.Rproj location)
# regardless of the current working directory
set.seed(123)

# ==============================================================================
# Pilot-data-derived prior hyperparameters (logit scale)
# ==============================================================================
# These come from the empirical detection estimates and are used to build
# informative priors. SD of 0.3 on the logit scale reflects moderate
# uncertainty around each point estimate (roughly ±10 percentage points
# at p = 0.76, ±5 percentage points at p = 0.93).
#
# For vague priors: mu = 0 (p = 0.5 on probability scale), sd = 5
# (so wide that the prior is essentially flat from 0 to 1).

make_prior <- function(p_central, sd_logit) {
  list(mu = qlogis(p_central), sd = sd_logit)
}

priors_informative <- list(
  p_trib  = make_prior(0.93, 0.30),
  p_gem   = make_prior(0.76, 0.30),
  p_ws    = make_prior(0.76, 0.30),
  p_sag   = make_prior(0.76, 0.30),
  p_efish = make_prior(0.25, 0.40)   # wider SD: more year-to-year variation
)

priors_vague <- list(
  p_trib  = make_prior(0.50, 5.0),
  p_gem   = make_prior(0.50, 5.0),
  p_ws    = make_prior(0.50, 5.0),
  p_sag   = make_prior(0.50, 5.0),
  p_efish = make_prior(0.50, 5.0)
)


# ==============================================================================
# True parameter values for simulation
# ==============================================================================

n_occ <- 5   # seasonal occasions (spring, early summer, late summer, fall, spring)

# Survival (converted from annual to per-occasion inside the scenario loop)
# Baseline: 20% annual -> 0.20^(1/5) = 0.725 per occasion

# Movement (per occasion, conditional on surviving)
psi_trib_out_true  <- 0.06    # trib emigration to mainstem
psi_main_back_true <- 0.04    # return from mainstem to trib
psi_main_down_true <- 0.04    # downstream emigration past Gembeck
psi_main_up_true   <- 0.01    # upstream emigration past Warm Springs (probably rare)
psi_sager_true     <- 0.40    # passage past Sager given alive in Gem-Sag section

stopifnot(psi_main_back_true + psi_main_down_true + psi_main_up_true <= 1)

# Detection (true values used in simulation)
p_trib_true  <- 0.93
p_gem_true   <- 0.76
p_ws_true    <- 0.76
p_sag_true   <- 0.76
p_efish_true <- 0.25

# Initial state distribution (mostly in tribs)
z_init_probs <- c(0.85, 0.15, 0.00, 0.00, 0.00, 0.00)


# ==============================================================================
# Scenario grid
# ==============================================================================
#
# Primary comparison: n_tagged (500 to 4000)
# Secondary: annual survival (10%, 20%, 30%)
# Tertiary: p_efish (bad year, typical, good year)
# Prior toggle: informative vs vague detection priors
#
# Scenarios 1-6: informative priors (using pilot data)
# Scenarios 7-9: vague priors on same baseline conditions
#               -> shows how much pilot data improves precision
#
# Scenario | n_tagged | phi_ann | p_efish | Priors       | Notes
# ---------|----------|---------|---------|--------------|-------------------
#    1     |    500   |   0.20  |  0.25   | informative  | small effort
#    2     |   1000   |   0.20  |  0.25   | informative  | baseline
#    3     |   2000   |   0.20  |  0.25   | informative  | 2x baseline
#    4     |   4000   |   0.20  |  0.25   | informative  | near-max effort
#    5     |   1000   |   0.10  |  0.25   | informative  | low survival
#    6     |   1000   |   0.30  |  0.25   | informative  | higher survival
#    7     |   1000   |   0.20  |  0.15   | informative  | poor efish year
#    8     |   1000   |   0.20  |  0.25   | vague        | baseline, no pilot data
#    9     |   2000   |   0.20  |  0.25   | vague        | 2x tags, no pilot data
#   10     |   1000   |   0.10  |  0.25   | vague        | low survival, no pilot data

scenarios <- tibble(
  scenario            = 1:10,
  n_tagged            = c( 500, 1000, 2000, 4000, 1000, 1000, 1000, 1000, 2000, 1000),
  phi_annual          = c(0.20, 0.20, 0.20, 0.20, 0.10, 0.30, 0.20, 0.20, 0.20, 0.10),
  p_efish_true_val    = c(0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.15, 0.25, 0.25, 0.25),
  use_inform_priors   = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE),
  label               = c(
    "500 tagged, inform. priors",
    "1000 tagged, inform. priors (baseline)",
    "2000 tagged, inform. priors",
    "4000 tagged, inform. priors",
    "Low survival (10%), inform. priors",
    "Higher survival (30%), inform. priors",
    "Poor efish (p=0.15), inform. priors",
    "1000 tagged, VAGUE priors",
    "2000 tagged, VAGUE priors",
    "Low survival (10%), VAGUE priors"
  )
)

n_reps <- 50


# ==============================================================================
# Simulation function
# ==============================================================================

simulate_fish <- function(n_tagged, phi_trib, phi_main, phi_lower,
                          psi_trib_out, psi_main_back, psi_main_down,
                          psi_main_up, psi_sager,
                          p_trib, p_gem, p_ws, p_sag, p_efish,
                          n_occ, z_init_probs) {

  n_states <- 6

  # Build transition matrix
  # States: 1=trib, 2=main WS-Gem, 3=main Gem-Sag, 4=past Sager, 5=past WS, 6=dead
  tr <- matrix(0, nrow = n_states, ncol = n_states)

  tr[1, 1] <- phi_trib * (1 - psi_trib_out)
  tr[1, 2] <- phi_trib * psi_trib_out
  tr[1, 6] <- 1 - phi_trib

  tr[2, 1] <- phi_main * psi_main_back
  tr[2, 2] <- phi_main * (1 - psi_main_back - psi_main_down - psi_main_up)
  tr[2, 3] <- phi_main * psi_main_down
  tr[2, 5] <- phi_main * psi_main_up
  tr[2, 6] <- 1 - phi_main

  tr[3, 3] <- phi_lower * (1 - psi_sager)
  tr[3, 4] <- phi_lower * psi_sager
  tr[3, 6] <- 1 - phi_lower

  tr[4, 4] <- 1   # past Sager: absorbing
  tr[5, 5] <- 1   # past WS: absorbing
  tr[6, 6] <- 1   # dead: absorbing

  # Simulate initial states
  z_init <- sample(1:n_states, size = n_tagged, replace = TRUE,
                   prob = z_init_probs)

  # Simulate hidden state trajectories (vectorized across fish)
  z <- matrix(NA_integer_, nrow = n_tagged, ncol = n_occ)
  z[, 1] <- z_init

  # Inverse-CDF categorical sampler: draws one category per row of a
  # probability matrix without any R-level fish loop
  draw_states <- function(prob_matrix) {
    u        <- runif(nrow(prob_matrix))
    cum_prob <- t(apply(prob_matrix, 1, cumsum))
    rowSums(cum_prob < u) + 1L
  }

  for(t in 2:n_occ) {
    z[, t] <- draw_states(tr[z[, t-1], ])
  }

  # Simulate observations (vectorized across fish)
  y_trib  <- matrix(0L, nrow = n_tagged, ncol = n_occ)
  y_gem   <- matrix(0L, nrow = n_tagged, ncol = n_occ)
  y_ws    <- matrix(0L, nrow = n_tagged, ncol = n_occ)
  y_sag   <- matrix(0L, nrow = n_tagged, ncol = n_occ)
  y_efish <- matrix(0L, nrow = n_tagged, ncol = n_occ)

  for(t in 2:n_occ) {
    prev <- z[, t-1]
    curr <- z[, t]

    # Identify which fish triggered each antenna this occasion
    trib_cross <- (prev == 1L & curr == 2L) | (prev == 2L & curr == 1L)
    gem_cross  <-  prev == 2L & curr == 3L
    ws_cross   <-  prev == 2L & curr == 5L
    sag_cross  <-  prev == 3L & curr == 4L
    in_main    <-  curr == 2L

    # Draw detections only for fish that could trigger each detector
    if(any(trib_cross)) y_trib[trib_cross, t] <- rbinom(sum(trib_cross), 1, p_trib)
    if(any(gem_cross))  y_gem[gem_cross,   t] <- rbinom(sum(gem_cross),  1, p_gem)
    if(any(ws_cross))   y_ws[ws_cross,     t] <- rbinom(sum(ws_cross),   1, p_ws)
    if(any(sag_cross))  y_sag[sag_cross,   t] <- rbinom(sum(sag_cross),  1, p_sag)
    if(any(in_main))    y_efish[in_main,   t] <- rbinom(sum(in_main),    1, p_efish)
  }

  list(z       = z,
       z_init  = z_init,
       y_trib  = y_trib,
       y_gem   = y_gem,
       y_ws    = y_ws,
       y_sag   = y_sag,
       y_efish = y_efish)
}


# ==============================================================================
# Power analysis loop
# ==============================================================================

out_list <- vector("list", nrow(scenarios))

# Wrap one scenario into a function so future_map can distribute across workers
run_scenario <- function(s) {

  n_tagged  <- scenarios$n_tagged[s]
  phi_occ   <- scenarios$phi_annual[s]^(1 / n_occ)
  p_efish_s <- scenarios$p_efish_true_val[s]
  priors    <- if(scenarios$use_inform_priors[s]) priors_informative else priors_vague

  out_tibble <- tibble(
    scenario      = s,
    n_tagged      = n_tagged,
    use_inform    = scenarios$use_inform_priors[s],
    rep           = 1:n_reps,
    phi_trib_mean = NA_real_,  phi_trib_sd  = NA_real_,
    phi_main_mean = NA_real_,  phi_main_sd  = NA_real_,
    psi_down_mean = NA_real_,  psi_down_sd  = NA_real_,
    psi_back_mean = NA_real_,  psi_back_sd  = NA_real_,
    p_gem_mean    = NA_real_,  p_gem_sd     = NA_real_,
    p_efish_mean  = NA_real_,  p_efish_sd   = NA_real_,
    converged     = logical(n_reps)
  )

  for(i in 1:n_reps) {

    sim <- simulate_fish(
      n_tagged      = n_tagged,
      phi_trib      = phi_occ,   phi_main  = phi_occ,  phi_lower = phi_occ,
      psi_trib_out  = psi_trib_out_true,
      psi_main_back = psi_main_back_true,
      psi_main_down = psi_main_down_true,
      psi_main_up   = psi_main_up_true,
      psi_sager     = psi_sager_true,
      p_trib        = p_trib_true,   p_gem = p_gem_true,
      p_ws          = p_ws_true,     p_sag = p_sag_true,
      p_efish       = p_efish_s,
      n_occ         = n_occ,
      z_init_probs  = z_init_probs
    )

    jags_data <- list(
      n_fish     = n_tagged,   n_occ   = n_occ,
      z_init     = as.integer(sim$z_init),
      y_trib     = sim$y_trib, y_gem   = sim$y_gem,
      y_ws       = sim$y_ws,   y_sag   = sim$y_sag,
      y_efish    = sim$y_efish,
      dir_alpha  = c(2, 2, 1),
      p_trib_mu  = priors$p_trib$mu,   p_trib_sd  = priors$p_trib$sd,
      p_gem_mu   = priors$p_gem$mu,    p_gem_sd   = priors$p_gem$sd,
      p_ws_mu    = priors$p_ws$mu,     p_ws_sd    = priors$p_ws$sd,
      p_sag_mu   = priors$p_sag$mu,    p_sag_sd   = priors$p_sag$sd,
      p_efish_mu = priors$p_efish$mu,  p_efish_sd = priors$p_efish$sd
    )

    z_inits      <- sim$z
    z_inits[, 1] <- NA

    jags_inits <- function() {
      list(
        z             = z_inits,
        phi_trib      = runif(1, 0.50, 0.95),
        phi_main      = runif(1, 0.50, 0.95),
        phi_lower     = runif(1, 0.50, 0.95),
        psi_trib_out  = runif(1, 0.02, 0.15),
        psi_main_move = runif(1, 0.05, 0.20),
        dir_main      = c(0.40, 0.40, 0.20),
        psi_sager     = runif(1, 0.20, 0.60),
        p_trib_logit  = rnorm(1, priors$p_trib$mu,  0.1),
        p_gem_logit   = rnorm(1, priors$p_gem$mu,   0.1),
        p_ws_logit    = rnorm(1, priors$p_ws$mu,    0.1),
        p_sag_logit   = rnorm(1, priors$p_sag$mu,   0.1),
        p_efish_logit = rnorm(1, priors$p_efish$mu, 0.1)
      )
    }

    jags_out <- tryCatch(
      jagsUI::jags(
        data       = jags_data,
        inits      = jags_inits,
        parameters = c("phi_trib", "phi_main", "phi_lower",
                       "psi_trib_out", "psi_main_back", "psi_main_down",
                       "psi_main_up", "psi_sager",
                       "p_trib", "p_gem", "p_ws", "p_sag", "p_efish"),
        model.file = here("survival_CFBrown.txt"),
        n.chains   = n_jags_chains,
        n.iter     = 12000,
        n.burnin   = 6000,
        n.thin     = 2,
        parallel   = TRUE,
        verbose    = FALSE
      ),
      error = function(e) NULL
    )

    if(!is.null(jags_out)) {
      out_tibble$phi_trib_mean[i] <- jags_out$mean$phi_trib
      out_tibble$phi_trib_sd[i]   <- jags_out$sd$phi_trib
      out_tibble$phi_main_mean[i] <- jags_out$mean$phi_main
      out_tibble$phi_main_sd[i]   <- jags_out$sd$phi_main
      out_tibble$psi_down_mean[i] <- jags_out$mean$psi_main_down
      out_tibble$psi_down_sd[i]   <- jags_out$sd$psi_main_down
      out_tibble$psi_back_mean[i] <- jags_out$mean$psi_main_back
      out_tibble$psi_back_sd[i]   <- jags_out$sd$psi_main_back
      out_tibble$p_gem_mean[i]    <- jags_out$mean$p_gem
      out_tibble$p_gem_sd[i]      <- jags_out$sd$p_gem
      out_tibble$p_efish_mean[i]  <- jags_out$mean$p_efish
      out_tibble$p_efish_sd[i]    <- jags_out$sd$p_efish
      out_tibble$converged[i]     <- all(unlist(jags_out$Rhat) < 1.1, na.rm = TRUE)
    }

  } # end rep loop

  # Checkpoint: save after each scenario completes
  saveRDS(out_tibble, here(paste0("checkpoint_scenario_", s, ".rds")))
  cat("Scenario", s, "complete.\n")

  out_tibble
}

# Run all scenarios in parallel across workers
# With 10 cores and 3 chains/rep: 3 scenarios run simultaneously
out_list <- future_map(1:nrow(scenarios), run_scenario,
                       .options = furrr_options(seed = TRUE))

all_data <- bind_rows(out_list)
write_rds(all_data, here("simulation_output_CFBrown.rds"))
cat("\nDone. Output saved to simulation_output_CFBrown.rds\n")
