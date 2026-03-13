#******************************************************************************
# Brown trout simulation output summary
#******************************************************************************

# Load data and libraries -------------------------------------------------
library(tidyverse)
library(forcats)
library(cowplot)
sim_out <- read_rds("./simulation_output.rds")
sim_data <- read.csv(file = "./simulation_input.csv")

n_scenarios <- nrow(sim_data)

# Calculate summary stats -------------------------------------------------
sim_data$bias_surv <- NA
sim_data$bias_mig <- NA
sim_data$mean_se_surv <- NA
sim_data$mean_se_mig <- NA
sim_data$rmse_surv <- NA
sim_data$rmse_mig <- NA

for( s in 1:n_scenarios){
  out_tmp <- sim_out |> filter(scenario == s)
  
  sim_data$bias_surv[s] <- mean(out_tmp$surv_mu - sim_data$survival[s])
  sim_data$bias_mig[s] <- mean(out_tmp$mig_mu - sim_data$migration[s])
  sim_data$mean_se_surv[s] <- mean(out_tmp$surv_sd)
  sim_data$mean_se_mig[s] <- mean(out_tmp$mig_sd)
  sim_data$rmse_surv[s] <- 
    sqrt(mean((out_tmp$surv_mu - sim_data$survival[s])^2))
  sim_data$rmse_mig[s] <- 
    sqrt(mean((out_tmp$mig_mu - sim_data$migration[s])^2))
  
}

write_csv(sim_data, "./Simulation_summary.csv")

mean(sim_data$bias_mig)

sim_out <- 
  sim_out |> 
  mutate(scenario = factor(scenario)) |> 
  mutate(scenario = fct_recode(scenario,
                               "1) Baseline" = "1",
                               "2) 100 tagged" = "2",
                               "3) 1000 tagged" = "3",
                               "4) Survival = 0.1" = "4",
                               "5) Migration = 0.5" = "5",
                               "6) Efish p = 0.45" = "6"))



surv_plot <- 
  sim_out |> 
  ggplot(aes(x = scenario, y = surv_mu)) +
  geom_boxplot(fill = "grey70") +
  xlab("Scenario") + 
  ylab("Survival probability") + 
  ylim(c(0, 0.6)) + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 11)) +
  geom_segment(aes(x = 0.5, xend = 3.5, y = 0.3), 
               col = "blue", linewidth = 1.5) +
  geom_segment(aes(x = 3.5, xend = 4.5, y = 0.1), 
               col = "blue", linewidth = 1.5) +
  geom_segment(aes(x = 4.5, xend = 6.5, y = 0.3), 
               col = "blue", linewidth = 1.5) +
  geom_boxplot(fill = "grey70")
  

mig_plot <- 
  sim_out |> 
  ggplot(aes(x = scenario, y = mig_mu)) +
  geom_boxplot(fill = "grey70") +
  xlab("Scenario") + 
  ylab("Migration probability") + 
  ylim(c(0, 0.6)) + 
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 1, 
                                   hjust = 1, 
                                   size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 13)) +
  geom_segment(aes(x = 0.5, xend = 4.5, y = 0.05), 
               col = "blue", linewidth = 1.5) +
  geom_segment(aes(x = 4.5, xend = 5.5, y = 0.5), 
               col = "blue", linewidth = 1.5) +
  geom_segment(aes(x = 5.5, xend = 6.5, y = 0.05), 
               col = "blue", linewidth = 1.5) +
  geom_boxplot(fill = "grey70")


combined_plot <- 
  plot_grid(surv_plot, mig_plot, 
            ncol = 1, align = "v", rel_heights = c(0.7, 1))

ggsave("precision_plot.tiff", plot = combined_plot, width = 5, height = 8, dpi = 300)




