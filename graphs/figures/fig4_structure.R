### ============================================================================
### Disentangling length-dependent and length-independent variations
### in survival and maturation in Atlantic salmon
### ----------------------------------------------------------------------------
### File:    fig4_structure.R
### Purpose: Builds Figure 4.
### Author:  ©BOULAIRE Eliot, NEVOUX Marie & RIVOT Etienne
### Version: 01/09/2026
### ============================================================================

## -----------------------------------------------------------------------------
## 0. Setup
## -----------------------------------------------------------------------------
rm(list = ls())
set.seed(123)

pkginstall <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    } else {
      library(pkg, character.only = TRUE)
    }
  }
}
pkginstall(c("dplyr", "tidyr", "ggplot2", "coda", "nimble", "ggridges", "qs"))

## -----------------------------------------------------------------------------
## 1. Define project settings
## -----------------------------------------------------------------------------
projects <- "M1"

## -----------------------------------------------------------------------------
## 2. Import MCMC object and functions
## -----------------------------------------------------------------------------
source("functions/nimblefunctions/nf_l.R")
Cnf_l <- compileNimble(nf_l)
source("functions/nimblefunctions/nf_pi.R")
Cnf_pi <- compileNimble(nf_pi)

load_mcmc <- function(project) {
  MCMC <- qread(paste0("saves/", project, "/MCMC.qs"))
  MCMC %>%
    seq_along() %>%
    lapply(function(i) mcmc(MCMC[[i]]$samples)) %>%
    mcmc.list()
}
MCMC_samples <- load_mcmc(projects)
MCMC_matrix <- as.matrix(MCMC_samples)

const <- qread(file = "data/realdata/const.qs")

## -----------------------------------------------------------------------------
## 3. Survival
## -----------------------------------------------------------------------------
# Parameters
lmin <- const$Min_L1
lmax <- const$Max_L1
MCMC_Mu_L1 <- MCMC_matrix[, grep("^Mu_L1\\[\\d+\\]$", colnames(MCMC_matrix))]
MCMC_Sd_L1 <- MCMC_matrix[, grep("^Sd_L1\\[\\d+\\]$", colnames(MCMC_matrix))]

# Calculation
Lt <- array(NA, dim = c(1000, 24))
Pit <- array(NA, dim = c(1000, 24))
for (t in seq_len(24)) {
  mu_t <- quantile(MCMC_Mu_L1[, t], probs = 0.5)
  sd_t <- quantile(MCMC_Sd_L1[, t], probs = 0.5)
  
  Lt[, t] <- Cnf_l(mu = mu_t, sd = sd_t, nbclass = 1000, Lmin = lmin, Lmax = lmax)
  Pit[, t] <- Cnf_pi(mu = mu_t, sd = sd_t, nbclass = 1000, Lmin = lmin, Lmax = lmax)
}

# Dataframe
data <- data.frame(
  Cohort = rep(as.factor(1996:2019), each = 1000),
  Lt = as.vector(Lt[, 1:24]),
  Pit = as.vector(Pit[, 1:24])
)

truc <- data %>%
  group_by(Cohort) %>%
  summarize(weighted_mean = weighted.mean(Lt, Pit, na.rm = TRUE))

# Plot
fig4_A <- ggplot() + 
  geom_density_ridges(data = data, aes(x = Lt, y = Cohort, height = Pit, group = Cohort),
                      stat = "identity", fill = "green", color = "green", linewidth = 0.8, scale = 0.9, alpha = 0.3) +
  geom_vline(data = truc, aes(xintercept = mean(weighted_mean)), color = "black", linewidth = 0.8, linetype = "dashed") +
  geom_point(data = truc, aes(x = weighted_mean, y = Cohort), color = "darkgreen", size = 1.8) +
  labs(x = "Scale length (mm)", y = "Year of smolt migration") + 
  scale_x_continuous(expand = c(0, 0), breaks = c(0.4, 0.8, 1.2, 1.6, 2)) +
  scale_y_discrete(breaks = unique(data$Cohort)[c(TRUE, FALSE)]) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey95"),
    axis.title  = element_text(size = 9, face = "bold"),
    axis.text   = element_text(size = 8),
    legend.position = "none",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text  = element_text(size = 8),
    strip.placement = "outside",
    strip.text  = element_text(size = 9, face = "bold")
  ) +
  coord_flip(xlim = c(0.38, 2.02), ylim = c(0, 25), expand = FALSE)
fig4_A

ggsave(
  filename = paste0("results/figures/", projects, "/fig4_A.pdf"),
  plot = fig4_A,
  width = 1961/300,
  height = 1440/300,
  units = "in"
)

## -----------------------------------------------------------------------------
## 3. Maturation
## -----------------------------------------------------------------------------
# Parameters
lmin <- const$Min_L2
lmax <- const$Max_L2
MCMC_Mu_L2 <- MCMC_matrix[, grep("^Mu_L2\\[\\d+\\]$", colnames(MCMC_matrix))]
MCMC_Sd_L2 <- MCMC_matrix[, grep("^Sd_L2\\[\\d+\\]$", colnames(MCMC_matrix))]

# Calculation
Lt <- array(NA, dim = c(1000, 24))
Pit <- array(NA, dim = c(1000, 24))
for (t in seq_len(24)) {
  mu_t <- quantile(MCMC_Mu_L2[, t], probs = 0.5)
  sd_t <- quantile(MCMC_Sd_L2[, t], probs = 0.5)
  
  Lt[, t] <- Cnf_l(mu = mu_t, sd = sd_t, nbclass = 1000, Lmin = lmin, Lmax = lmax)
  Pit[, t] <- Cnf_pi(mu = mu_t, sd = sd_t, nbclass = 1000, Lmin = lmin, Lmax = lmax)
}

# Dataframe
data2 <- data.frame(
  Cohort = rep(as.factor(1996:2019), each = 1000),
  Lt = as.vector(Lt[, 1:24]),
  Pit = as.vector(Pit[, 1:24])
)

truc2 <- data2 %>%
  group_by(Cohort) %>%
  summarize(weighted_mean = weighted.mean(Lt, Pit, na.rm = TRUE))

# Plot
fig4_B <- ggplot() + 
  geom_density_ridges(data = data2, aes(x = Lt, y = Cohort, height = Pit, group = Cohort),
                      stat = "identity", fill = "purple", color = "purple", linewidth = 0.8, scale = 0.9, alpha = 0.3) +
  geom_vline(data = truc2, aes(xintercept = mean(weighted_mean)), color = "black", linewidth = 0.8, linetype = "dashed") +
  geom_point(data = truc2, aes(x = weighted_mean, y = Cohort), color = "purple4", size = 1.8) +
  labs(x = "Scale length (mm)", y = "Year of smolt migration") + 
  scale_x_continuous(expand = c(0, 0), breaks = c(1.8, 2.2, 2.6, 3, 3.4, 3.8)) +
  scale_y_discrete(breaks = unique(data$Cohort)[c(TRUE, FALSE)]) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey95"),
    axis.title  = element_text(size = 9, face = "bold"),
    axis.text   = element_text(size = 8),
    legend.position = "none",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text  = element_text(size = 8),
    strip.placement = "outside",
    strip.text  = element_text(size = 9, face = "bold")
  ) +
  coord_flip(xlim = c(1.78, 4.02), ylim = c(0, 25), expand = FALSE)
fig4_B

ggsave(
  filename = paste0("results/figures/", projects, "/fig4_B.pdf"),
  plot = fig4_B,
  width = 1961/300,
  height = 1440/300,
  units = "in"
)