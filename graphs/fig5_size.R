### ============================================================================
### Disentangling length-dependent and length-independent variations
### in survival and maturation in Atlantic salmon
### ----------------------------------------------------------------------------
### File:    fig5_size.R
### Purpose: Produce figure 5 Annual length-dependent vital rates
### Author:  ©BOULAIRE Eliot, NEVOUX Marie & RIVOT Etienne
### Version: 17/08/2026
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
pkginstall(c("dplyr", "tidyr", "ggplot2", "coda", "qs", "ggnewscale"))

## -----------------------------------------------------------------------------
## 1. Define project settings
## -----------------------------------------------------------------------------
projects <- "M1"

## -----------------------------------------------------------------------------
## 2. Import MCMC object and functions
## -----------------------------------------------------------------------------
source("functions/nf_l.R")
Cnf_l <- compileNimble(nf_l)
source("functions/nf_pi.R")
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

data <- qread(file = "data/realdata/data.qs")
const <- qread(file = "data/realdata/const.qs")

## -----------------------------------------------------------------------------
## 3. Survival
## -----------------------------------------------------------------------------
# Parameters
lmin <- const$Min_L1
lmax <- const$Max_L1
lmean <- const$Mean_L1
num_lseq <- 1000
lseq <- seq(from = lmin, to = lmax, length.out = num_lseq)

MCMC_beta1 <- MCMC_matrix[, grep("^Beta1", colnames(MCMC_matrix))]
MCMC_alpha1 <- MCMC_matrix[, grep("^Alpha1", colnames(MCMC_matrix))]
num_years <- ncol(MCMC_alpha1)

S1 <- as.vector(na.omit(data$S1))
mu_S1 <- mean(S1)

quantiles <- c(0.25, 0.50, 0.75)
num_quantiles <- length(quantiles)

# Calculation
surv <- array(data = NA, dim = c(num_lseq, num_years, num_quantiles))
for (l in 1:num_lseq) {
  for (t in 1:num_years) {
    surv_temp <- (1 / (1 + exp(-(MCMC_alpha1[,t] + (MCMC_beta1 * (lseq[l] - lmean))))))
    surv[l, t, ] <- quantile(surv_temp, probs = quantiles, na.rm = TRUE)
  }
}

# Dataframe
df_surv <- data.frame(
  Cohort = as.character(rep(c(1996:2019), each = num_lseq)),
  Taillecailles = rep(lseq, times = num_years),
  Q25 = as.vector(surv[,,1]),
  Q50 = as.vector(surv[,,2]),
  Q75 = as.vector(surv[,,3])
)
df_surv2 <- df_surv %>%
  filter(Cohort %in% c("2002","2008"))

# Plot
fig5_A <- ggplot() +
  geom_line(data = df_surv, aes(x = Taillecailles, y = Q50, group = Cohort), color = "darkgreen", linewidth = 0.4) +
  geom_ribbon(data = df_surv2, aes(x = Taillecailles, ymin = Q25, ymax = Q75, group = Cohort,  fill = Cohort), alpha = 0.3) +
  geom_line(data = df_surv2, aes(x = Taillecailles, y = Q50, group = Cohort, color = Cohort), linewidth = 0.8) +
  geom_point(aes(x = S1, y = -0.05), color = "grey", size = 0.9, position = position_jitter(width = 0, height = 0.04)) +
  geom_vline(aes(xintercept = mu_S1), color = "black", linewidth = 0.8, linetype = "dashed") +
  geom_hline(aes(yintercept = 0), color = "black", linewidth = 0.8) +
  labs(x="Scale length (mm)", y="Survival rate") +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0.4, 0.8, 1.2, 1.6, 2)) +
  coord_cartesian(ylim = c(-0.1, 0.82), xlim = c(0.38, 2.02), expand = FALSE) +
  scale_color_manual(values = c("2002" = "yellowgreen", "2008" = "green"),
                     name = "Year of smolt migration") +
  scale_fill_manual(values = c("2002" = "yellowgreen", "2008" = "green"),
                    name = "Year of smolt migration") +
theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey95"),
    axis.title  = element_text(size = 9, face = "bold"),
    axis.text   = element_text(size = 8),
    legend.position = "inside",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text  = element_text(size = 8),
    legend.position.inside = c(0.19,0.65),
    strip.placement = "outside",
    strip.text  = element_text(size = 9, face = "bold")
  )
fig5_A

ggsave(
  filename = "results/fig5_A.pdf",
  plot = fig5_A,
  width = 1961/300,
  height = 1440/300,
  units = "in"
)

## -----------------------------------------------------------------------------
## 3. Maturation
## -----------------------------------------------------------------------------
# Parameters
lmin2 <- const$Min_L2
lmax2 <- const$Max_L2
lmean2 <- const$Mean_L2
num_lseq <- 1000
lseq2 <- seq(from = lmin2, to = lmax2, length.out = num_lseq)

MCMC_beta2 <- MCMC_matrix[, grep("^Beta2", colnames(MCMC_matrix))]
MCMC_alpha2F <- MCMC_matrix[, grep("^Alpha2\\[\\d+, 1\\]$", colnames(MCMC_matrix))]
MCMC_alpha2M <- MCMC_matrix[, grep("^Alpha2\\[\\d+, 2\\]$", colnames(MCMC_matrix))]
num_years <- ncol(MCMC_alpha2F)

S2 <- as.vector(na.omit(data$S2))
mu_S2 <- mean(S2)

quantiles <- c(0.25, 0.50, 0.75)
num_quantiles <- length(quantiles)

# Calculation
matF <- array(data = NA, dim = c(num_lseq, num_years, num_quantiles))
matM <- array(data = NA, dim = c(num_lseq, num_years, num_quantiles))
for (l in 1:num_lseq) {
  for (t in 1:num_years) {
    matF_temp <- (1 / (1 + exp(-(MCMC_alpha2F[,t] + (MCMC_beta2 * (lseq2[l] - lmean2))))))
    matF[l, t, ] <- quantile(matF_temp, probs = quantiles, na.rm = TRUE)
    
    matM_temp <- (1 / (1 + exp(-(MCMC_alpha2M[,t] + (MCMC_beta2 * (lseq2[l] - lmean2))))))
    matM[l, t, ] <- quantile(matM_temp, probs = quantiles, na.rm = TRUE)
  }
}

# Dataframe
# Females
df_matF <- data.frame(
  Cohort = as.character(rep(c(1996:2019), each = num_lseq)),
  Taillecailles = rep(lseq2, times = num_years),
  Q25 = as.vector(matF[,,1]),
  Q50 = as.vector(matF[,,2]),
  Q75 = as.vector(matF[,,3])
)
df_matF2 <- df_matF %>%
  filter(Cohort %in% c("2004","2014"))

# Males
df_matM <- data.frame(
  Cohort = as.character(rep(c(1996:2019), each = num_lseq)),
  Taillecailles = rep(lseq2, times = num_years),
  Q25 = as.vector(matM[,,1]),
  Q50 = as.vector(matM[,,2]),
  Q75 = as.vector(matM[,,3])
)
df_matM2 <- df_matM %>%
  filter(Cohort %in% c("2000","2007"))

# Plot
fig5_B <- ggplot() +
  geom_line(data = df_matM, aes(x = Taillecailles, y = Q50, group = Cohort), color = "darkblue", linewidth = 0.4) +
  geom_ribbon(data = df_matM2, aes(x = Taillecailles, ymin = Q25, ymax = Q75, group = Cohort, fill = Cohort), alpha = 0.3) +
  geom_line(data = df_matM2, aes(x = Taillecailles, y = Q50, group = Cohort, color = Cohort), linewidth = 0.8) +
  scale_color_manual(values = c("2000" = "deepskyblue", "2007" = "blue"),
                     name = "Year of smolt migration (M)",
                     guide = guide_legend(order = 1)) +
  scale_fill_manual(values = c("2000" = "deepskyblue", "2007" = "blue"),
                    name = "Year of smolt migration (M)",
                    guide = guide_legend(order = 1)) +
  new_scale_color() +
  new_scale_fill() +
  geom_line(data = df_matF, aes(x = Taillecailles, y = Q50, group = Cohort), color = "darkred", linewidth = 0.4) +
  geom_ribbon(data = df_matF2, aes(x = Taillecailles, ymin = Q25, ymax = Q75, group = Cohort, fill = Cohort), alpha = 0.3) +
  geom_line(data = df_matF2, aes(x = Taillecailles, y = Q50, group = Cohort, color = Cohort), linewidth = 0.8) +
  scale_color_manual(values = c("2004" = "red", "2014" = "coral"),
                     name = "Year of smolt migration (F)",
                     guide = guide_legend(order = 2)) +
  scale_fill_manual(values = c("2004" = "red", "2014" = "coral"),
                    name = "Year of smolt migration (F)",
                    guide = guide_legend(order = 2)) +
  geom_point(aes(x = S2, y = -0.05), color = "grey", size = 0.9, position = position_jitter(width = 0, height = 0.04)) +
  geom_vline(aes(xintercept = mu_S2), color = "black", linewidth = 0.8, linetype = "dashed") +
  geom_hline(aes(yintercept = 0), color = "black", linewidth = 0.8) +
  labs(x="Scale length (mm)", y="Maturation rate") +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(1.8, 2.2, 2.6, 3, 3.4, 3.8)) +
  coord_cartesian(ylim = c(-0.1, 1.02), xlim = c(1.78, 4.02)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey95"),
    axis.title  = element_text(size = 9, face = "bold"),
    axis.text   = element_text(size = 8),
    legend.position = "inside",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text  = element_text(size = 8),
    legend.position.inside = c(0.8,0.4),
    strip.placement = "outside",
    strip.text  = element_text(size = 9, face = "bold")
  )
fig5_B

ggsave(
  filename = "results/fig5_B.pdf",
  plot = fig5_B,
  width = 1961/300,
  height = 1440/300,
  units = "in"
)
