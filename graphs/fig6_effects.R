### ============================================================================
### Disentangling length-dependent and length-independent variations
### in survival and maturation in Atlantic salmon
### ----------------------------------------------------------------------------
### File:    fig6_effects.R
### Purpose: Produce figure 6 Annual variation in population-level vital rates
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
pkginstall(c("dplyr", "tidyr", "ggplot2", "coda", "qs", "ggnewscale", "relaimpo"))

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
num_cohorts <- const$C
num_iterations <- nrow(MCMC_matrix)
num_class <- const$L

## -----------------------------------------------------------------------------
## 3. Survival
## -----------------------------------------------------------------------------
# Parameters
l1min <- const$Min_L1
l1max <- const$Max_L1
l1mean <- const$Mean_L1

MCMC_beta1 <- MCMC_matrix[, grep("^Beta1", colnames(MCMC_matrix))]
MCMC_alpha1 <- MCMC_matrix[, grep("^Alpha1", colnames(MCMC_matrix))]
MCMC_mualpha1 <- MCMC_matrix[, grep("^Mu_Alpha1", colnames(MCMC_matrix))]
MCMC_Mu_L1 <- MCMC_matrix[, grep("^Mu_L1\\[\\d+\\]$", colnames(MCMC_matrix))]
MCMC_Sd_L1 <- MCMC_matrix[, grep("^Sd_L1\\[\\d+\\]$", colnames(MCMC_matrix))]

# Calculation
L1_base <- array(NA, dim = c(num_class, num_cohorts, num_iterations))
Pi1_base <- array(NA, dim = c(num_class, num_cohorts, num_iterations))

surv_base <- array(NA, dim = c(num_class, num_cohorts,num_iterations))
surv_size <- array(NA, dim = c(num_class, num_cohorts,num_iterations))
surv_time <- array(NA, dim = c(num_class, num_cohorts,num_iterations))
surv_moy <- array(NA, dim = c(num_class, num_cohorts,num_iterations))

survM_base <- array(NA, dim = c(num_cohorts,num_iterations))
survM_size <- array(NA, dim = c(num_cohorts,num_iterations))
for (y in seq_len(num_cohorts)) {
  for (i in seq_len(num_iterations)) {
    L1_base[,y,i] <- Cnf_l(mu = MCMC_Mu_L1[i,y], sd = MCMC_Sd_L1[i,y], nbclass = num_class, Lmin = l1min, Lmax = l1max)
    Pi1_base[,y,i] <- Cnf_pi(mu = MCMC_Mu_L1[i,y], sd = MCMC_Sd_L1[i,y], nbclass = num_class, Lmin = l1min, Lmax = l1max)
    
    surv_base[,y,i] <- 1 / (1 + exp(- (MCMC_alpha1[i,y] + (MCMC_beta1[i] * (L1_base[,y,i] - l1mean)))))
    survM_base[y,i] <- weighted.mean(surv_base[,y,i], Pi1_base[,y,i])
    
    surv_size[,y,i] <- 1 / (1 + exp(- (MCMC_mualpha1[i] + (MCMC_beta1[i] * (L1_base[,y,i] - l1mean)))))
    survM_size[y,i] <- weighted.mean(surv_size[,y,i], Pi1_base[,y,i])
    
    surv_time[,y,i] <- 1 / (1 + exp(- (MCMC_alpha1[i,y])))
    surv_moy[,y,i] <- 1 / (1 + exp(- MCMC_mualpha1[i]))
    
  }
}
survQ_base <- apply(survM_base, 1, function(x) quantile(x, probs = c(0.25,0.5,0.75), na.rm = TRUE))
survQ_size <- apply(survM_size, 1, function(x) quantile(x, probs = 0.5, na.rm = TRUE))

survM_time <- apply(surv_time, c(2,3), function(x) mean(x, na.rm = TRUE))
survQ_time <- apply(survM_time, 1, function(x) quantile(x, probs = 0.5, na.rm = TRUE))

survM_moy <- apply(surv_moy, c(2,3), function(x) mean(x, na.rm = TRUE))
survQ_moy <- apply(survM_moy, 1, function(x) quantile(x, probs = 0.5, na.rm = TRUE))

# Dataframe
data1 <- data.frame(
  cohort = rep(1996:2019),
  
  survQ25_base = as.vector(survQ_base[1,]),
  survQ50_base = as.vector(survQ_base[2,]),
  survQ75_base = as.vector(survQ_base[3,]),
  
  survQ50_size = as.vector(survQ_size),
  
  survQ50_time = as.vector(survQ_time),
  
  survQ50_moy = as.vector(survQ_moy)
)

model <- lm(survQ50_base ~ survQ50_size + survQ50_time, data = data1)
summary(model)
adjr2_size <- round(calc.relimp(model, type = "lmg", rela = TRUE)$lmg[1], digits = 2)
adjr2_time <- round(calc.relimp(model, type = "lmg", rela = TRUE)$lmg[2], digits = 2)

data2 <- data.frame(
  cohort = rep(1996:2019),

  survQ50_size = as.vector(survQ_size),

  survQ50_time = as.vector(survQ_time)
)
data2_long <- data2 %>%
  pivot_longer(
    cols = starts_with("surv"),
    names_to = c("stat", "type"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = "stat",
    values_from = "value"
  ) %>%
  mutate(type = case_when(
    type == "time" ~ "Size-independent",
    type == "size" ~ "Size-dependent"
  ))

# Plot
fig6_A <- ggplot(data = data1, aes(x = cohort)) + 
  geom_linerange(data = data2_long, aes(ymin = mean(survQ_moy), ymax = survQ50, x = cohort, color = type), position = position_dodge(width = 0.7), linewidth = 3) +
  geom_hline(aes(yintercept = mean(survQ50_moy)), color = "black", linewidth = 0.8) +
  geom_ribbon(aes(ymin = survQ25_base, ymax = survQ75_base), fill = "green3", alpha = 0.3) +
  geom_line(aes(y = survQ50_base), color = "green3", linewidth = 0.8) + 
  geom_point(aes(y = survQ50_base), color = "green3", size = 1.8) + 
  coord_cartesian(xlim = c(1995.5,2019.5), ylim = c(-0.01,0.36), expand = FALSE) +
  scale_x_continuous(breaks = unique(data2_long$cohort)[c(TRUE, FALSE)]) +
  scale_y_continuous(breaks = c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35))+
  scale_color_manual(
    values = c("Size-independent" = "darkgreen", "Size-dependent" = "yellowgreen"),
    labels = c(
      paste0("Size-dependent\nR²lmg = ", adjr2_size),
      paste0("Size-independent\nR²lmg = ", adjr2_time)
    ), name = "Temporal variation") +
  labs(x = "Year of smolt migration", y = "Survival rate") + 
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey95"),
    axis.title  = element_text(size = 9, face = "bold"),
    axis.text   = element_text(size = 8),
    legend.position = "inside",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text  = element_text(size = 8),
    legend.position.inside = c(0.85,0.85),
    strip.placement = "outside",
    strip.text  = element_text(size = 9, face = "bold")
  )
fig6_A

ggsave(
  filename = "results/fig6_A.pdf",
  plot = fig6_A,
  width = 1961/300,
  height = 1440/300,
  units = "in"
)


## -----------------------------------------------------------------------------
## 4. Maturation
## -----------------------------------------------------------------------------
# Parameters
l2min <- const$Min_L2
l2max <- const$Max_L2
l2mean <- const$Mean_L2

MCMC_beta2 <- MCMC_matrix[, grep("^Beta2", colnames(MCMC_matrix))]
MCMC_alpha2F <- MCMC_matrix[, grep("^Alpha2\\[\\d+, 1\\]$", colnames(MCMC_matrix))]
MCMC_mualpha2F <- MCMC_matrix[, grep("^Mu_Alpha2\\[1\\]$", colnames(MCMC_matrix))]
MCMC_alpha2M <- MCMC_matrix[, grep("^Alpha2\\[\\d+, 2\\]$", colnames(MCMC_matrix))]
MCMC_mualpha2M <- MCMC_matrix[, grep("^Mu_Alpha2\\[2\\]$", colnames(MCMC_matrix))]
MCMC_Mu_L2 <- MCMC_matrix[, grep("^Mu_L2\\[\\d+\\]$", colnames(MCMC_matrix))]
MCMC_Sd_L2 <- MCMC_matrix[, grep("^Sd_L2\\[\\d+\\]$", colnames(MCMC_matrix))]

# Calculation
L2_base <- array(NA, dim = c(num_class, num_cohorts, num_iterations))
Pi2_base <- array(NA, dim = c(num_class, num_cohorts, num_iterations))

mat_baseF <- array(NA, dim = c(num_class, num_cohorts,num_iterations))
mat_sizeF <- array(NA, dim = c(num_class, num_cohorts,num_iterations))
mat_timeF <- array(NA, dim = c(num_class, num_cohorts,num_iterations))
mat_moyF <- array(NA, dim = c(num_class, num_cohorts,num_iterations))
mat_baseM <- array(NA, dim = c(num_class, num_cohorts,num_iterations))
mat_sizeM <- array(NA, dim = c(num_class, num_cohorts,num_iterations))
mat_timeM <- array(NA, dim = c(num_class, num_cohorts,num_iterations))
mat_moyM <- array(NA, dim = c(num_class, num_cohorts,num_iterations))

matM_baseF <- array(NA, dim = c(num_cohorts,num_iterations))
matM_sizeF <- array(NA, dim = c(num_cohorts,num_iterations))
matM_baseM <- array(NA, dim = c(num_cohorts,num_iterations))
matM_sizeM <- array(NA, dim = c(num_cohorts,num_iterations))
for (y in seq_len(num_cohorts)) {
  for (i in seq_len(num_iterations)) {
    L2_base[,y,i] <- Cnf_l(mu = MCMC_Mu_L2[i,y], sd = MCMC_Sd_L2[i,y], nbclass = num_class, Lmin = l2min, Lmax = l2max)
    Pi2_base[,y,i] <- Cnf_pi(mu = MCMC_Mu_L2[i,y], sd = MCMC_Sd_L2[i,y], nbclass = num_class, Lmin = l2min, Lmax = l2max)
    
    mat_baseF[,y,i] <- 1 / (1 + exp(- (MCMC_alpha2F[i,y] + (MCMC_beta2[i] * (L2_base[,y,i] - l2mean)))))
    matM_baseF[y,i] <- weighted.mean(mat_baseF[,y,i], Pi2_base[,y,i])
    mat_baseM[,y,i] <- 1 / (1 + exp(- (MCMC_alpha2M[i,y] + (MCMC_beta2[i] * (L2_base[,y,i] - l2mean)))))
    matM_baseM[y,i] <- weighted.mean(mat_baseM[,y,i], Pi2_base[,y,i])
    
    mat_sizeF[,y,i] <- 1 / (1 + exp(- (MCMC_mualpha2F[i] + (MCMC_beta2[i] * (L2_base[,y,i] - l2mean)))))
    matM_sizeF[y,i] <- weighted.mean(mat_sizeF[,y,i], Pi2_base[,y,i])
    mat_sizeM[,y,i] <- 1 / (1 + exp(- (MCMC_mualpha2M[i] + (MCMC_beta2[i] * (L2_base[,y,i] - l2mean)))))
    matM_sizeM[y,i] <- weighted.mean(mat_sizeM[,y,i], Pi2_base[,y,i])
    
    mat_timeF[,y,i] <- 1 / (1 + exp(- (MCMC_alpha2F[i,y])))
    mat_timeM[,y,i] <- 1 / (1 + exp(- (MCMC_alpha2M[i,y])))
    
    mat_moyF[,y,i] <- 1 / (1 + exp(- MCMC_mualpha2F[i]))
    mat_moyM[,y,i] <- 1 / (1 + exp(- MCMC_mualpha2M[i]))
  }
}
matQ_baseF <- apply(matM_baseF, 1, function(x) quantile(x, probs = c(0.25,0.5,0.75), na.rm = TRUE))
matQ_baseM <- apply(matM_baseM, 1, function(x) quantile(x, probs = c(0.25,0.5,0.75), na.rm = TRUE))

matQ_sizeF <- apply(matM_sizeF, 1, function(x) quantile(x, probs = 0.5, na.rm = TRUE))
matQ_sizeM <- apply(matM_sizeM, 1, function(x) quantile(x, probs = 0.5, na.rm = TRUE))

matM_timeF <- apply(mat_timeF, c(2,3), function(x) mean(x, na.rm = TRUE))
matQ_timeF <- apply(matM_timeF, 1, function(x) quantile(x, probs = 0.5, na.rm = TRUE))
matM_timeM <- apply(mat_timeM, c(2,3), function(x) mean(x, na.rm = TRUE))
matQ_timeM <- apply(matM_timeM, 1, function(x) quantile(x, probs = 0.5, na.rm = TRUE))

matM_moyF <- apply(mat_moyF, c(2,3), function(x) mean(x, na.rm = TRUE))
matQ_moyF <- apply(matM_moyF, 1, function(x) quantile(x, probs = 0.5, na.rm = TRUE))
matM_moyM <- apply(mat_moyM, c(2,3), function(x) mean(x, na.rm = TRUE))
matQ_moyM <- apply(matM_moyM, 1, function(x) quantile(x, probs = 0.5, na.rm = TRUE))

# Dataframe
data3 <- data.frame(
  cohort = rep(1996:2019),
  
  matQ25_baseF = as.vector(matQ_baseF[1,]),
  matQ50_baseF = as.vector(matQ_baseF[2,]),
  matQ75_baseF = as.vector(matQ_baseF[3,]),
  matQ50_sizeF = as.vector(matQ_sizeF),
  matQ50_timeF = as.vector(matQ_timeF),
  
  matQ25_baseM = as.vector(matQ_baseM[1,]),
  matQ50_baseM = as.vector(matQ_baseM[2,]),
  matQ75_baseM = as.vector(matQ_baseM[3,]),
  matQ50_sizeM = as.vector(matQ_sizeM),
  matQ50_timeM = as.vector(matQ_timeM),
  
  matQ50_moyF = as.vector(matQ_moyF),
  matQ50_moyM = as.vector(matQ_moyM)
)

modelF <- lm(matQ50_baseF ~ matQ50_sizeF + matQ50_timeF, data = data3)
summary(modelF)
adjr2_sizeF <- round(calc.relimp(modelF, type = "lmg", rela = TRUE)$lmg[1], digits = 2)
adjr2_timeF <- round(calc.relimp(modelF, type = "lmg", rela = TRUE)$lmg[2], digits = 2)

modelM <- lm(matQ50_baseM ~ matQ50_sizeM + matQ50_timeM, data = data3)
adjr2_sizeM <- round(calc.relimp(modelM, type = "lmg", rela = TRUE)$lmg[1], digits = 2)
adjr2_timeM <- round(calc.relimp(modelM, type = "lmg", rela = TRUE)$lmg[2], digits = 2)

data4 <- data.frame(
  cohort = rep(1996:2019),
  
  matQ50_sizeF = as.vector(matQ_sizeF),
  
  matQ50_timeF = as.vector(matQ_timeF)
)
data4_long <- data4 %>%
  pivot_longer(
    cols = starts_with("mat"),
    names_to = c("stat", "type"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = "stat",
    values_from = "value"
  ) %>%
  mutate(type = case_when(
    type == "timeF" ~ "Size-independent",
    type == "sizeF" ~ "Size-dependent"
  ))

data5 <- data.frame(
  cohort = rep(1996:2019),
  
  matQ50_sizeM = as.vector(matQ_sizeM),
  
  matQ50_timeM = as.vector(matQ_timeM)
)
data5_long <- data5 %>%
  pivot_longer(
    cols = starts_with("mat"),
    names_to = c("stat", "type"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = "stat",
    values_from = "value"
  ) %>%
  mutate(type = case_when(
    type == "timeM" ~ "Size-independent",
    type == "sizeM" ~ "Size-dependent"
  ))

# Plot
fig6_B <- ggplot(data = data3, aes(x = cohort)) + 
  geom_linerange(data = data4_long, aes(ymin=mean(matQ_moyF), ymax=matQ50, x=cohort, color=type), position = position_dodge(width = 0.7), linewidth = 3) +
  geom_hline(aes(yintercept = mean(matQ50_moyF)), color = "black", linewidth = 0.8) +
  geom_ribbon(aes(ymin = matQ25_baseF, ymax = matQ75_baseF), fill = "red", alpha = 0.3) +
  geom_line(aes(y = matQ50_baseF), color = "red", linewidth = 0.8) + 
  geom_point(aes(y = matQ50_baseF), color = "red", size = 1.8) + 
  scale_color_manual(
    values = c("Size-independent" = "darkred", "Size-dependent" = "coral"),
    labels = c(
      paste0("Size-dependent\nR²lmg = ", adjr2_sizeF),
      paste0("Size-independent\nR²lmg = ", adjr2_timeF)
    ),
    name = "Temporal variation (F)",
    guide = guide_legend(order = 2)) +
  new_scale_color() +
  geom_linerange(data = data5_long, aes(ymin=mean(matQ_moyM), ymax=matQ50, x=cohort, color=type), position = position_dodge(width = 0.7), linewidth = 3) +
  geom_hline(aes(yintercept = mean(matQ50_moyM)), color = "black", linewidth = 0.8) +
  geom_ribbon(aes(ymin = matQ25_baseM, ymax = matQ75_baseM), fill = "blue", alpha = 0.3) +
  geom_line(aes(y = matQ50_baseM), color = "blue", linewidth = 0.8) + 
  geom_point(aes(y = matQ50_baseM), color = "blue", size = 1.8) + 
  scale_color_manual(
    values = c("Size-independent" = "darkblue", "Size-dependent" = "deepskyblue"),
    labels = c(
      paste0("Size-dependent\nR²lmg = ", adjr2_sizeM),
      paste0("Size-independent\nR²lmg = ", adjr2_timeM)
    ),
    name = "Temporal variation (M)",
    guide = guide_legend(order = 1)) +
  coord_cartesian(xlim = c(1995.5,2019.5), ylim = c(0.28,1.02), expand = FALSE) +
  scale_x_continuous(breaks = unique(data2_long$cohort)[c(TRUE, FALSE)]) +
  scale_y_continuous(breaks = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))+
  labs(x = "Year of smolt migration", y = "Maturation rate") + 
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey95"),
    axis.title  = element_text(size = 9, face = "bold"),
    axis.text   = element_text(size = 8),
    legend.position = "inside",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text  = element_text(size = 8),
    legend.position.inside = c(0.33,0.18),
    legend.box = "horizontal",
    strip.placement = "outside",
    strip.text  = element_text(size = 9, face = "bold")
  )
fig6_B

ggsave(
  filename = "results/fig6_B.pdf",
  plot = fig6_B,
  width = 1961/300,
  height = 1440/300,
  units = "in"
)
