### ============================================================================
### Disentangling length-dependent and length-independent variations
### in survival and maturation in Atlantic salmon
### ----------------------------------------------------------------------------
### File:    figA7_cor.R
### Purpose: Builds Appendix 7.1 Figures (A7.1.1 & A7.1.2).
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
pkginstall(c("ggplot2", "dplyr", "tidyr", "coda", "nimble", "qs"))

## -----------------------------------------------------------------------------
## 1. Define project settings
## -----------------------------------------------------------------------------
projects <- c("M1")

## -----------------------------------------------------------------------------
## 2. Import MCMC object and data
## -----------------------------------------------------------------------------
load_mcmc <- function(project) {
  MCMC <- qread(paste0("saves/", project, "/MCMC.qs"))
  MCMC %>%
    seq_along() %>%
    lapply(function(i) mcmc(MCMC[[i]]$samples)) %>%
    mcmc.list()
}
MCMC_samples <- load_mcmc(projects)
MCMC_matrix <- as.matrix(MCMC_samples)

source("functions/nimblefunctions/nf_l.R")
Cnf_l <- compileNimble(nf_l)
source("functions/nimblefunctions/nf_pi.R")
Cnf_pi <- compileNimble(nf_pi)

const <- qread(file = "data/realdata/const.qs")
num_cohorts <- const$C
num_iterations <- nrow(MCMC_matrix)
num_class <- const$L

## -----------------------------------------------------------------------------
## 1. Alpha values
## -----------------------------------------------------------------------------
# Correlation test
MCMC_alpha2F <- MCMC_matrix[, grep("^Alpha2", colnames(MCMC_matrix))][,1:24]
MCMC_alpha2M <- MCMC_matrix[, grep("^Alpha2", colnames(MCMC_matrix))][,25:48]

cor_test1 <- mapply(
  function(x, y) cor.test(x, y, method = "pearson"),
  split(MCMC_alpha2F, row(MCMC_alpha2F)),
  split(MCMC_alpha2M, row(MCMC_alpha2M)),
  SIMPLIFY = FALSE
)
cor_values1 <- sapply(cor_test1, function(t) t$estimate)
data1 <- data.frame(
  samples = 1:30000,
  pearson = cor_values1
)
p_neg1 <- mean(data1$pearson <= 0)

# Plot
FigA7.1.1 <- ggplot(data1, aes(x = pearson)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 100,
                 fill = "white",
                 color = "black") +
  geom_histogram(data = subset(data1, pearson <= 0),
                 aes(y = after_stat(density * p_neg1)),
                 bins = 100,
                 fill = "red",
                 alpha = 0.4) +
  labs(x = "Pearson correlation (r)", y = "Density") + 
  geom_vline(xintercept = 0, color = "red", linewidth = 1.2, linetype = "dashed") +
  geom_vline(xintercept = quantile(data1$pearson, probs = c(0.5)), color = "black", linewidth = 1.2) +
  scale_x_continuous(breaks = c(-1,-0.8,-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1)) +
  coord_cartesian(xlim = c(-1,1), ylim = c(0,3)) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 35, face = "bold"),
    axis.text.x = element_text(size = 25,face = "italic"),
    axis.text.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey75"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )
FigA7.1.1

ggsave(
  filename = paste0("results/figures/", projects, "/FigA7.1.1.pdf"),
  plot = FigA7.1.1,
  width = 1961/300,
  height = 1440/300,
  units = "in"
)

## -----------------------------------------------------------------------------
## 3. Maturation rates
## -----------------------------------------------------------------------------
# Correlation test
l2mean <- const$Mean_L2
l2min <- 1.5
l2max <- 4
MCMC_Mu_L2 <- MCMC_matrix[, grep("^Mu_L2\\[\\d+\\]$", colnames(MCMC_matrix))]
MCMC_Sd_L2 <- MCMC_matrix[, grep("^Sd_L2\\[\\d+\\]$", colnames(MCMC_matrix))]
MCMC_beta2 <- MCMC_matrix[, grep("^Beta2", colnames(MCMC_matrix))]

L2_base <- array(NA, dim = c(num_class, num_cohorts, num_iterations))
Pi2_base <- array(NA, dim = c(num_class, num_cohorts, num_iterations))
mat_baseF <- array(NA, dim = c(num_class, num_cohorts,num_iterations))
mat_baseM <- array(NA, dim = c(num_class, num_cohorts,num_iterations))
matM_baseF <- array(NA, dim = c(num_cohorts,num_iterations))
matM_baseM <- array(NA, dim = c(num_cohorts,num_iterations))
for (y in seq_len(num_cohorts)) {
  for (i in seq_len(num_iterations)) {
    L2_base[,y,i] <- Cnf_l(mu = MCMC_Mu_L2[i,y], sd = MCMC_Sd_L2[i,y], nbclass = num_class, Lmin = l2min, Lmax = l2max)
    Pi2_base[,y,i] <- Cnf_pi(mu = MCMC_Mu_L2[i,y], sd = MCMC_Sd_L2[i,y], nbclass = num_class, Lmin = l2min, Lmax = l2max)
    
    mat_baseF[,y,i] <- 1 / (1 + exp(- (MCMC_alpha2F[i,y] + (MCMC_beta2[i] * (L2_base[,y,i] - l2mean)))))
    matM_baseF[y,i] <- weighted.mean(mat_baseF[,y,i], Pi2_base[,y,i])
    mat_baseM[,y,i] <- 1 / (1 + exp(- (MCMC_alpha2M[i,y] + (MCMC_beta2[i] * (L2_base[,y,i] - l2mean)))))
    matM_baseM[y,i] <- weighted.mean(mat_baseM[,y,i], Pi2_base[,y,i])
  }
}

cor_test2 <- mapply(
  function(x, y) cor.test(x, y, method = "pearson"),
  split(t(matM_baseF), row(t(matM_baseF))),
  split(t(matM_baseM), row(t(matM_baseM))),
  SIMPLIFY = FALSE
)
cor_values2 <- sapply(cor_test2, function(t) t$estimate)
data2 <- data.frame(
  samples = 1:30000,
  pearson = cor_values2
)
p_neg2 <- mean(data2$pearson <= 0)

FigA7.1.2 <- ggplot(data2, aes(x = pearson)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 100,
                 fill = "white",
                 color = "black") +
  geom_histogram(data = subset(data2, pearson <= 0),
                 aes(y = after_stat(density * p_neg2)),
                 bins = 100,
                 fill = "red",
                 alpha = 0.4) +
  labs(x = "Pearson correlation (r)", y = "Density") + 
  geom_vline(xintercept = 0, color = "red", linewidth = 1.2, linetype = "dashed") +
  geom_vline(xintercept = quantile(data2$pearson, probs = c(0.5)), color = "black", linewidth = 1.2) +
  scale_x_continuous(breaks = c(-1,-0.8,-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1)) +
  coord_cartesian(xlim = c(-1,1), ylim = c(0,3)) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 35, face = "bold"),
    axis.text.x = element_text(size = 25,face = "italic"),
    axis.text.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey75"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )
FigA7.1.2

ggsave(
  filename = paste0("results/figures/", projects, "/FigA7.1.2.pdf"),
  plot = FigA7.1.2,
  width = 1961/300,
  height = 1440/300,
  units = "in"
)