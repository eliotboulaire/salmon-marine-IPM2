### ============================================================================
### Disentangling length-dependent and length-independent variations
### in survival and maturation in Atlantic salmon
### ----------------------------------------------------------------------------
### File:    figA7_scatter.R
### Purpose: Builds Appendix 7.2 Figures (A7.2.1 & A7.2.2).
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
pkginstall(c("ggplot2", "dplyr", "tidyr", "coda", "nimble", "qs", "patchwork"))

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
## 3. Survival
## -----------------------------------------------------------------------------
# Parameters
l1mean <- const$Mean_L1
l1min <- 0
l1max <- 2.5

MCMC_beta1 <- MCMC_matrix[, grep("^Beta1", colnames(MCMC_matrix))]
MCMC_alpha1 <- MCMC_matrix[, grep("^Alpha1", colnames(MCMC_matrix))]
MCMC_mualpha1 <- MCMC_matrix[, grep("^Mu_Alpha1", colnames(MCMC_matrix))]
MCMC_Mu_L1 <- MCMC_matrix[, grep("^Mu_L1\\[\\d+\\]$", colnames(MCMC_matrix))]
MCMC_Sd_L1 <- MCMC_matrix[, grep("^Sd_L1\\[\\d+\\]$", colnames(MCMC_matrix))]

# Calcultion
L1_base <- array(NA, dim = c(num_class, num_cohorts, num_iterations))
Pi1_base <- array(NA, dim = c(num_class, num_cohorts, num_iterations))

surv_size <- array(NA, dim = c(num_class, num_cohorts,num_iterations))
survM_size <- array(NA, dim = c(num_cohorts,num_iterations))
for (y in seq_len(num_cohorts)) {
  for (i in seq_len(num_iterations)) {
    L1_base[,y,i] <- Cnf_l(mu = MCMC_Mu_L1[i,y], sd = MCMC_Sd_L1[i,y], nbclass = num_class, Lmin = l1min, Lmax = l1max)
    Pi1_base[,y,i] <- Cnf_pi(mu = MCMC_Mu_L1[i,y], sd = MCMC_Sd_L1[i,y], nbclass = num_class, Lmin = l1min, Lmax = l1max)
    
    surv_size[,y,i] <- 1 / (1 + exp(- (MCMC_mualpha1[i] + (MCMC_beta1[i] * (L1_base[,y,i] - l1mean)))))
    survM_size[y,i] <- weighted.mean(surv_size[,y,i], Pi1_base[,y,i])
  }
}
survQ50_size <- apply(survM_size, 1, function(x) quantile(x, probs = c(0.5), na.rm = TRUE))
L1mean_size <- colMeans(MCMC_Mu_L1)

# Dataframe
data1 <- data.frame(
  cohort = 1996:2019,
  size = L1mean_size,
  surv = survQ50_size
)

r_size <- range(data1$size, na.rm = TRUE)
r_surv <- range(data1$surv, na.rm = TRUE)
offset <- 0.15 * diff(r_size)
to_size <- function(x) (x - r_surv[1]) / diff(r_surv) * diff(r_size) + r_size[1] + offset
to_surv <- function(y) (y - offset - r_size[1]) / diff(r_size) * diff(r_surv) + r_surv[1]
yrs <- range(data1$cohort, na.rm = TRUE)

# Plot
FigA7.2.1 <- ggplot(data1, aes(x = cohort)) +
  geom_line(aes(y = size), colour = "darkgreen", linewidth = 1.5) +
  geom_point(aes(y = size), colour = "darkgreen", size = 5) +
  geom_line(aes(y = to_size(surv)), colour = "green3", linewidth = 1.5) +
  geom_point(aes(y = to_size(surv)), colour = "green3", size = 5) +
  scale_y_continuous(
    name     = "Mean smolt scale length (mm)",
    sec.axis = sec_axis(~ to_surv(.), name = "Length-dependent survival")
  ) +
  scale_x_continuous(
    name         = "Year of smolt migration",
    breaks       = seq(yrs[1], yrs[2], by = 2)
  ) +
  theme_minimal() +
  theme(
    axis.title        = element_text(size = 30, face = "bold"),
    axis.text         = element_text(size = 20, face = "bold"),
    axis.title.y.left  = element_text(colour = "darkgreen"),
    axis.text.y.left   = element_text(colour = "darkgreen"),
    axis.title.y.right = element_text(colour = "green3"),
    axis.text.y.right  = element_text(colour = "green3"),
    panel.grid.major.x = element_line(linewidth = 0.6),
    panel.grid.minor.x = element_line(linewidth = 0.3)
  )
FigA7.2.1

ggsave(
  filename = paste0("results/figures/", projects, "/FigA7.2.1.pdf"),
  plot = FigA7.2.1,
  width = 1961/300,
  height = 1440/300,
  units = "in"
)

## -----------------------------------------------------------------------------
## 4. Maturation
## -----------------------------------------------------------------------------
# Parameters
l2mean <- const$Mean_L2
l2min <- 1.5
l2max <- 4
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

mat_sizeF <- array(NA, dim = c(num_class, num_cohorts,num_iterations))
matM_sizeF <- array(NA, dim = c(num_cohorts,num_iterations))
mat_sizeM <- array(NA, dim = c(num_class, num_cohorts,num_iterations))
matM_sizeM <- array(NA, dim = c(num_cohorts,num_iterations))
for (y in seq_len(num_cohorts)) {
  for (i in seq_len(num_iterations)) {
    L2_base[,y,i] <- Cnf_l(mu = MCMC_Mu_L2[i,y], sd = MCMC_Sd_L2[i,y], nbclass = num_class, Lmin = l2min, Lmax = l2max)
    Pi2_base[,y,i] <- Cnf_pi(mu = MCMC_Mu_L2[i,y], sd = MCMC_Sd_L2[i,y], nbclass = num_class, Lmin = l2min, Lmax = l2max)
    
    mat_sizeF[,y,i] <- 1 / (1 + exp(- (MCMC_mualpha2F[i] + (MCMC_beta2[i] * (L2_base[,y,i] - l2mean)))))
    matM_sizeF[y,i] <- weighted.mean(mat_sizeF[,y,i], Pi2_base[,y,i])
    mat_sizeM[,y,i] <- 1 / (1 + exp(- (MCMC_mualpha2M[i] + (MCMC_beta2[i] * (L2_base[,y,i] - l2mean)))))
    matM_sizeM[y,i] <- weighted.mean(mat_sizeM[,y,i], Pi2_base[,y,i])
  }
}
matQ50_sizeF <- apply(matM_sizeF, 1, function(x) quantile(x, probs = c(0.5), na.rm = TRUE))
matQ50_sizeM <- apply(matM_sizeM, 1, function(x) quantile(x, probs = c(0.5), na.rm = TRUE))
L2mean_size <- colMeans(MCMC_Mu_L2)

# Dataframe
data2 <- data.frame(
  cohort = 1996:2019,
  size = L2mean_size,
  matF = matQ50_sizeF,
  matM = matQ50_sizeM
)
yrs <- range(data2$cohort, na.rm = TRUE)

# Plot females
r_size <- range(data2$size, na.rm = TRUE)
r_mat  <- range(data2$matF, na.rm = TRUE)
offset <- 0.15 * diff(r_size)
  
to_size <- function(x) (x - r_mat[1]) / diff(r_mat) * diff(r_size) + r_size[1] + offset
to_mat  <- function(y) (y - offset - r_size[1]) / diff(r_size) * diff(r_mat) + r_mat[1]
  
pF <- ggplot(data2, aes(x = cohort)) +
  geom_line(aes(y = size), colour = "purple", linewidth = 1.5) +
  geom_point(aes(y = size), colour = "purple", size = 5) +
  geom_line(aes(y = to_size(matF)), colour = "red", linewidth = 1.5) +
  geom_point(aes(y = to_size(matF)), colour = "red", size = 5) +
  scale_y_continuous(
    name     = "Mean smolt\nscale length (mm)",
    sec.axis = sec_axis(~ to_mat(.), name = "Length-dependent\nmaturation")
  ) +
  scale_x_continuous(
    name         = "Year of smolt migration",
    breaks       = seq(yrs[1], yrs[2], by = 2),
    minor_breaks = seq(yrs[1], yrs[2], by = 1)
  ) +
  labs(title = "Females") +
  theme_minimal() +
  theme(
    plot.title         = element_text(size = 30, face = "bold", colour = "red"),
    axis.title         = element_text(size = 30, face = "bold"),
    axis.text          = element_text(size = 20, face = "bold"),
    axis.title.y.left  = element_text(colour = "purple"),
    axis.text.y.left   = element_text(colour = "purple"),
    axis.title.y.right = element_text(colour = "red"),
    axis.text.y.right  = element_text(colour = "red"),
    panel.grid.major.x = element_line(linewidth = 0.6),
    panel.grid.minor.x = element_line(linewidth = 0.3)
  )
pF

# Plot males
r_size <- range(data2$size, na.rm = TRUE)
r_mat  <- range(data2$matM, na.rm = TRUE)
offset <- 0.15 * diff(r_size)

to_size <- function(x) (x - r_mat[1]) / diff(r_mat) * diff(r_size) + r_size[1] + offset
to_mat  <- function(y) (y - offset - r_size[1]) / diff(r_size) * diff(r_mat) + r_mat[1]

pM <- ggplot(data2, aes(x = cohort)) +
  geom_line(aes(y = size), colour = "purple", linewidth = 1.5) +
  geom_point(aes(y = size), colour = "purple", size = 5) +
  geom_line(aes(y = to_size(matM)), colour = "blue", linewidth = 1.5) +
  geom_point(aes(y = to_size(matM)), colour = "blue", size = 5) +
  scale_y_continuous(
    name     = "Mean smolt\nscale length (mm)",
    sec.axis = sec_axis(~ to_mat(.), name = "Length-dependent\nmaturation")
  ) +
  scale_x_continuous(
    name         = "Year of smolt migration",
    breaks       = seq(yrs[1], yrs[2], by = 2),
    minor_breaks = seq(yrs[1], yrs[2], by = 1)
  ) +
  labs(title = "Males") +
  theme_minimal() +
  theme(
    plot.title         = element_text(size = 30, face = "bold", colour = "blue"),
    axis.title         = element_text(size = 30, face = "bold"),
    axis.text          = element_text(size = 20, face = "bold"),
    axis.title.y.left  = element_text(colour = "purple"),
    axis.text.y.left   = element_text(colour = "purple"),
    axis.title.y.right = element_text(colour = "blue"),
    axis.text.y.right  = element_text(colour = "blue"),
    panel.grid.major.x = element_line(linewidth = 0.6),
    panel.grid.minor.x = element_line(linewidth = 0.3)
  )
pM

FigA7.2.2 <- pF / pM
FigA7.2.2

ggsave(
  filename = paste0("results/figures/", projects, "/FigA7.2.2.pdf"),
  plot = FigA7.2.2,
  width = 1961/300,
  height = 1440/300,
  units = "in"
)
