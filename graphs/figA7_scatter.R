### ============================================================================
### Creating figure of annual variation in length-dependent effect on vital
### rates and annual mean corresponding length (Appendix 7)
### ============================================================================
rm(list = ls())

pkginstall <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      message(paste("Installation du package manquant:", pkg))
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
      message(paste(pkg, "chargé avec succès"))
    } else {
      message(paste("Le package", pkg, "est déjà installé."))
      library(pkg, character.only = TRUE)
      message(paste(pkg, "chargé avec succès"))
    }
  }
}

packages <- c("dplyr", "tidyr", "ggplot2", "mgcv", "patchwork", "nimble", "ggridges", "plotly", "qs", "ggpubr", "coda", "ggnewscale", "relaimpo", "lm.beta")
pkginstall(packages)

HOME <- rstudioapi::getActiveProject()

# Charger les données
project <- c("M1")
MCMC <- qread(file = paste0(HOME,"/saves/", project, "/MCMC.qs"))
n_chains <- seq_along(MCMC)

MCMC_samples <- lapply(n_chains, function(i) {
  mcmc(MCMC[[i]]$samples)
})
names(MCMC_samples) <- paste0("chain", n_chains)
MCMC_samples <- mcmc.list(MCMC_samples)
qsave(MCMC_samples, file = paste0(HOME,"/saves/", project, "/MCMC_samples.qs"))
MCMC_matrix <- as.matrix(MCMC_samples)

source(paste0(HOME, "/functions/nf_l.R"))
Cnf_l <- compileNimble(nf_l)
source(paste0(HOME, "/functions/nf_pi.R"))
Cnf_pi <- compileNimble(nf_pi)

const <- qread(paste0(HOME, "/data/realdata/const.qs"))
num_cohorts <- const$C
num_iterations <- nrow(MCMC_matrix)
num_class <- const$L

### SURVIVAL
# Paramètres
l1mean <- const$Mean_L1
l1min <- 0
l1max <- 2.5

MCMC_beta1 <- MCMC_matrix[, grep("^Beta1", colnames(MCMC_matrix))]
MCMC_alpha1 <- MCMC_matrix[, grep("^Alpha1", colnames(MCMC_matrix))]
MCMC_mualpha1 <- MCMC_matrix[, grep("^Mu_Alpha1", colnames(MCMC_matrix))]
MCMC_Mu_L1 <- MCMC_matrix[, grep("^Mu_L1\\[\\d+\\]$", colnames(MCMC_matrix))]
MCMC_Sd_L1 <- MCMC_matrix[, grep("^Sd_L1\\[\\d+\\]$", colnames(MCMC_matrix))]

# Calcul
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

data1 <- data.frame(
  cohort = 1996:2019,
  size = L1mean_size,
  surv = survQ50_size
)

r_size <- range(data1$size, na.rm = TRUE)
r_surv <- range(data1$surv, na.rm = TRUE)

offset <- 0.15 * diff(r_size)   # raise the survival series; increase to lift further

to_size <- function(x) (x - r_surv[1]) / diff(r_surv) * diff(r_size) + r_size[1] + offset
to_surv <- function(y) (y - offset - r_size[1]) / diff(r_size) * diff(r_surv) + r_surv[1]

yrs <- range(data1$cohort, na.rm = TRUE)

ggplot(data1, aes(x = cohort)) +
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

### MATURATION
# Paramètres
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

# Calcul
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

data2 <- data.frame(
  cohort = 1996:2019,
  size = L2mean_size,
  matF = matQ50_sizeF,
  matM = matQ50_sizeM
)

library(patchwork)

yrs <- range(data2$cohort, na.rm = TRUE)

make_panel <- function(mat, mat_col, sex_label, offset_frac = 0.15) {
  
  r_size <- range(data2$size, na.rm = TRUE)
  r_mat  <- range(mat, na.rm = TRUE)
  offset <- offset_frac * diff(r_size)
  
  to_size <- function(x) (x - r_mat[1]) / diff(r_mat) * diff(r_size) + r_size[1] + offset
  to_mat  <- function(y) (y - offset - r_size[1]) / diff(r_size) * diff(r_mat) + r_mat[1]
  
  df <- data.frame(cohort = data2$cohort, size = data2$size, mat = mat)
  
  ggplot(df, aes(x = cohort)) +
    geom_line(aes(y = size), colour = "purple", linewidth = 1.5) +
    geom_point(aes(y = size), colour = "purple", size = 5) +
    geom_line(aes(y = to_size(mat)), colour = mat_col, linewidth = 1.5) +
    geom_point(aes(y = to_size(mat)), colour = mat_col, size = 5) +
    scale_y_continuous(
      name     = "Mean smolt\nscale length (mm)",
      sec.axis = sec_axis(~ to_mat(.), name = "Length-dependent\nmaturation")
    ) +
    scale_x_continuous(
      name         = "Year of smolt migration",
      breaks       = seq(yrs[1], yrs[2], by = 2),
      minor_breaks = seq(yrs[1], yrs[2], by = 1)
    ) +
    labs(title = sex_label) +
    theme_minimal() +
    theme(
      plot.title         = element_text(size = 30, face = "bold", colour = mat_col),
      axis.title         = element_text(size = 30, face = "bold"),
      axis.text          = element_text(size = 20, face = "bold"),
      axis.title.y.left  = element_text(colour = "purple"),
      axis.text.y.left   = element_text(colour = "purple"),
      axis.title.y.right = element_text(colour = mat_col),
      axis.text.y.right  = element_text(colour = mat_col),
      panel.grid.major.x = element_line(linewidth = 0.6),
      panel.grid.minor.x = element_line(linewidth = 0.3)
    )
}

pF <- make_panel(data2$matF, "red", "Females") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )

pM <- make_panel(data2$matM, "blue", "Males")

pF / pM
