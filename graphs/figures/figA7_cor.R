### ============================================================================
### Creating figure of Pearson correlations distribution acorss posterior
### samples between the sex in maturation (Appendix 7)
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

packages <- c("dplyr", "tidyr", "ggplot2", "mgcv", "patchwork", "nimble", "ggridges", "plotly", "qs", "coda")
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

MCMC_alpha1 <- MCMC_matrix[, grep("^Alpha1", colnames(MCMC_matrix))]
MCMC_alpha2 <- MCMC_matrix[, grep("^Alpha2", colnames(MCMC_matrix))]
MCMC_alpha2F <- MCMC_matrix[, grep("^Alpha2", colnames(MCMC_matrix))][,1:24]
MCMC_alpha2M <- MCMC_matrix[, grep("^Alpha2", colnames(MCMC_matrix))][,25:48]

truc <- mapply(
  function(x, y) cor.test(x, y, method = "pearson"),
  split(MCMC_alpha2F, row(MCMC_alpha2F)),
  split(MCMC_alpha2M, row(MCMC_alpha2M)),
  SIMPLIFY = FALSE
)
cor_values <- sapply(truc, function(t) t$estimate)
truc2 <- data.frame(
  samples = 1:30000,
  pearson = cor_values
)
p_neg <- mean(truc2$pearson <= 0)

plot_truc1 <- ggplot(truc2, aes(x = pearson)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 100,
                 fill = "white",
                 color = "black") +
  geom_histogram(data = subset(truc2, pearson <= 0),
                 aes(y = after_stat(density * p_neg)),
                 bins = 100,
                 fill = "red",
                 alpha = 0.4) +
  # annotate(
  #   "text",
  #   x = -0.2,
  #   y = 1,
  #   label = paste0("P(r <= 0) = ", round(100 * p_neg, 1), "%"),
  #   color = "red",
  #   size = 10
  # ) +
  labs(x = "Pearson correlation (r)", y = "Density") + 
  geom_vline(xintercept = 0, color = "red", linewidth = 1.2, linetype = "dashed") +
  # geom_rect(aes(xmin = quantile(truc2$pearson, probs = c(0.25)), xmax = quantile(truc2$pearson, probs = c(0.75)), ymin = -1, ymax = 4), fill = NA, color = "green", linewidth = 1) +
  geom_vline(xintercept = quantile(truc2$pearson, probs = c(0.5)), color = "black", linewidth = 1.2) +
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

### MATURATION
source(paste0(HOME, "/functions/nf_l.R"))
Cnf_l <- compileNimble(nf_l)
source(paste0(HOME, "/functions/nf_pi.R"))
Cnf_pi <- compileNimble(nf_pi)

const <- qread(paste0(HOME, "/data/realdata/const.qs"))
num_cohorts <- const$C
num_iterations <- nrow(MCMC_matrix)
num_class <- const$L

# Paramètres
l2mean <- const$Mean_L2
l2min <- 1.5
l2max <- 4
MCMC_Mu_L2 <- MCMC_matrix[, grep("^Mu_L2\\[\\d+\\]$", colnames(MCMC_matrix))]
MCMC_Sd_L2 <- MCMC_matrix[, grep("^Sd_L2\\[\\d+\\]$", colnames(MCMC_matrix))]
MCMC_beta2 <- MCMC_matrix[, grep("^Beta2", colnames(MCMC_matrix))]

# Calcul
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

truc <- mapply(
  function(x, y) cor.test(x, y, method = "pearson"),
  split(t(matM_baseF), row(t(matM_baseF))),
  split(t(matM_baseM), row(t(matM_baseM))),
  SIMPLIFY = FALSE
)
cor_values <- sapply(truc, function(t) t$estimate)
truc2 <- data.frame(
  samples = 1:30000,
  pearson = cor_values
)
p_neg <- mean(truc2$pearson <= 0)

plot_truc2 <- ggplot(truc2, aes(x = pearson)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 100,
                 fill = "white",
                 color = "black") +
  geom_histogram(data = subset(truc2, pearson <= 0),
                 aes(y = after_stat(density * p_neg)),
                 bins = 100,
                 fill = "red",
                 alpha = 0.4) +
  # annotate(
  #   "text",
  #   x = -0.2,
  #   y = 1,
  #   label = paste0("P(r <= 0) = ", round(100 * p_neg, 1), "%"),
  #   color = "red",
  #   size = 10
  # ) +
  labs(x = "Pearson correlation (r)", y = "Density") + 
  geom_vline(xintercept = 0, color = "red", linewidth = 1.2, linetype = "dashed") +
  # geom_rect(aes(xmin = quantile(truc2$pearson, probs = c(0.25)), xmax = quantile(truc2$pearson, probs = c(0.75)), ymin = -1, ymax = 4), fill = NA, color = "green", linewidth = 1) +
  geom_vline(xintercept = quantile(truc2$pearson, probs = c(0.5)), color = "black", linewidth = 1.2) +
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

save_plot <- function(filename, plot, width_px = 1920, height_px = 1177, ori_dpi = 72, new_dpi = 300, bg = "white") {
  rate_dpi <- (new_dpi/ori_dpi)
  
  ggsave(
    filename,
    plot = plot,
    width  = (width_px*rate_dpi),
    height = (height_px*rate_dpi),
    units = "px",
    dpi = new_dpi,
    bg = bg
  )
}
save_plot(filename = "plot_truc1.png", plot = plot_truc1, width_px = 1920, height_px = 1177, ori_dpi = 72, new_dpi = 300, bg = "transparent")
save_plot(filename = "plot_truc2.png", plot = plot_truc2, width_px = 1920, height_px = 1177, ori_dpi = 72, new_dpi = 300, bg = "transparent")



MCMC_Mu_L1 <- MCMC_matrix[, grep("^Mu_L1\\[\\d+\\]$", colnames(MCMC_matrix))]
MCMC_Sd_L1 <- MCMC_matrix[, grep("^Sd_L1\\[\\d+\\]$", colnames(MCMC_matrix))]
Mu_L1 <- t(apply(MCMC_Mu_L1, 2, quantile, probs = c(0.25, 0.5, 0.75)))
Sd_L1 <- t(apply(MCMC_Sd_L1, 2, quantile, probs = c(0.25, 0.5, 0.75)))
Mu_L2 <- t(apply(MCMC_Mu_L2, 2, quantile, probs = c(0.25, 0.5, 0.75)))
Sd_L2 <- t(apply(MCMC_Sd_L2, 2, quantile, probs = c(0.25, 0.5, 0.75)))
Alpha2F <- t(apply(MCMC_alpha2, 2, quantile, probs = c(0.25, 0.5, 0.75)))[1:24,]
Alpha2M <- t(apply(MCMC_alpha2, 2, quantile, probs = c(0.25, 0.5, 0.75)))[25:48,]

ggplot() + 
  geom_point(aes(x = 1996:2019, y = Mu_L1[,2]), color = "darkgreen", size = 4) + 
  geom_smooth(aes(x = 1996:2019, y = Mu_L1[,2]), method = "lm", color = "green", fill = "green", alpha = 0.2, linewidth = 1.2) +
  theme_minimal()

ggplot() + 
  geom_point(aes(x = 1996:2019, y = Sd_L1[,2]), color = "darkgreen", size = 4) + 
  geom_smooth(aes(x = 1996:2019, y = Sd_L1[,2]), method = "lm", color = "green", fill = "green", alpha = 0.2, linewidth = 1.2) +
  theme_minimal()

ggplot() + 
  geom_point(aes(x = 1996:2019, y = Mu_L2[,2]), color = "purple4", size = 4) + 
  geom_smooth(aes(x = 1996:2019, y = Mu_L2[,2]), method = "lm", color = "purple", fill = "purple", alpha = 0.2, linewidth = 1.2) +
  theme_minimal()

ggplot() + 
  geom_point(aes(x = 1996:2019, y = Sd_L2[,2]), color = "purple4", size = 4) + 
  geom_smooth(aes(x = 1996:2019, y = Sd_L2[,2]), method = "lm", color = "purple", fill = "purple", alpha = 0.2, linewidth = 1.2) +
  theme_minimal()

ggplot() + 
  geom_point(aes(x = Alpha2F[,2], y = Alpha2M[,2]), color = "purple4", size = 4) + 
  geom_smooth(aes(x = Alpha2F[,2], y = Alpha2M[,2]), method = "lm", color = "purple", fill = "purple", alpha = 0.2, linewidth = 1.2) +
  theme_minimal()
