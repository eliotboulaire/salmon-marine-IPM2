### ============================================================================
### Disentangling length-dependent and length-independent variations
### in survival and maturation in Atlantic salmon
### ----------------------------------------------------------------------------
### File:    figA5_fit.R
### Purpose: Create figure of Posterior densities vs data (Appendix 5)
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
pkginstall(c("ggplot2", "dplyr", "tidyr"))

## -----------------------------------------------------------------------------
## 1. Define project settings
## -----------------------------------------------------------------------------
projects <- c("M0")

load_mcmc <- function(project) {
  MCMC <- qread(paste0("saves/", project, "/MCMC.qs"))
  MCMC %>%
    seq_along() %>%
    lapply(function(i) mcmc(MCMC[[i]]$samples)) %>%
    mcmc.list()
}
MCMC_samples <- load_mcmc(projects)
MCMC_matrix <- as.matrix(MCMC_samples)

## -----------------------------------------------------------------------------
## 1. CMR abundance estimates
## -----------------------------------------------------------------------------
# Posterior density
CMR3 <- t(apply(MCMC_matrix[,grep("^Mu_N1|^Mu_N3|^Mu_N4", colnames(MCMC_matrix))],
                2, function(x) quantile(x, c(0.025, 0.25, 0.50, 0.75, 0.975))))
est_N <- data.frame(
  stade = rep(c("Smolt", "1SW", "2SW"), each = 24),
  cohort = as.numeric(rep(1996:2019, times = 3)),
  Q2.5 = as.vector(CMR3[,1]),
  Q25 = as.vector(CMR3[,2]),
  Q50 = as.vector(CMR3[,3]),
  Q75 = as.vector(CMR3[,4]),
  Q97.5 = as.vector(CMR3[,5])
) %>%
  mutate(stade = factor(stade, levels = c("Smolt", "1SW", "2SW")))

# Data
CMR <- read.csv2(file = "data/rawdata/data_Scorff_CMR.csv")
CMR2 <- CMR %>%
  mutate(
    meanlog = log(mean^2/sqrt(mean^2 + sd^2)),
    sdlog = sqrt(log(1+((sd^2)/(mean^2)))),
    stade = factor(stade, levels = c("Smolt", "1SW", "2SW")),
    cohort = case_when(
      grepl("Smolt", stade, ignore.case = TRUE) ~ year,
      grepl("1SW", stade, ignore.case = TRUE) ~ year-1,
      grepl("2SW", stade, ignore.case = TRUE) ~ year-2,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(cohort >= 1996 & cohort <= 2019)

# Dataframe
data_N <- CMR2 %>%
  rowwise() %>%
  mutate(
    Q2.5  = qlnorm(0.025, meanlog, sdlog),
    Q25   = qlnorm(0.25,  meanlog, sdlog),
    Q50   = qlnorm(0.50,  meanlog, sdlog),
    Q75   = qlnorm(0.75,  meanlog, sdlog),
    Q97.5 = qlnorm(0.975, meanlog, sdlog)
  ) %>%
  ungroup()

# Plot
figA5.1.1 <- ggplot() + 
  geom_ribbon(data = data_N, aes(x = cohort, ymin = Q2.5, ymax = Q97.5), fill = "red", alpha = 0.2) +
  geom_line(data = data_N, aes(x = cohort, y = Q50), color = "red", linewidth = 1.2) +
  geom_point(data = data_N, aes(x = cohort, y = Q50, shape = stade), color = "red", size = 4) +
  geom_ribbon(data = est_N, aes(x = cohort, ymin = Q2.5, ymax = Q97.5), fill = "black", alpha = 0.2) +
  geom_line(data = est_N, aes(x = cohort, y = Q50), color = "black", linewidth = 1.2) +
  geom_point(data = est_N, aes(x = cohort, y = Q50, shape = stade), color = "black", size = 4) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
  facet_grid(stade~., scale = "free", switch = "y") +
  labs(x = "Years of migrating smolt", y = "Number of individuals") +
  coord_cartesian(xlim = c(1995, 2020), ylim = c(0,NA), expand = TRUE) +
  scale_x_continuous(
    breaks = seq(min(data_N$cohort), max(data_N$cohort), by = 2),
    minor_breaks = seq(min(data_N$cohort), max(data_N$cohort), by = 1)
  ) +
  scale_y_continuous(n.breaks = 6) +
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
  )
figA5.1.1

ggsave(
  filename = "results/figA5.1.1.pdf",
  plot = figA5.1.1,
  width = 1961/300,
  height = 1440/300,
  units = "in"
)

## -----------------------------------------------------------------------------
## 2. Scale sex data
## -----------------------------------------------------------------------------
# Posterior density
ORE4 <- t(apply(MCMC_matrix[,grep("^Gamma\\d+\\[\\d+, 1\\]", colnames(MCMC_matrix))],
                2, function(x) quantile(x, c(0.025, 0.25, 0.50, 0.75, 0.975))))
est_S <- data.frame(
  stage = rep(c("Smolt", "1SW", "2SW"), each = 24),
  cohort = as.numeric(rep(1996:2019, times = 3)),
  Q2.5 = as.vector(ORE4[,1]),
  Q25 = as.vector(ORE4[,2]),
  Q50 = as.vector(ORE4[,3]),
  Q75 = as.vector(ORE4[,4]),
  Q97.5 = as.vector(ORE4[,5])
) %>%
  mutate(stage = factor(stage, levels = c("Smolt", "1SW", "2SW")))

# Data
ORE <- read.csv2(file = "data/rawdata/data_Scorff_ORE.csv")
ORE2 <- ORE %>%
  mutate(
    stage = case_when(
      grepl("smolt", stade, ignore.case = TRUE) ~ "Smolt",
      grepl("1HM", stade, ignore.case = TRUE) ~ "1SW",
      grepl("2HM", stade, ignore.case = TRUE) ~ "2SW",
      TRUE ~ NA_character_
    ),
    cohort = case_when(
      grepl("Smolt", stage, ignore.case = TRUE) ~ annee,
      grepl("1SW", stage, ignore.case = TRUE) ~ annee-1,
      grepl("2SW", stage, ignore.case = TRUE) ~ annee-2,
      TRUE ~ NA_real_
    ),
    sex = case_when(
      grepl("F", sexegentique, ignore.case = TRUE) ~ "F",
      grepl("M", sexegentique, ignore.case = TRUE) ~ "M",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(cohort >= 1996 & cohort <= 2019) %>%
  filter(stage %in% c("Smolt", "1SW", "2SW")) %>%
  mutate(stage = factor(stage, levels = c("Smolt", "1SW", "2SW")))

# Dataframe
data_S <- ORE2 %>%
  group_by(stage, cohort) %>%
  summarise(
    Fem = sum(sexegentique %in% c("F"), na.rm = TRUE),
    Ma = sum(sexegentique %in% c("M"), na.rm = TRUE),
    Tot = sum(sexegentique %in% c("F", "M"), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(
    Q2.5  = qbeta(0.025, (Fem + 1), ((Tot - Fem) + 1)),
    Q25   = qbeta(0.25, (Fem + 1), ((Tot - Fem) + 1)),
    Q50   = qbeta(0.50, (Fem + 1), ((Tot - Fem) + 1)),
    Q75   = qbeta(0.75, (Fem + 1), ((Tot - Fem) + 1)),
    Q97.5 = qbeta(0.975, (Fem + 1), ((Tot - Fem) + 1))
  ) %>%
  ungroup()

# Plot
figA5.1.2 <- ggplot() + 
  geom_ribbon(data = data_S, aes(x = cohort, ymin = Q2.5, ymax = Q97.5), fill = "red", alpha = 0.2) +
  geom_line(data = data_S, aes(x = cohort, y = Q50), color = "red", linewidth = 1.2) +
  geom_point(data = data_S, aes(x = cohort, y = Q50, shape = stage), color = "red", size = 4) +
  geom_ribbon(data = est_S, aes(x = cohort, ymin = Q2.5, ymax = Q97.5), fill = "black", alpha = 0.2) +
  geom_line(data = est_S, aes(x = cohort, y = Q50), color = "black", linewidth = 1.2) +
  geom_point(data = est_S, aes(x = cohort, y = Q50, shape = stage), color = "black", size = 4) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
  facet_grid(stage~., scale = "free", switch = "y") +
  labs(x = "Years of migrating smolt", y = "Proportion of female") +
  coord_cartesian(xlim = c(1995, 2020), ylim = c(0,1), expand = TRUE) +
  scale_x_continuous(
    breaks = seq(min(data_S$cohort), max(data_S$cohort), by = 2),
    minor_breaks = seq(min(data_S$cohort), max(data_S$cohort), by = 1)
  ) +
  scale_y_continuous(n.breaks = 6) +
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
  )
figA5.1.2

ggsave(
  filename = "results/figA5.1.2.pdf",
  plot = figA5.1.2,
  width = 1961/300,
  height = 1440/300,
  units = "in"
)

## -----------------------------------------------------------------------------
## 3. Scale length data
## -----------------------------------------------------------------------------
# Data
ORE3 <- ORE2 %>%
  filter(!is.na(edge)) %>%
  mutate(
    riv = case_when(
      grepl("Smolt", stage, ignore.case = TRUE) ~ edge,
      grepl("SW", stage, ignore.case = TRUE) ~ Transition,
      TRUE ~ NA_real_),
    sum1 = case_when(
      grepl("SW", stage, ignore.case = TRUE) ~ X1SW_bg,
      TRUE ~ NA_real_),
    win = case_when(
      grepl("SW", stage, ignore.case = TRUE) ~ edge,
      TRUE ~ NA_real_)
  )
data_L <- bind_rows(
  ORE3 %>%
    filter(stage == "Smolt" & !is.na(riv)) %>%
    mutate(stage2 = "Smolt",
           length = riv),
  ORE3 %>%
    filter(grepl("SW", stage) & !is.na(riv)) %>%
    mutate(stage2 = "Surviving smolt",
           length = riv),
  ORE3 %>%
    filter(grepl("SW", stage) & !is.na(sum1)) %>%
    mutate(stage2 = "Post-smolt",
           length = sum1),
  ORE3 %>%
    filter(stage == "1SW" & !is.na(sum1)) %>%
    mutate(stage2 = "Maturing
post-smolt",
           length = sum1),
  ORE3 %>%
    filter(stage == "2SW" & !is.na(sum1)) %>%
    mutate(stage2 = "Non-maturing
post-smolt",
           length = sum1),
  ORE3 %>%
    filter(stage == "1SW" & !is.na(edge)) %>%
    mutate(stage2 = "1SW",
           length = edge),
  ORE3 %>%
    filter(stage == "2SW" & !is.na(edge)) %>%
    mutate(stage2 = "2SW",
           length = edge),
) %>%
  group_by(stage2, cohort) %>%
  summarise(
    Q2.5  = quantile(length, 0.025),
    Q25   = quantile(length, 0.25),
    Q50   = quantile(length, 0.5),
    Q75   = quantile(length, 0.75),
    Q97.5 = quantile(length, 0.975),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  filter(stage2 %in% c("Smolt", "Surviving smolt", "Post-smolt", "Maturing
post-smolt", "Non-maturing
post-smolt")) %>%
  mutate(stage2 = factor(stage2, levels = c("Smolt", "Surviving smolt", "Post-smolt", "Maturing
post-smolt", "Non-maturing
post-smolt")))

# Posterior density
ORE5 <- t(apply(cbind(
  colMeans(MCMC_matrix[,grep("^Mu_L", colnames(MCMC_matrix))]),
  colMeans(MCMC_matrix[,grep("^Sd_L", colnames(MCMC_matrix))])
), 1, function(x) qnorm(c(0.025, 0.25, 0.50, 0.75, 0.975), mean = x[1], sd = x[2])))

est_L <- data.frame(
  stage2 = rep(c("Smolt", "Surviving smolt", "Post-smolt", "Maturing
post-smolt", "Non-maturing
post-smolt"), each = 24),
  cohort = as.numeric(rep(1996:2019, times = 5)),
  Q2.5 = as.vector(ORE5[,1]),
  Q25 = as.vector(ORE5[,2]),
  Q50 = as.vector(ORE5[,3]),
  Q75 = as.vector(ORE5[,4]),
  Q97.5 = as.vector(ORE5[,5])
) %>%
  mutate(stage2 = factor(stage2, levels = c("Smolt", "Surviving smolt", "Post-smolt", "Maturing
post-smolt", "Non-maturing
post-smolt")))

# Plot
figA5.1.3 <- ggplot() + 
  geom_ribbon(data = data_L, aes(x = cohort, ymin = Q2.5, ymax = Q97.5), fill = "red", alpha = 0.2) +
  geom_line(data = data_L, aes(x = cohort, y = Q50), color = "red", linewidth = 1.2) +
  geom_point(data = data_L, aes(x = cohort, y = Q50, shape = stage2), color = "red", size = 4) +
  geom_ribbon(data = est_L, aes(x = cohort, ymin = Q2.5, ymax = Q97.5), fill = "black", alpha = 0.2) +
  geom_line(data = est_L, aes(x = cohort, y = Q50), color = "black", linewidth = 1.2) +
  geom_point(data = est_L, aes(x = cohort, y = Q50, shape = stage2), color = "black", size = 4) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
  facet_grid(stage2~., scale = "free", switch = "y") +
  labs(x = "Years of migrating smolt", y = "Scale length (mm)") +
  coord_cartesian(xlim = c(1995, 2020), ylim = c(0,NA), expand = TRUE) +
  scale_x_continuous(
    breaks = seq(min(data_L$cohort), max(data_L$cohort), by = 2),
    minor_breaks = seq(min(data_L$cohort), max(data_L$cohort), by = 1)
  ) +
  scale_y_continuous(n.breaks = 6) +
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
  )
figA5.1.3

ggsave(
  filename = "results/figA5.1.3.pdf",
  plot = figA5.1.3,
  width = 1961/300,
  height = 1440/300,
  units = "in"
)
