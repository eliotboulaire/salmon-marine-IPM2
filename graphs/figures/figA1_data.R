### ============================================================================
### Disentangling length-dependent and length-independent variations
### in survival and maturation in Atlantic salmon
### ----------------------------------------------------------------------------
### File:    figA1_data.R
### Purpose: Create data figure A1.1; A1.2 & A1.3
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
## 1. CMR abundance estimates
## -----------------------------------------------------------------------------
CMR <- read.csv2(file = "data/rawdata/data_Scorff_CMR.csv")
CMR2 <- CMR %>%
  mutate(
    meanlog = log(mean^2/sqrt(mean^2 + sd^2)),
    sdlog = sqrt(log(1+((sd^2)/(mean^2)))),
    stage = factor(stage, levels = c("Smolt", "1SW", "2SW"))
  )

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

figA1.1 <- ggplot() + 
  geom_ribbon(data = data_N, aes(x = cohort, ymin = Q2.5, ymax = Q97.5), fill = "red", alpha = 0.2) +
  geom_line(data = data_N, aes(x = cohort, y = Q50), color = "red", linewidth = 0.8) +
  geom_point(data = data_N, aes(x = cohort, y = Q50), color = "red", size = 1.8) +
  facet_grid(stage~., scale = "free", switch = "y") +
  labs(x = "Year of smolt migration", y = "Number of individuals") +
  coord_cartesian(xlim = c(1995, 2020), ylim = c(NA, NA), expand = TRUE) +
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
figA1.1

ggsave(
  filename = "results/figA1.1.pdf",
  plot = figA1.1,
  width = 1961/300,
  height = 1440/300,
  units = "in"
)

## -----------------------------------------------------------------------------
## 2. Scale sex data
## -----------------------------------------------------------------------------
ORE <- read.csv2(file = "data/rawdata/data_Scorff_ORE.csv")
ORE2 <- ORE %>%
  mutate(
    stage2 = case_when(
      grepl("smolt", stage, ignore.case = TRUE) ~ "Smolt",
      grepl("1SW", age, ignore.case = TRUE) ~ "1SW",
      grepl("2SW", age, ignore.case = TRUE) ~ "2SW",
      TRUE ~ NA_character_
    ),
    stage2 = factor(stage2, levels = c("Smolt", "1SW", "2SW"))) %>%
  filter(stage2 %in% c("Smolt", "1SW", "2SW"))

data_S <- ORE2 %>%
  group_by(stage2, cohort) %>%
  summarise(
    Fem = sum(sexe %in% c("F"), na.rm = TRUE),
    Ma = sum(sexe %in% c("M"), na.rm = TRUE),
    Tot = sum(sexe %in% c("F", "M"), na.rm = TRUE),
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

figA1.2 <- ggplot() + 
  geom_ribbon(data = data_S, aes(x = cohort, ymin = Q2.5, ymax = Q97.5), fill = "red", alpha = 0.2) +
  geom_line(data = data_S, aes(x = cohort, y = Q50), color = "red", linewidth = 0.8) +
  geom_point(data = data_S, aes(x = cohort, y = Q50), color = "red", size = 1.8) +
  facet_grid(stage2~., scale = "free", switch = "y") +
  labs(x = "Year of smolt migration", y = "Proportion of female") +
  coord_cartesian(xlim = c(1995, 2020), ylim = c(NA, NA), expand = TRUE) +
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
figA1.2

ggsave(
  filename = "results/figA1.2.pdf",
  plot = figA1.2,
  width = 1961/300,
  height = 1440/300,
  units = "in"
)

## -----------------------------------------------------------------------------
## 3. Scale length data
## -----------------------------------------------------------------------------
ORE <- read.csv2(file = "data/rawdata/data_Scorff_ORE.csv")

data_L <- bind_rows(
  ORE %>%
    filter(stage == "Smolt" & !is.na(migration)) %>%
    mutate(stage2 = "Smolt",
           length = migration),
  ORE %>%
    filter(stage == "Adult" & !is.na(migration)) %>%
    mutate(stage2 = "Surviving
smolt",
           length = migration),
  ORE %>%
    filter(stage == "Adult" & !is.na(end1sum)) %>%
    mutate(stage2 = "Post-smolt",
           length = end1sum),
  ORE %>%
    filter(age == "1SW" & !is.na(end1sum)) %>%
    mutate(stage2 = "Maturing
post-smolt",
           length = end1sum),
  ORE %>%
    filter(age == "2SW" & !is.na(end1sum)) %>%
    mutate(stage2 = "Non-maturing
post-smolt",
           length = end1sum)
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
  filter(stage2 %in% c("Smolt", "Surviving
smolt", "Post-smolt", "Maturing
post-smolt", "Non-maturing
post-smolt")) %>%
  mutate(stage2 = factor(stage2, levels = c("Smolt", "Surviving
smolt", "Post-smolt", "Maturing
post-smolt", "Non-maturing
post-smolt")))

figA1.3 <- ggplot() + 
  geom_ribbon(data = data_L, aes(x = cohort, ymin = Q2.5, ymax = Q97.5), fill = "red", alpha = 0.2) +
  geom_line(data = data_L, aes(x = cohort, y = Q50), color = "red", linewidth = 0.8) +
  geom_point(data = data_L, aes(x = cohort, y = Q50), color = "red", size = 1.8) +
  facet_grid(stage2~., scale = "free", switch = "y") +
  labs(x = "Year of smolt migration", y = "Scale length (mm)") +
  coord_cartesian(xlim = c(1995, 2020), ylim = c(NA, NA), expand = TRUE) +
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
figA1.3

ggsave(
  filename = "results/figA1.3.pdf",
  plot = figA1.3,
  width = 1961/300,
  height = 1440/300,
  units = "in"
)
