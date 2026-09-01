### ============================================================================
### Disentangling length-dependent and length-independent variations
### in survival and maturation in Atlantic salmon
### ----------------------------------------------------------------------------
### File:    tabA1_data.R
### Purpose: Builds Appendix 1 Tables (A1.2 & A1.3).
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
pkginstall(c("ggplot2", "dplyr", "tidyr"))

## -----------------------------------------------------------------------------
## 1. Import data
## -----------------------------------------------------------------------------
CMR <- read.csv2(file = "data/rawdata/data_Scorff_CMR.csv")
ORE <- read.csv2(file = "data/rawdata/data_Scorff_ORE.csv")

## -----------------------------------------------------------------------------
## 2. Scale sex data
## -----------------------------------------------------------------------------
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

tabA1.2 <- ORE2 %>%
  group_by(stage2, cohort) %>%
  summarise(
    Tot = sum(sexe %in% c("F", "M"), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(values_from = Tot, names_from = stage2) %>%
  rename(years = cohort)
write.csv2(tabA1.2, file = "results/tables/tabA1.2.csv", row.names = FALSE)

## -----------------------------------------------------------------------------
## 3. Scale length data
## -----------------------------------------------------------------------------
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

tabA1.3 <- ORE2 %>%
  group_by(stage2, cohort) %>%
  summarise(
    Tot = sum(!is.na(migration) | !is.na(end1sum), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(values_from = Tot, names_from = stage2) %>%
  rename(years = cohort)
write.csv2(tabA1.3, file = "results/tables/tabA1.3.csv", row.names = FALSE)
