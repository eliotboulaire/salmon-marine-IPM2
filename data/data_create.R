### ============================================================================
### Disentangling length-dependent and length-independent variations
### in survival and maturation in Atlantic salmon
### ----------------------------------------------------------------------------
### File:    data_create.R
### Purpose: Reads both data_Scorff files to build the datasets used
###          by the model (LogN, Ns, S).
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
pkginstall(c("dplyr", "tidyr"))

## -----------------------------------------------------------------------------
## 1. Import data
## -----------------------------------------------------------------------------
source("functions/rfunctions/f_fillNA.R")

ORE <- read.csv2(file = "data/rawdata/data_Scorff_ORE.csv")
CMR <- read.csv2(file = "data/rawdata/data_Scorff_CMR.csv")

## -----------------------------------------------------------------------------
## 2. Scales length data
## -----------------------------------------------------------------------------
# Data create
S1 <- ORE %>%
  mutate(index = row_number()) %>%
  select(index, stage, age, cohort, migration) %>%
  pivot_wider(names_from = cohort, values_from = migration) %>%
  filter(!if_all(-c(index, stage, age), is.na)) %>%
  filter(stage == "Smolt") %>%
  select(-c(index, stage, age)) %>%
  f_fillNA()
write.csv2(S1, file = "data/rawdata/scales/S1.csv", row.names = FALSE)

S1surv <- ORE %>%
  mutate(index = row_number()) %>%
  select(index, stage, age, cohort, migration) %>%
  pivot_wider(names_from = cohort, values_from = migration) %>%
  filter(!if_all(-c(index, stage, age), is.na)) %>%
  filter(stage == "Adult") %>%
  select(-c(index, stage, age)) %>%
  f_fillNA()
write.csv2(S1surv, file = "data/rawdata/scales/S1surv.csv", row.names = FALSE)

S2 <- ORE %>%
  mutate(index = row_number()) %>%
  select(index, stage, age, cohort, end1sum) %>%
  pivot_wider(names_from = cohort, values_from = end1sum) %>%
  filter(!if_all(-c(index, stage, age), is.na)) %>%
  filter(stage == "Adult") %>%
  select(-c(index, stage, age)) %>%
  f_fillNA()
write.csv2(S2, file = "data/rawdata/scales/S2.csv", row.names = FALSE)

S2m <- ORE %>%
  mutate(index = row_number()) %>%
  select(index, stage, age, cohort, end1sum) %>%
  pivot_wider(names_from = cohort, values_from = end1sum) %>%
  filter(!if_all(-c(index, stage, age), is.na)) %>%
  filter(age == "1SW") %>%
  select(-c(index, stage, age)) %>%
  f_fillNA()
write.csv2(S2m, file = "data/rawdata/scales/S2m.csv", row.names = FALSE)

S2nm <- ORE %>%
  mutate(index = row_number()) %>%
  select(index, stage, age, cohort, end1sum) %>%
  pivot_wider(names_from = cohort, values_from = end1sum) %>%
  filter(!if_all(-c(index, stage, age), is.na)) %>%
  filter(age == "2SW") %>%
  select(-c(index, stage, age)) %>%
  f_fillNA()
write.csv2(S2nm, file = "data/rawdata/scales/S2nm.csv", row.names = FALSE)

## -----------------------------------------------------------------------------
## 3. Scales sex data
## -----------------------------------------------------------------------------
# Data create
Ns1 <- ORE %>%
  filter(stage == "Smolt") %>%
  filter(sexe %in% c("F", "M")) %>%
  select(cohort, sexe) %>%
  group_by(cohort) %>%
  summarise(
    females = as.numeric(sum(sexe == "F")),
    males = as.numeric(sum(sexe == "M")),
    total = as.numeric(n()))
write.csv2(Ns1, file = "data/rawdata/sexes/Ns1.csv", row.names = FALSE)

Ns3 <- ORE %>%
  filter(age == "1SW") %>%
  filter(sexe %in% c("F", "M")) %>%
  select(cohort, sexe) %>%
  group_by(cohort) %>%
  summarise(
    females = as.numeric(sum(sexe == "F")),
    males = as.numeric(sum(sexe == "M")),
    total = as.numeric(n()))
write.csv2(Ns3, file = "data/rawdata/sexes/Ns3.csv", row.names = FALSE)

Ns4 <- ORE %>%
  filter(age == "2SW") %>%
  filter(sexe %in% c("F", "M")) %>%
  select(cohort, sexe) %>%
  group_by(cohort) %>%
  summarise(
    females = as.numeric(sum(sexe == "F")),
    males = as.numeric(sum(sexe == "M")),
    total = as.numeric(n()))
write.csv2(Ns4, file = "data/rawdata/sexes/Ns4.csv", row.names = FALSE)

## -----------------------------------------------------------------------------
## 4. CMR abundances
## -----------------------------------------------------------------------------
# Data create
LogN1 <- CMR %>%
  filter(stage == "Smolt") %>%
  select(cohort, mean, sd) %>%
  mutate(log_sd = sqrt(log(1+((sd^2)/(mean^2)))))
write.csv2(LogN1, file = "data/rawdata/abundances/LogN1.csv", row.names = FALSE)

LogN3 <- CMR %>%
  filter(stage == "1SW") %>%
  select(cohort, mean, sd)%>%
  mutate(log_sd = sqrt(log(1+((sd^2)/(mean^2)))))
write.csv2(LogN3, file = "data/rawdata/abundances/LogN3.csv", row.names = FALSE)

LogN4 <- CMR %>%
  filter(stage == "2SW") %>%
  select(cohort, mean, sd) %>%
  mutate(log_sd = sqrt(log(1+((sd^2)/(mean^2)))))
write.csv2(LogN4, file = "data/rawdata/abundances/LogN4.csv", row.names = FALSE)