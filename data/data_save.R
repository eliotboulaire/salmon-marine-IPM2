### ============================================================================
### Disentangling length-dependent and length-independent variations
### in survival and maturation in Atlantic salmon
### ----------------------------------------------------------------------------
### File:    data_save.R
### Purpose: Read each csv files created with data_create.R 
###          Assemble objects used by the NIMBLE model: 
###          (1) the "data" list (observations)
###          (2) the `const` list (dimensions, sample sizes,
###              length bounds, fixed survivals).
###          data and const are saved in .qs to data/realdata/
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
pkginstall(c("dplyr", "tidyr", "qs"))

## -----------------------------------------------------------------------------
## 1. Import data
## -----------------------------------------------------------------------------
# Import Scales length data
S1 <- read.csv2(file = "data/rawdata/scales/S1.csv") %>%
  as.matrix()

S1surv <- read.csv2(file = "data/rawdata/scales/S1surv.csv") %>%
  as.matrix()

S2 <- read.csv2(file = "data/rawdata/scales/S2.csv") %>%
  as.matrix()

S2m <- read.csv2(file = "data/rawdata/scales/S2m.csv") %>%
  as.matrix()

S2nm <- read.csv2(file = "data/rawdata/scales/S2nm.csv") %>%
  as.matrix()

# Import Scales sex data
Ns1 <- read.csv2(file = "data/rawdata/sexes/Ns1.csv") %>%
  mutate(females = as.numeric(females),
         males = as.numeric(males),
         total = as.numeric(total)) %>%
  select(-cohort) %>%
  as.matrix()

Ns3 <- read.csv2(file = "data/rawdata/sexes/Ns3.csv") %>%
  mutate(females = as.numeric(females),
         males = as.numeric(males),
         total = as.numeric(total)) %>%
  select(-cohort) %>%
  as.matrix()

Ns4 <- read.csv2(file = "data/rawdata/sexes/Ns4.csv") %>%
  mutate(females = as.numeric(females),
         males = as.numeric(males),
         total = as.numeric(total)) %>%
  select(-cohort) %>%
  as.matrix()

# Import CMR abundance estimates
LogN1 <- read.csv2(file = "data/rawdata/abundances/LogN1.csv") %>%
  select(-cohort) %>%
  as.matrix()

LogN3 <- read.csv2(file = "data/rawdata/abundances/LogN3.csv") %>%
  select(-cohort) %>%
  as.matrix()

LogN4 <- read.csv2(file = "data/rawdata/abundances/LogN4.csv") %>%
  select(-cohort) %>%
  as.matrix()

## -----------------------------------------------------------------------------
## 2. Define "data" file
## -----------------------------------------------------------------------------
data <- list(S1 = S1, LogN1 = LogN1[, 1], Ns1 = Ns1[, 1:2],
             S1surv = S1surv,
             S2 = S2,
             S2m = S2m,
             S2nm = S2nm,
             LogN3 = LogN3[, 1], Ns3 = Ns3[, 1:2],
             LogN4 = LogN4[, 1], Ns4 = Ns4[, 1:2])
qs::qsave(data, file = paste0(HOME, "/data/realdata/data.qs"))

## -----------------------------------------------------------------------------
## 3. Define "const" file
## -----------------------------------------------------------------------------
source(paste0(HOME, "/functions/f_round.R"))

I1 <- as.numeric(apply(S1, 2, function(x) sum(!is.na(x))))
Min_L1 <- f_round(min(S1, na.rm = TRUE), digit = 1, method = 1)
Max_L1 <- f_round(max(S1, na.rm = TRUE), digit = 1, method = 2)
Mean_L1 <- round(mean(S1, na.rm = TRUE), digits = 3)

I1surv <- as.numeric(apply(S1surv, 2, function(x) sum(!is.na(x))))

I2 <- as.numeric(apply(S2, 2, function(x) sum(!is.na(x))))
Min_L2 <- f_round(min(S2, na.rm = TRUE), digit = 1, method = 1)
Max_L2 <- f_round(max(S2, na.rm = TRUE), digit = 1, method = 2)
Mean_L2 <- round(mean(S2, na.rm = TRUE), digits = 3)

I2m <- as.numeric(apply(S2m, 2, function(x) sum(!is.na(x))))

I2nm <- as.numeric(apply(S2nm, 2, function(x) sum(!is.na(x))))

C <- nrow(LogN1)
L <- 60
S <- 2

Theta3 <- exp(-0.03*9)
Theta4 <- exp(-0.03*17)

const <- list(Sd_LogN1 = LogN1[, 3], Nt1 = Ns1[, 3],
              Sd_LogN3 = LogN3[, 3], Nt3 = Ns3[, 3],
              Sd_LogN4 = LogN4[, 3], Nt4 = Ns4[, 3],
              I1 = I1, I1surv = I1surv, I2 = I2, I2m = I2m, I2nm = I2nm, L = L, C = C, S = S,
              Theta3 = Theta3, Theta4 = Theta4,
              Min_L1 = Min_L1, Max_L1 = Max_L1, Mean_L1 = Mean_L1,
              Min_L2 = Min_L2, Max_L2 = Max_L2, Mean_L2 = Mean_L2)
qs::qsave(const, file = paste0(HOME, "/data/realdata/const.qs"))