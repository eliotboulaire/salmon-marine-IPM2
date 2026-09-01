### ============================================================================
### Disentangling length-dependent and length-independent variations
### in survival and maturation in Atlantic salmon
### ----------------------------------------------------------------------------
### File:    tabA4_conv.R
### Purpose: Builds Appendix 4 Table A4.1.
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
pkginstall(c("ggplot2", "dplyr", "tidyr", "ggmcmc", "coda", "qs", "purrr"))

## -----------------------------------------------------------------------------
## 1. Define project settings
## -----------------------------------------------------------------------------
projects <- paste0("M", 0:9)

## -----------------------------------------------------------------------------
## 2. Import MCMC object and functions
## -----------------------------------------------------------------------------
load_mcmc <- function(project) {
  MCMC <- qread(paste0("saves/", project, "/MCMC.qs"))
  MCMC %>%
    seq_along() %>%
    lapply(function(i) mcmc(MCMC[[i]]$samples)) %>%
    mcmc.list()
}

compute_diag <- function(project) {
  MCMC_samples <- load_mcmc(project)
  S <- ggs(MCMC_samples)
  
  diag1 <- ggs_Rhat(S, version_rhat = "BG98", plot = FALSE)
  diag2 <- ggs_geweke(S, plot = FALSE)
  diag3 <- ggs_effective(S, proportion = FALSE, version_effective = "BDA3", plot = FALSE)
  
  diag <- diag2 %>%
    pivot_wider(names_from = Chain, values_from = z, names_prefix = "zchain") %>%
    left_join(diag1, by = "Parameter") %>%
    left_join(diag3, by = "Parameter") %>%
    dplyr::select(Parameter, Rhat, starts_with("zchain"), Effective)
  
  write.csv2(x = diag, file = paste0("results/tables/", project, "/tabA4.1.csv"), row.names = FALSE)
}

## -----------------------------------------------------------------------------
## 3. Run MCMC diagnostic
## -----------------------------------------------------------------------------
walk(projects, compute_diag)