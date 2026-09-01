### ============================================================================
### Disentangling length-dependent and length-independent variations
### in survival and maturation in Atlantic salmon
### ----------------------------------------------------------------------------
### File:    tabA4_sum.R
### Purpose: Builds Appendix 4 Table A4.2.
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
pkginstall(c("ggplot2", "dplyr", "tidyr", "coda", "qs", "purrr"))

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

compute_sum <- function(project) {
  MCMC_samples <- load_mcmc(project)
  s <- summary(MCMC_samples)
  
  res <- as.data.frame(s[["statistics"]]) %>%
    tibble::rownames_to_column("Parameter") %>%
    left_join(
      as.data.frame(s[["quantiles"]]) %>% tibble::rownames_to_column("Parameter"),
      by = "Parameter"
    )
  
  write.csv2(x = res, file = paste0("results/tables/", project, "/tabA4.2.csv"), row.names = FALSE)
}

## -----------------------------------------------------------------------------
## 3. Run MCMC diagnostic
## -----------------------------------------------------------------------------
walk(projects, compute_sum)
