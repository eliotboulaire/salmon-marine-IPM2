### ============================================================================
### Disentangling length-dependent and length-independent variations
### in survival and maturation in Atlantic salmon
### ----------------------------------------------------------------------------
### File:    tab1_loo.R
### Purpose: Builds Table 1.
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
pkginstall(c("coda", "loo", "qs", "dplyr"))

## -----------------------------------------------------------------------------
## 1. Define project settings
## -----------------------------------------------------------------------------
projects <- paste0("M", 0:9)

## -----------------------------------------------------------------------------
## 2. Run loo comparison
## -----------------------------------------------------------------------------
compute_loo <- function(project) {
  MCMC <- qread(paste0("saves/", project, "/MCMC.qs"))
  
  loglik_matrix <- MCMC %>%
    seq_along() %>%
    lapply(function(i) mcmc(MCMC[[i]]$samples2)) %>%
    mcmc.list() %>%
    as.matrix() %>%
    (\(m) m[complete.cases(m), ])()
  
  loo_result <- loo(loglik_matrix)
  
  # Return a named vector or small data frame
  as.data.frame(t(loo_result$estimates[, 1])) %>%
    cbind(model = project)
}

df_loo <- map(projects, compute_loo) %>%
  list_rbind() %>%
  mutate(dlooic = looic - looic[1])
write.csv2(df_loo, "results/tables/tab1.csv", row.names = FALSE)
