### ============================================================================
### Disentangling length-dependent and length-independent variations
### in survival and maturation in Atlantic salmon
### ----------------------------------------------------------------------------
### File:    data_save.R
### Purpose: Read "data" and "const" object from data_save.R 
###          Assemble objects used by the NIMBLE model: 
###          (1) the "inits" list (initial values for MCMC chains to start)
###              "inits_chain" are initial values only for certain parameters
###              "inits_nchain" are initial values for all parameters and chains
###          (2) the `monitor` list (monitored values to keep for results)
###              "monitor1" are monitored parameters to compute results
###              "monitor2" are monitored parameters to compute loo comparisons
###          "inits" and "monitor" are saved in .qs to data/realdata/[M0-M9]/
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
pkginstall(c("nimble", "parallel", "coda", "qs"))

## -----------------------------------------------------------------------------
## 1. Define project settings
## -----------------------------------------------------------------------------
projects <- paste0("M", 0:9)
n_chains <- 6
seeds <- sample(1:1e9, n_chains)

## -----------------------------------------------------------------------------
## 2. Import model function, data and const
## -----------------------------------------------------------------------------
# Data and const
data <- qread("data/realdata/data.qs")
const <- qread("data/realdata/const.qs")

# Model functions
source("functions/nf_res.R")
source("functions/nf_l.R")
source("functions/nf_pi.R")

## -----------------------------------------------------------------------------
## 3. Define inits
## -----------------------------------------------------------------------------

## Inits chains
source(paste0(HOME, "/functions/f_geninit.R"))

for (i in seq_len(n_chain)) {
  set.seed(seeds[i])
  inits <- f_geninit(n = c(24,24,1,1,1,1,1,1,1,1,1,1,1,1,2,2))
  qsave(inits, file = file.path("data", "realdata", paste0("inits_chain", i, ".qs")))
}

## Inits allchains
for (project in projects) {
  if (file.exists(file.path("data","realdata", project, "inits_nchains.qs"))) {
    cat("\n", "Already done", project, "\n")
  } else {
    cat("\n", "Doing", project, "\n")
    
    source(file.path("models", project, "model_code.R"))
    
    inits_chains <- lapply(seq_len(n_chains), function(c) {
      qread(file.path("data", "realdata", project, paste0("inits_chain", c, ".qs")))
    })
    
    inits_nchains <- lapply(inits_chains, function(inits_chain) {
      model <- nimbleModel(
        code = model_code,
        name = 'model_nimble',
        constants = const,
        data = data,
        inits = inits_chain
      )
      
      sim_nodes <- model$getDependencies(
        names(inits_chain),
        self = FALSE,
        downstream = TRUE
      )
      sim_nodes <- as.matrix(sim_nodes)
      model$simulate(nodes = sim_nodes)
      
      sim_var <- model$getVarNames()
      inits <- lapply(sim_var, function(var_name)
        model[[var_name]])
      names(inits) <- sim_var
      
      return(inits)
    })
    
    inits_1chain <- inits_nchains[[1]]
    
    qsave(inits_nchains, file = file.path("data", "realdata", project, "inits_nchains.qs"))
    qsave(inits_1chain, file = file.path("data", "realdata", project, "inits_1chain.qs"))
    cat("\n", "Done", project, "\n")
  }
}

## -----------------------------------------------------------------------------
## 4. Define monitors
## -----------------------------------------------------------------------------

## Monitor 1
build_monitor1 <- function(mualpha1, sdalpha1, beta1, mualpha2, sdalpha2, beta2) {
  c(
    monitor1_params <- c(
      "Mu_L1", "Mu_L1surv", "Mu_L2", "Mu_L2m", "Mu_L2nm",
      "Sd_L1", "Sd_L1surv", "Sd_L2", "Sd_L2m", "Sd_L2nm",
      "Mu_Mu_L1", "Mu_Mu_L1surv", "Mu_Mu_L2", "Mu_Mu_L2m", "Mu_Mu_L2nm",
      "Mu_Sd_L1", "Mu_Sd_L1surv", "Mu_Sd_L2", "Mu_Sd_L2m", "Mu_Sd_L2nm",
      "Mu_N1", "Mu_N2", "Mu_N2m", "Mu_N2nm", "Mu_N3", "Mu_N4",
      "Gamma1_F", "Gamma1", "Gamma3", "Gamma4"),
    "Alpha1",
    if (mualpha1) "Mu_Alpha1",
    if (sdalpha1) "Sd_Alpha1",
    if (beta1) "Beta1",
    "Alpha2",
    if (mualpha2) "Mu_Alpha2",
    if (sdalpha2) "Sd_Alpha2",
    if (beta2) "Beta2"
  )
}

model_monitor1 <- tibble::tribble(
  ~model, ~mualpha1, ~sdalpha1, ~beta1, ~mualpha2, ~sdalpha2, ~beta2,
  "M0",   TRUE,   TRUE,   TRUE,   TRUE,   TRUE,   TRUE,
  "M1",   TRUE,   TRUE,   TRUE,   TRUE,   TRUE,   TRUE,
  "M2",   TRUE,   TRUE,   TRUE,   TRUE,   TRUE,   TRUE,
  "M3",   TRUE,   TRUE,   TRUE,   TRUE,   TRUE,   TRUE,
  "M4",   FALSE,  FALSE,  TRUE,   TRUE,   TRUE,   TRUE,
  "M5",   TRUE,   TRUE,   TRUE,   FALSE,  FALSE,  TRUE,
  "M6",   FALSE,  FALSE,  TRUE,   FALSE,  FALSE,  TRUE,
  "M7",   TRUE,   TRUE,   FALSE,  TRUE,   TRUE,   TRUE,
  "M8",   TRUE,   TRUE,   TRUE,   TRUE,   TRUE,   FALSE,
  "M9",   TRUE,   TRUE,   FALSE,  TRUE,   TRUE,   FALSE
)

monitor1 <- purrr::pmap(
  model_monitor1 %>% select(mualpha1, sdalpha1, beta1, mualpha2, sdalpha2, beta2),
  build_monitor1
) %>% setNames(model_monitor1$model)
  
purrr::iwalk(monitor1, function(monitor1, model) {
  qs::qsave(monitor1, file = file.path(file.path("data", "realdata", model), "monitor1.qs"))
})

## Monitor 2
for (project in projects) {
  if (file.exists(file.path("data","realdata", project, "monitor2.qs"))) {
    cat("\n", "Already done", project, "\n")
  } else {
    cat("\n", "Doing", project, "\n")
    
    source(file.path("models", project, "model_code.R"))
    
    inits_1chain <- qread(file = file.path("data", "realdata", project, "inits_1chain.qs"))
    
    model <- nimbleModel(
      code = model_code,
      name = 'model_nimble',
      constants = const,
      data = data,
      inits = inits_1chain
    )
    
    monitor2 <- paste0("logProb_", model$getNodeNames(dataOnly = TRUE))
    qsave(monitor2, file = file.path("data", "realdata", project, "monitor2.qs"))
    cat("\n", "Done", project, "\n")
  }
}