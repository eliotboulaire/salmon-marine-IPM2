### ============================================================================
### Disentangling length-dependent and length-independent variations
### in survival and maturation in Atlantic salmon
### ----------------------------------------------------------------------------
### File:    debug_model.R
### Purpose: Read all NIMBLE model objects from data_save.R and model_save.R
###          with model_code.R to run a specific  NIMBLE models.
###          Allows to check where are the problems encountered in run_model.R
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
project <- paste0("M", 0)
n_chains <- 6
seeds <- sample(1:1e9, n_chains)

## -----------------------------------------------------------------------------
## 2. Import model, data and const
## -----------------------------------------------------------------------------
# Data and const
data <- qread("data/realdata/data.qs")
const <- qread("data/realdata/const.qs")

# Model
source(file.path("models", project, "model_code.R"))

# Model functions
source("functions/nf_res.R")
source("functions/nf_l.R")
source("functions/nf_pi.R")

## -----------------------------------------------------------------------------
## 3. Define inits and monitors
## -----------------------------------------------------------------------------
# Model functions
inits_nchains <- qread(file.path("data","realdata", project, "inits_nchains.qs"))
inits_1chain <- qread(file.path("data","realdata", project, "inits_1chain.qs"))

# --- Define monitor2 ---
manitor1 <- qread(file.path("data","realdata", project, "monitor1.qs"))
manitor2 <- qread(file.path("data","realdata", project, "monitor2.qs"))

## -----------------------------------------------------------------------------
## 4. Define MCMC settings
## -----------------------------------------------------------------------------
n_thin <- 1
n_keep <- 1000
n_burnin <- 100
n_iter <- n_keep * n_thin + n_burnin

## -----------------------------------------------------------------------------
## 5. Run model
## -----------------------------------------------------------------------------
model_nimble <- nimbleModel(
  code = model_code,
  name = 'model_nimble',
  constants = const,
  data = data,
  inits = inits_1chain
)
compiled_model <- compileNimble(model_nimble)
    
model_conf <- configureMCMC(
  model = model_nimble,
  thin = 1,
  monitors = monitor1,
  thin2 = 1,
  monitors2 = monitor2,
  inits = inits_1chain,
  enableWAIC = FALSE
)
model_MCMC <- buildMCMC(model_conf)
compiled_MCMC <- compileNimble(model_MCMC, project = model_nimble)
  
MCMC <- runMCMC(
  mcmc = compiled_MCMC,
  niter = n_iter,
  nburnin = n_burnin,
  nchains = 1,
  thin = n_thin,
  thin2 = n_thin,
  inits = inits_nchains,
  progressBar = FALSE,
  samples = TRUE,
  samplesAsCodaMCMC = TRUE,
  setSeed = seeds,
  summary = FALSE,
  WAIC = FALSE
)
