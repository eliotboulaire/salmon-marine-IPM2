### ============================================================================
### Disentangling length-dependent and length-independent variations
### in survival and maturation in Atlantic salmon
### ----------------------------------------------------------------------------
### File:    run_model.R
### Purpose: Read all NIMBLE model objects from data_save.R and model_save.R
###          with model_code.R to run the NIMBLE models.
###          (1) "MCMC" output are saved in .qs to saves/[M0-M9]/
###          (2) run_info are saved in .txt to saves/[M0-M9]/
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
## 2. Import model, data and const
## -----------------------------------------------------------------------------
# Data and const
data <- qread("data/realdata/data.qs")
const <- qread("data/realdata/const.qs")

# Model functions
source("functions/nf_res.R")
source("functions/nf_l.R")
source("functions/nf_pi.R")

## -----------------------------------------------------------------------------
## 3. Define MCMC settings
## -----------------------------------------------------------------------------
n_thin <- 300
n_keep <- 30000/n_chains
n_burnin <- 500000
n_iter <- n_keep * n_thin + n_burnin

## -----------------------------------------------------------------------------
## 4. Run models
## -----------------------------------------------------------------------------
for (project in projects) {
  inits_nchains <- qread(file = file.path("data", "realdata", project, "inits_nchains.qs"))
  inits_1chain <- qread(file = file.path("data", "realdata", project, "inits_1chain.qs"))
  source(file.path("models", project, "model_code.R"))
  monitor1 <- qread(file = file.path("data", "realdata", project, "monitor1.qs"))
  monitor2 <- qread(file = file.path("data", "realdata", project, "monitor2.qs"))
  
  cl <- makeCluster(n_chains)
  on.exit({
    clusterCall(cl, function() {
      rm(list = ls())
      gc(full = TRUE)
    })
    
    stopCluster(cl)
    gc(full = TRUE)
  })
  
  clusterExport(cl, c("projects", "seeds",
                      "data", "const", "inits_nchains", "inits_1chain",
                      "nf_res", "nf_l", "nf_pi", "nf1", "nf2",
                      "model_code", "monitor1", "monitor2",
                      "n_chains", "n_thin", "n_burnin", "n_iter"), envir = environment())
  
  clusterEvalQ(cl, {
    pkginstall <- function(packages) {
      for (pkg in packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
          message("Installing: ", pkg)
          install.packages(pkg, dependencies = TRUE)
        }
        library(pkg, character.only = TRUE)
        message(pkg, " loaded")
      }
    }
    
    packages <- c("nimble", "parallel", "coda", "qs")
    pkginstall(packages)
    
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
  })
  
  start_time <- Sys.time()
  MCMC <- parLapply(cl, seq_len(n_chains), function(i) {
    runMCMC(
      mcmc = compiled_MCMC,
      niter = n_iter,
      nburnin = n_burnin,
      nchains = 1,
      thin = n_thin,
      thin2 = n_thin,
      inits = inits_nchains[[i]],
      progressBar = FALSE,
      samples = TRUE,
      samplesAsCodaMCMC = TRUE,
      setSeed = seeds[i],
      summary = FALSE,
      WAIC = FALSE
    )
  })
  end_time <- Sys.time()
  
  dir.create(file.path(path, "saves", project), recursive = TRUE, showWarnings = FALSE)
  qsave(MCMC, file = file.path(path, "saves", project, "MCMC.qs"))
  
  run_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  writeLines(c(paste(project),
               paste("seeds", seeds),
               paste("run info", n_chains, n_iter, n_burnin, n_thin),
               paste("run time",  run_time, "secs")),
             file.path("saves", project, "run_info.txt"))
  
  gc(full = TRUE)
}