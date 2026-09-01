### ============================================================================
### Disentangling length-dependent and length-independent variations
### in survival and maturation in Atlantic salmon
### ----------------------------------------------------------------------------
### File:    figA5_ppc-alpha.R
### Purpose: Builds Appendix 5.2 Table A5.2.6.
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
pkginstall(c("dplyr", "tidyr", "coda", "nimble", "qs"))

## -----------------------------------------------------------------------------
## 1. Define project settings
## -----------------------------------------------------------------------------
projects <- c("M1")

## -----------------------------------------------------------------------------
## 2. Import MCMC object and data
## -----------------------------------------------------------------------------
load_mcmc <- function(project) {
  MCMC <- qread(paste0("saves/", project, "/MCMC.qs"))
  MCMC %>%
    seq_along() %>%
    lapply(function(i) mcmc(MCMC[[i]]$samples)) %>%
    mcmc.list()
}
MCMC_samples <- load_mcmc(projects)
MCMC_matrix <- as.matrix(MCMC_samples)

data  <- qread("data/realdata/data.qs")
const <- qread("data/realdata/const.qs")
inits <- qread(file = file.path("data", "realdata", project, "inits_1chain.qs"))

source("functions/nf_l.R")
source("functions/nf_pi.R")
source("functions/nf_res.R")

source(file.path("models", project, "model_code.R"))

model_nimble <- nimbleModel(
  code = model_code,
  name = 'model_nimble',
  constants = const,
  data = data,
  inits = inits
)
compiled_model <- compileNimble(model_nimble)
theta_order  <- model_nimble$topologicallySortNodes(colnames(MCMC_matrix))
theta_order  <- model_nimble$expandNodeNames(theta_order, returnScalarComponents = TRUE)
MCMC_matrix2 <- MCMC_matrix[, theta_order]

## -----------------------------------------------------------------------------
## 3. PPC test 1: chi-squared discrepancies
## -----------------------------------------------------------------------------
nf_ppc1_alpha <- nimbleFunction(
  
  setup = function(model, samples, alphaNode, muNode, sdNode) {
    
    theta <- colnames(samples)
    theta <- model$topologicallySortNodes(theta)
    theta <- model$expandNodeNames(theta, returnScalarComponents = TRUE)
    deps  <- model$getDependencies(theta, self = FALSE)
    
    alphaNodesExpanded <- model$expandNodeNames(alphaNode, returnScalarComponents = TRUE)
    nT <- length(alphaNodesExpanded)
    
    muNodeExpanded <- model$expandNodeNames(muNode, returnScalarComponents = TRUE)
    sdNodeExpanded <- model$expandNodeNames(sdNode, returnScalarComponents = TRUE)
    
    T_obs_store <- matrix(0, nrow = 1, ncol = 1)
    T_rep_store <- matrix(0, nrow = 1, ncol = 1)
  },
  
  run = function(samplesSub = double(2)) {
    
    nSamp <- dim(samplesSub)[1]
    setSize(T_obs_store, nSamp, 1)
    setSize(T_rep_store, nSamp, 1)
    
    for (i in 1:nSamp) {
      
      values(model, theta) <<- samplesSub[i, ]
      model$calculate(deps)
      
      muFitVec <- values(model, muNodeExpanded)
      sdFitVec <- values(model, sdNodeExpanded)
      muFit <- muFitVec[1]
      sdFit <- sdFitVec[1]
      
      alphaFit <- values(model, alphaNodesExpanded)
      
      alphaRep <- numeric(nT)
      for (t in 1:nT) alphaRep[t] <- rnorm(1, muFit, sdFit)
      
      chiFit <- 0
      chiRep <- 0
      for (t in 1:nT) {
        chiFit <- chiFit + (alphaFit[t] - muFit)^2 / sdFit^2
        chiRep <- chiRep + (alphaRep[t] - muFit)^2 / sdFit^2
      }
      
      T_obs_store[i, 1] <<- chiFit
      T_rep_store[i, 1] <<- chiRep
    }
  },
  
  methods = list(
    getTobs = function() { returnType(double(2)); return(T_obs_store) },
    getTrep = function() { returnType(double(2)); return(T_rep_store) }
  )
)

## -----------------------------------------------------------------------------
## 4. PPC test 2: lag-1 autocorrelation (AR1-type) statistic
## -----------------------------------------------------------------------------
nf_ppc2_alpha <- nimbleFunction(
  
  setup = function(model, samples, alphaNode, muNode, sdNode) {
    
    theta <- colnames(samples)
    theta <- model$topologicallySortNodes(theta)
    theta <- model$expandNodeNames(theta, returnScalarComponents = TRUE)
    deps  <- model$getDependencies(theta, self = FALSE)
    
    alphaNodesExpanded <- model$expandNodeNames(alphaNode, returnScalarComponents = TRUE)
    nT <- length(alphaNodesExpanded)
    
    muNodeExpanded <- model$expandNodeNames(muNode, returnScalarComponents = TRUE)
    sdNodeExpanded <- model$expandNodeNames(sdNode, returnScalarComponents = TRUE)
    
    T_obs_store <- matrix(0, nrow = 1, ncol = 1)
    T_rep_store <- matrix(0, nrow = 1, ncol = 1)
  },
  
  run = function(samplesSub = double(2)) {
    
    nSamp <- dim(samplesSub)[1]
    setSize(T_obs_store, nSamp, 1)
    setSize(T_rep_store, nSamp, 1)
    
    for (i in 1:nSamp) {
      
      values(model, theta) <<- samplesSub[i, ]
      model$calculate(deps)
      
      muFitVec <- values(model, muNodeExpanded)
      sdFitVec <- values(model, sdNodeExpanded)
      muFit <- muFitVec[1]
      sdFit <- sdFitVec[1]
      
      alphaFit <- values(model, alphaNodesExpanded)
      
      alphaRep <- numeric(nT)
      for (t in 1:nT) alphaRep[t] <- rnorm(1, muFit, sdFit)
      
      ar1Fit <- 0
      ar1Rep <- 0
      for (t in 1:(nT - 1)) {
        ar1Fit <- ar1Fit + (alphaFit[t] - muFit) * (alphaFit[t + 1] - muFit) / sdFit^2
        ar1Rep <- ar1Rep + (alphaRep[t] - muFit) * (alphaRep[t + 1] - muFit) / sdFit^2
      }
      
      T_obs_store[i, 1] <<- ar1Fit
      T_rep_store[i, 1] <<- ar1Rep
    }
  },
  
  methods = list(
    getTobs = function() { returnType(double(2)); return(T_obs_store) },
    getTrep = function() { returnType(double(2)); return(T_rep_store) }
  )
)

## -----------------------------------------------------------------------------
## 4. Bayesian p-values for test 1
## -----------------------------------------------------------------------------
compute_pB <- function(T_obs, T_rep) {
  mean(T_rep >= T_obs)
}

# Survival
ppc_R1 <- nf_ppc1_alpha(model_nimble, MCMC_matrix2, 'Alpha1', 'Mu_Alpha1', 'Sd_Alpha1')
ppc_C1 <- compileNimble(ppc_R1, project = compiled_model)

ppc_C1$run(MCMC_matrix2)
T_obs1 <- ppc_C1$getTobs()
T_rep1 <- ppc_C1$getTrep()
res1 <- compute_pB(T_obs1, T_rep1)

# Maturation female
ppc_R2f <- nf_ppc1_alpha(model_nimble, MCMC_matrix2, 'Alpha2[,1]', 'Mu_Alpha2[1]', 'Sd_Alpha2[1]')
ppc_C2f <- compileNimble(ppc_R2f, project = compiled_model)

ppc_C2f$run(MCMC_matrix2)
T_obs2f <- ppc_C2f$getTobs()
T_rep2f <- ppc_C2f$getTrep()
res2f <- compute_pB(T_obs2f, T_rep2f)

# Maturation male
ppc_R2m <- nf_ppc1_alpha(model_nimble, MCMC_matrix2, 'Alpha2[,2]', 'Mu_Alpha2[2]', 'Sd_Alpha2[2]')
ppc_C2m <- compileNimble(ppc_R2m, project = compiled_model)

ppc_C2m$run(MCMC_matrix2)
T_obs2m <- ppc_C2m$getTobs()
T_rep2m <- ppc_C2m$getTrep()
res2m <- compute_pB(T_obs2m, T_rep2m)

# Dataframe
pB_df1 <- data.frame(
  Survival  = res1,
  `Maturation female` = res2f,
  `Maturation male` = res2m,
  check.names = FALSE
)
print(pB_df1, row.names = FALSE)
write.csv2(pB_df1, file.path("results",project,"ppc_alpha1.csv"))

## -----------------------------------------------------------------------------
## 5. Bayesian p-values for test 2
## -----------------------------------------------------------------------------
compute_pB <- function(T_obs, T_rep) {
  mean(T_rep >= T_obs)
}

# Survival
ppc_R1 <- nf_ppc2_alpha(model_nimble, MCMC_matrix2, 'Alpha1', 'Mu_Alpha1', 'Sd_Alpha1')
ppc_C1 <- compileNimble(ppc_R1, project = compiled_model)

ppc_C1$run(MCMC_matrix2)
T_obs1 <- ppc_C1$getTobs()
T_rep1 <- ppc_C1$getTrep()
res1 <- compute_pB(T_obs1, T_rep1)

# Maturation female
ppc_R2f <- nf_ppc2_alpha(model_nimble, MCMC_matrix2, 'Alpha2[,1]', 'Mu_Alpha2[1]', 'Sd_Alpha2[1]')
ppc_C2f <- compileNimble(ppc_R2f, project = compiled_model)

ppc_C2f$run(MCMC_matrix2)
T_obs2f <- ppc_C2f$getTobs()
T_rep2f <- ppc_C2f$getTrep()
res2f <- compute_pB(T_obs2f, T_rep2f)

# Maturation male
ppc_R2m <- nf_ppc2_alpha(model_nimble, MCMC_matrix2, 'Alpha2[,2]', 'Mu_Alpha2[2]', 'Sd_Alpha2[2]')
ppc_C2m <- compileNimble(ppc_R2m, project = compiled_model)

ppc_C2m$run(MCMC_matrix2)
T_obs2m <- ppc_C2m$getTobs()
T_rep2m <- ppc_C2m$getTrep()
res2m <- compute_pB(T_obs2m, T_rep2m)

# Dataframe
pB_df1 <- data.frame(
  Survival  = res1,
  `Maturation female` = res2f,
  `Maturation male` = res2m,
  check.names = FALSE
)
print(pB_df1, row.names = FALSE)
write.csv2(pB_df1, file.path("results",project,"ppc_alpha2.csv"))
