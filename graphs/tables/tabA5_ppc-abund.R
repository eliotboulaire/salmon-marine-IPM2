### ============================================================================
### Disentangling length-dependent and length-independent variations
### in survival and maturation in Atlantic salmon
### ----------------------------------------------------------------------------
### File:    figA5_ppc-abund.R
### Purpose: Builds Appendix 5.2 Table A5.2.1.
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
## 3. PPC test: chi-squared discrepancies
## -----------------------------------------------------------------------------
nf_ppc_abund <- nimbleFunction(
  
  setup = function(model, samples, obsNodes, muNodes, sdFixed) {
    
    theta <- colnames(samples)
    theta <- model$topologicallySortNodes(theta)
    theta <- model$expandNodeNames(theta, returnScalarComponents = TRUE)
    deps  <- model$getDependencies(theta, self = FALSE)
    
    obsNodesExpanded <- model$expandNodeNames(obsNodes, returnScalarComponents = TRUE)
    muNodesExpanded  <- model$expandNodeNames(muNodes,  returnScalarComponents = TRUE)
    n <- length(obsNodesExpanded)
    
    obsData <- log(values(model, obsNodesExpanded))
    sdFixed <- as.numeric(sdFixed)
    
    T_obs_store <- matrix(0, nrow = 1, ncol = 1)
    T_rep_store <- matrix(0, nrow = 1, ncol = 1)
  },
  
  run = function(samplesSub = double(2)) {
    
    nSamp <- dim(samplesSub)[1]
    
    setSize(T_obs_store, nSamp, n)
    setSize(T_rep_store, nSamp, n)
    
    for (i in 1:nSamp) {
      
      values(model, theta) <<- samplesSub[i, ]
      model$calculate(deps)
      muFit <- values(model, muNodesExpanded)
      
      for (j in 1:n) {T_obs_store[i, j] <<- (obsData[j] - muFit[j])^2 / sdFixed[j]^2}
      
      model$simulate(obsNodes, includeData = TRUE)
      repVals <- log(values(model, obsNodesExpanded))
      
      for (j in 1:n) {T_rep_store[i, j] <<- (repVals[j] - muFit[j])^2 / sdFixed[j]^2}
    }
  },
  
  methods = list(
    getTobs = function() {
      returnType(double(2))
      return(T_obs_store)
    },
    getTrep = function() {
      returnType(double(2))
      return(T_rep_store)
    }
  )
)

## -----------------------------------------------------------------------------
## 4. Bayesian p-values
## -----------------------------------------------------------------------------
compute_pB <- function(T_obs, T_rep) {
  list(
    pB_by_year = colMeans(T_rep >= T_obs),
    pB_pooled  = mean(rowSums(T_rep) >= rowSums(T_obs))
  )
}

# Smolts
ppc_R1 <- nf_ppc_abund(model_nimble, MCMC_matrix2, 'LogN1', 'Mu_LogN1', const$Sd_LogN1)
ppc_C1 <- compileNimble(ppc_R1, project = compiled_model)

ppc_C1$run(MCMC_matrix2)
T_obs1 <- ppc_C1$getTobs()
T_rep1 <- ppc_C1$getTrep()
res1 <- compute_pB(T_obs1, T_rep1)

# 1SW
ppc_R3 <- nf_ppc_abund(model_nimble, MCMC_matrix2, 'LogN3', 'Mu_LogN3', const$Sd_LogN3)
ppc_C3 <- compileNimble(ppc_R3, project = compiled_model)

ppc_C3$run(MCMC_matrix2)
T_obs3 <- ppc_C3$getTobs()
T_rep3 <- ppc_C3$getTrep()
res3 <- compute_pB(T_obs3, T_rep3)

# 2SW
ppc_R4 <- nf_ppc_abund(model_nimble, MCMC_matrix2, 'LogN4', 'Mu_LogN4', const$Sd_LogN4)
ppc_C4 <- compileNimble(ppc_R4, project = compiled_model)

ppc_C4$run(MCMC_matrix2)
T_obs4 <- ppc_C4$getTobs()
T_rep4 <- ppc_C4$getTrep()
res4 <- compute_pB(T_obs4, T_rep4)

# Dataframe
pB_df1 <- data.frame(
  year   = 1996:2019,
  Smolts = res1$pB_by_year,
  `1SW`  = res3$pB_by_year,
  `2SW`  = res4$pB_by_year,
  check.names = FALSE
)
print(pB_df1, row.names = FALSE)
write.csv2(pB_df1, file.path("results",project,"ppc_abundance.csv"))
