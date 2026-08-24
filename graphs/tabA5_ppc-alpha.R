### ============================================================================
### Creating table of Posterior Predictive Check on alpha estimate assumption
### ============================================================================
rm(list = ls())
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
pkginstall(c("coda", "nimble", "ggplot2", "qs", "dplyr", "tidyr"))
set.seed(123)
project <- paste0("M", 1)

## =============================================================================
## 1. Import model, data and const
## =============================================================================
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

## =============================================================================
## 2. Import posterior MCMC samples and reorder to match NIMBLE's internal node ordering
## =============================================================================
load_mcmc <- function(model) {
  MCMC <- qread(file.path("saves", model, "MCMC.qs"))
  MCMC %>%
    seq_along() %>%
    lapply(function(i) mcmc(MCMC[[i]]$samples)) %>%
    mcmc.list()
}
MCMC_matrix <- as.matrix(load_mcmc(project))

theta_order  <- model_nimble$topologicallySortNodes(colnames(MCMC_matrix))
theta_order  <- model_nimble$expandNodeNames(theta_order, returnScalarComponents = TRUE)
MCMC_matrix2 <- MCMC_matrix[, theta_order]

## =============================================================================
## 3. NIMBLE function: compute chi-squared discrepancies T_obs and T_rep
##    for one data source (e.g. LogN1, LogN3, or LogN4), for every posterior draw
## =============================================================================
## =============================================================================
## PPC test 1 -- chi2 on fitted year-effects (alpha) vs their own fitted hyper-Normal
## =============================================================================
nf_ppc1_alpha <- nimbleFunction(
  
  setup = function(model, samples, alphaNode, muNode, sdNode) {
    
    theta <- colnames(samples)
    theta <- model$topologicallySortNodes(theta)
    theta <- model$expandNodeNames(theta, returnScalarComponents = TRUE)
    deps  <- model$getDependencies(theta, self = FALSE)
    
    alphaNodesExpanded <- model$expandNodeNames(alphaNode, returnScalarComponents = TRUE)
    nT <- length(alphaNodesExpanded)   # number of years, e.g. 24
    
    muNodeExpanded <- model$expandNodeNames(muNode, returnScalarComponents = TRUE)
    sdNodeExpanded <- model$expandNodeNames(sdNode, returnScalarComponents = TRUE)
    
    T_obs_store <- matrix(0, nrow = 1, ncol = 1)   # one column: pooled T per iteration
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
      
      alphaFit <- values(model, alphaNodesExpanded)   # length nT
      
      # --- draw ONE replicate year-effect vector from the same fitted hyper-Normal ---
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

## =============================================================================
## PPC test 2 -- lag-1 autocorrelation (AR1-type) statistic on fitted year-effects
## =============================================================================
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

## =============================================================================
## 4. Build, compile, and run for each stage: smolts (1), 1SW (3), 2SW (4)
## =============================================================================

## =============================================================================
## PPC test 1 -- 
## =============================================================================

## ---- Survival ----
ppc_R1 <- nf_ppc1_alpha(model_nimble, MCMC_matrix2, 'Alpha1', 'Mu_Alpha1', 'Sd_Alpha1')
ppc_C1 <- compileNimble(ppc_R1, project = compiled_model)

ppc_C1$run(MCMC_matrix2)
T_obs1 <- ppc_C1$getTobs()
T_rep1 <- ppc_C1$getTrep()

## ---- Maturation female ----
ppc_R2f <- nf_ppc1_alpha(model_nimble, MCMC_matrix2, 'Alpha2[,1]', 'Mu_Alpha2[1]', 'Sd_Alpha2[1]')
ppc_C2f <- compileNimble(ppc_R2f, project = compiled_model)

ppc_C2f$run(MCMC_matrix2)
T_obs2f <- ppc_C2f$getTobs()
T_rep2f <- ppc_C2f$getTrep()

## ---- Maturation male ----
ppc_R2m <- nf_ppc1_alpha(model_nimble, MCMC_matrix2, 'Alpha2[,2]', 'Mu_Alpha2[2]', 'Sd_Alpha2[2]')
ppc_C2m <- compileNimble(ppc_R2m, project = compiled_model)

ppc_C2m$run(MCMC_matrix2)
T_obs2m <- ppc_C2m$getTobs()
T_rep2m <- ppc_C2m$getTrep()

compute_pB <- function(T_obs, T_rep) {
  mean(T_rep >= T_obs)
}

res1 <- compute_pB(T_obs1, T_rep1)
res2f <- compute_pB(T_obs2f, T_rep2f)
res2m <- compute_pB(T_obs2m, T_rep2m)

pB_df1 <- data.frame(
  Survival  = res1,
  `Maturation female` = res2f,
  `Maturation male` = res2m,
  check.names = FALSE
)
print(pB_df1, row.names = FALSE)
write.csv2(pB_df1, file.path("results",project,"ppc_alpha1.csv"))

pB_df2 <- pB_df1 %>%
  pivot_longer(cols = everything(), names_to = "stage", values_to = "pB") %>%
  mutate(
    stage = factor(stage, levels = c("Survival", "Maturation female", "Maturation male")),
    dev   = pB - 0.5
  )
ggplot(pB_df2) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.05, ymax = 0.95,
           alpha = 0.2, fill = "grey", colour = "black") +
  geom_hline(yintercept = 0.5, linewidth = 1.1) +
  geom_point(aes(x = stage, y = pB), colour = "black", size = 5) +
  geom_point(aes(x = stage, y = pB, colour = dev), size = 3.8) +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-0.5, 0.5),
                         name = "Deviation\nfrom 0.5") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Stage", y = "Bayesian p-value") +
  theme_minimal()

## =============================================================================
## PPC test 2 -- 
## =============================================================================

## ---- Survival ----
ppc_R1 <- nf_ppc2_alpha(model_nimble, MCMC_matrix2, 'Alpha1', 'Mu_Alpha1', 'Sd_Alpha1')
ppc_C1 <- compileNimble(ppc_R1, project = compiled_model)

ppc_C1$run(MCMC_matrix2)
T_obs1 <- ppc_C1$getTobs()
T_rep1 <- ppc_C1$getTrep()

## ---- Maturation female ----
ppc_R2f <- nf_ppc2_alpha(model_nimble, MCMC_matrix2, 'Alpha2[,1]', 'Mu_Alpha2[1]', 'Sd_Alpha2[1]')
ppc_C2f <- compileNimble(ppc_R2f, project = compiled_model)

ppc_C2f$run(MCMC_matrix2)
T_obs2f <- ppc_C2f$getTobs()
T_rep2f <- ppc_C2f$getTrep()

## ---- Maturation male ----
ppc_R2m <- nf_ppc2_alpha(model_nimble, MCMC_matrix2, 'Alpha2[,2]', 'Mu_Alpha2[2]', 'Sd_Alpha2[2]')
ppc_C2m <- compileNimble(ppc_R2m, project = compiled_model)

ppc_C2m$run(MCMC_matrix2)
T_obs2m <- ppc_C2m$getTobs()
T_rep2m <- ppc_C2m$getTrep()

compute_pB <- function(T_obs, T_rep) {
  mean(T_rep >= T_obs)
}

res1 <- compute_pB(T_obs1, T_rep1)
res2f <- compute_pB(T_obs2f, T_rep2f)
res2m <- compute_pB(T_obs2m, T_rep2m)

pB_df1 <- data.frame(
  Survival  = res1,
  `Maturation female` = res2f,
  `Maturation male` = res2m,
  check.names = FALSE
)
print(pB_df1, row.names = FALSE)
write.csv2(pB_df1, file.path("results",project,"ppc_alpha2.csv"))

pB_df2 <- pB_df1 %>%
  pivot_longer(cols = everything(), names_to = "stage", values_to = "pB") %>%
  mutate(
    stage = factor(stage, levels = c("Survival", "Maturation female", "Maturation male")),
    dev   = pB - 0.5
  )
ggplot(pB_df2) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.05, ymax = 0.95,
           alpha = 0.2, fill = "grey", colour = "black") +
  geom_hline(yintercept = 0.5, linewidth = 1.1) +
  geom_point(aes(x = stage, y = pB), colour = "black", size = 5) +
  geom_point(aes(x = stage, y = pB, colour = dev), size = 3.8) +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-0.5, 0.5),
                         name = "Deviation\nfrom 0.5") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Stage", y = "Bayesian p-value") +
  theme_minimal()