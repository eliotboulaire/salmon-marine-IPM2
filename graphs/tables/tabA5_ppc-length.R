### ============================================================================
### Disentangling length-dependent and length-independent variations
### in survival and maturation in Atlantic salmon
### ----------------------------------------------------------------------------
### File:    figA5_ppc-length.R
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
## 3. PPC test 1: chi-squared discrepancies of mean scale lengths
## -----------------------------------------------------------------------------
nf_ppc1_scale <- nimbleFunction(
  
  setup = function(model, samples, muNode, sdNode, obsMat, nIndivPerCohort) {
    
    theta <- colnames(samples)
    theta <- model$topologicallySortNodes(theta)
    theta <- model$expandNodeNames(theta, returnScalarComponents = TRUE)
    deps  <- model$getDependencies(theta, self = FALSE)
    
    muNodesExpanded <- model$expandNodeNames(muNode, returnScalarComponents = TRUE)
    sdNodesExpanded <- model$expandNodeNames(sdNode, returnScalarComponents = TRUE)
    nC <- length(muNodesExpanded)
    
    nIndivPerCohort <- as.integer(nIndivPerCohort)
    
    obsMat <- as.matrix(obsMat)
    
    obsMean <- numeric(nC)
    for (y in 1:nC) {
      n <- nIndivPerCohort[y]
      x <- obsMat[1:n, y]
      obsMean[y] <- mean(x)
    }
    
    T_obs_store <- matrix(0, nrow = 1, ncol = 1)
    T_rep_store <- matrix(0, nrow = 1, ncol = 1)
  },
  
  run = function(samplesSub = double(2)) {
    
    nSamp <- dim(samplesSub)[1]
    setSize(T_obs_store, nSamp, nC)
    setSize(T_rep_store, nSamp, nC)
    
    for (i in 1:nSamp) {
      
      values(model, theta) <<- samplesSub[i, ]
      model$calculate(deps)
      
      for (y in 1:nC) {
        
        muFitVec <- values(model, muNodesExpanded[y])
        sdFitVec <- values(model, sdNodesExpanded[y])
        muFit <- muFitVec[1]
        sdFit <- sdFitVec[1]
        
        n <- nIndivPerCohort[y]
        
        repVals <- numeric(n)
        for (j in 1:n) repVals[j] <- rnorm(1, muFit, sdFit)
        meanRep <- mean(repVals[1:n])
        
        T_obs_store[i, y] <<- (obsMean[y] - muFit)^2 / sdFit^2
        T_rep_store[i, y] <<- (meanRep    - muFit)^2 / sdFit^2
      }
    }
  },
  
  methods = list(
    getTobs = function() { returnType(double(2)); return(T_obs_store) },
    getTrep = function() { returnType(double(2)); return(T_rep_store) }
  )
)

## -----------------------------------------------------------------------------
## 4. PPC test 2: chi-squared discrepancies of 90% IQR
## -----------------------------------------------------------------------------
nf_ppc2_scale <- nimbleFunction(
  
  setup = function(model, samples, muNode, sdNode, obsMat, nIndivPerCohort, alpha) {
    
    theta <- colnames(samples)
    theta <- model$topologicallySortNodes(theta)
    theta <- model$expandNodeNames(theta, returnScalarComponents = TRUE)
    deps  <- model$getDependencies(theta, self = FALSE)
    
    muNodesExpanded <- model$expandNodeNames(muNode, returnScalarComponents = TRUE)
    sdNodesExpanded <- model$expandNodeNames(sdNode, returnScalarComponents = TRUE)
    nC <- length(muNodesExpanded)
    
    nIndivPerCohort <- as.integer(nIndivPerCohort)
    alpha <- as.numeric(alpha)
    
    obsMat <- as.matrix(obsMat)
    
    obsIQR <- numeric(nC)
    for (y in 1:nC) {
      n <- nIndivPerCohort[y]
      x <- obsMat[1:n, y]
      q <- quantile(x, probs = c(alpha, 1 - alpha), type = 7, names = FALSE)
      obsIQR[y] <- q[2] - q[1]
    }
    
    T_obs_store <- matrix(0, nrow = 1, ncol = 1)
    T_rep_store <- matrix(0, nrow = 1, ncol = 1)
  },
  
  run = function(samplesSub = double(2)) {
    
    nSamp <- dim(samplesSub)[1]
    setSize(T_obs_store, nSamp, nC)
    setSize(T_rep_store, nSamp, nC)
    
    for (i in 1:nSamp) {
      
      values(model, theta) <<- samplesSub[i, ]
      model$calculate(deps)
      
      for (y in 1:nC) {
        
        muFitVec <- values(model, muNodesExpanded[y])
        sdFitVec <- values(model, sdNodesExpanded[y])
        muFit <- muFitVec[1]
        sdFit <- sdFitVec[1]
        
        n <- nIndivPerCohort[y]
        repVals <- numeric(n)
        for (j in 1:n) repVals[j] <- rnorm(1, muFit, sdFit)
        
        iqrFit <- sdFit * (qnorm(1 - alpha, 0, 1) - qnorm(alpha, 0, 1))
        repSorted <- insertionSort(repVals[1:n], n)
        qLoRep <- quantileInterp(repSorted, n, alpha)
        qHiRep <- quantileInterp(repSorted, n, 1 - alpha)
        iqrRep <- qHiRep - qLoRep
        
        T_obs_store[i, y] <<- (obsIQR[y] - iqrFit)^2 / sdFit^2
        T_rep_store[i, y] <<- (iqrRep    - iqrFit)^2 / sdFit^2
      }
    }
  },
  
  methods = list(
    
    insertionSort = function(x = double(1), n = integer(0)) {
      xs <- x
      for (k in 2:n) {
        key <- xs[k]
        j <- k - 1
        while (j >= 1 & xs[j] > key) {
          xs[j + 1] <- xs[j]
          j <- j - 1
        }
        xs[j + 1] <- key
      }
      returnType(double(1))
      return(xs)
    },
    
    quantileInterp = function(sortedX = double(1), n = integer(0), p = double(0)) {
      h <- (n - 1) * p + 1
      lo <- floor(h)
      hi <- ceiling(h)
      if (lo < 1) lo <- 1
      if (hi > n) hi <- n
      frac <- h - lo
      val <- sortedX[lo] + frac * (sortedX[hi] - sortedX[lo])
      returnType(double(0))
      return(val)
    },
    
    getTobs = function() { returnType(double(2)); return(T_obs_store) },
    getTrep = function() { returnType(double(2)); return(T_rep_store) }
  )
)

## -----------------------------------------------------------------------------
## 5. PPC test 3: chi-squared discrepancies of length-class proportions
## -----------------------------------------------------------------------------
nf_ppc3_scale <- nimbleFunction(
  
  setup = function(model, samples, muNode, sdNode, obsMat, nIndivPerCohort,
                   breaks, midpoints) {
    
    theta <- colnames(samples)
    theta <- model$topologicallySortNodes(theta)
    theta <- model$expandNodeNames(theta, returnScalarComponents = TRUE)
    deps  <- model$getDependencies(theta, self = FALSE)
    
    muNodesExpanded <- model$expandNodeNames(muNode, returnScalarComponents = TRUE)
    sdNodesExpanded <- model$expandNodeNames(sdNode, returnScalarComponents = TRUE)
    nC <- length(muNodesExpanded)
    
    nIndivPerCohort <- as.integer(nIndivPerCohort)
    
    breaks    <- as.numeric(breaks)
    midpoints <- as.numeric(midpoints)
    nbins <- length(midpoints)
    
    obsMat <- as.matrix(obsMat)
    
    obsProp <- matrix(0, nrow = nC, ncol = nbins)
    for (y in 1:nC) {
      n <- nIndivPerCohort[y]
      x <- obsMat[1:n, y]
      cutBins <- cut(x, breaks = breaks, right = FALSE, include.lowest = TRUE)
      obsProp[y, ] <- as.numeric(table(cutBins)) / n
    }
    
    T_obs_store <- matrix(0, nrow = 1, ncol = 1)
    T_rep_store <- matrix(0, nrow = 1, ncol = 1)
  },
  
  run = function(samplesSub = double(2)) {
    
    nSamp <- dim(samplesSub)[1]
    setSize(T_obs_store, nSamp, nC)
    setSize(T_rep_store, nSamp, nC)
    
    for (i in 1:nSamp) {
      
      values(model, theta) <<- samplesSub[i, ]
      model$calculate(deps)
      
      for (y in 1:nC) {
        
        muFitVec <- values(model, muNodesExpanded[y])
        sdFitVec <- values(model, sdNodesExpanded[y])
        muFit <- muFitVec[1]
        sdFit <- sdFitVec[1]
        n <- nIndivPerCohort[y]
        
        repVals <- numeric(n)
        for (j in 1:n) repVals[j] <- rnorm(1, muFit, sdFit)
        
        pFit <- numeric(nbins)
        for (k in 1:nbins) pFit[k] <- pnorm(breaks[k+1], muFit, sdFit) - pnorm(breaks[k], muFit, sdFit)
        
        cntRep <- numeric(nbins)
        for (j in 1:n) {
          for (k in 1:nbins) {
            if (repVals[j] >= breaks[k] & repVals[j] < breaks[k+1]) cntRep[k] <- cntRep[k] + 1
          }
        }
        propRep <- cntRep / n
        
        chiO <- 0; chiR <- 0
        for (k in 1:nbins) {
          if (obsProp[y, k] > 0) {
            chiO <- chiO + (obsProp[y, k] - pFit[k])^2
          }
          if (propRep[k] > 0) {
            chiR <- chiR + (propRep[k] - pFit[k])^2
          }
        }
        T_obs_store[i, y] <<- chiO
        T_rep_store[i, y] <<- chiR
      }
    }
  },
  
  methods = list(
    getTobs = function() { returnType(double(2)); return(T_obs_store) },
    getTrep = function() { returnType(double(2)); return(T_rep_store) }
  )
)

## -----------------------------------------------------------------------------
## 6. Bayesian p-values for test 1
## -----------------------------------------------------------------------------
compute_pB <- function(T_obs, T_rep) {
  list(
    pB_by_year = colMeans(T_rep >= T_obs), 
    pB_pooled  = mean(rowSums(T_rep) >= rowSums(T_obs))
  )
}

# Smolt
ppc_R1 <- nf_ppc1_scale(model_nimble, MCMC_matrix2, 'Mu_L1', 'Sd_L1', data$S1, const$I1)
ppc_C1 <- compileNimble(ppc_R1, project = compiled_model)

ppc_C1$run(MCMC_matrix2)
T_obs1 <- ppc_C1$getTobs()
T_rep1 <- ppc_C1$getTrep()
res1 <- compute_pB(T_obs1, T_rep1)

# Surviving smolt
ppc_R1surv <- nf_ppc1_scale(model_nimble, MCMC_matrix2, 'Mu_L1surv', 'Sd_L1surv', data$S1surv, const$I1surv)
ppc_C1surv <- compileNimble(ppc_R1surv, project = compiled_model)

ppc_C1surv$run(MCMC_matrix2)
T_obs1surv <- ppc_C1surv$getTobs()
T_rep1surv <- ppc_C1surv$getTrep()
res1surv <- compute_pB(T_obs1surv, T_rep1surv)

# Post-smolt
ppc_R2 <- nf_ppc1_scale(model_nimble, MCMC_matrix2, 'Mu_L2', 'Sd_L2', data$S2, const$I2)
ppc_C2 <- compileNimble(ppc_R2, project = compiled_model)

ppc_C2$run(MCMC_matrix2)
T_obs2 <- ppc_C2$getTobs()
T_rep2 <- ppc_C2$getTrep()
res2 <- compute_pB(T_obs2, T_rep2)

# Maturing post-smolt
ppc_R2m <- nf_ppc1_scale(model_nimble, MCMC_matrix2, 'Mu_L2m', 'Sd_L2m', data$S2m, const$I2m)
ppc_C2m <- compileNimble(ppc_R2m, project = compiled_model)

ppc_C2m$run(MCMC_matrix2)
T_obs2m <- ppc_C2m$getTobs()
T_rep2m <- ppc_C2m$getTrep()
res2m <- compute_pB(T_obs2m, T_rep2m)

# Non-maturing post-smolt
ppc_R2nm <- nf_ppc1_scale(model_nimble, MCMC_matrix2, 'Mu_L2nm', 'Sd_L2nm', data$S2nm, const$I2nm)
ppc_C2nm <- compileNimble(ppc_R2nm, project = compiled_model)

ppc_C2nm$run(MCMC_matrix2)
T_obs2nm <- ppc_C2nm$getTobs()
T_rep2nm <- ppc_C2nm$getTrep()
res2nm <- compute_pB(T_obs2nm, T_rep2nm)

# Dataframe
pB_df1 <- data.frame(
  year = 1996:2019,
  Smolt  = res1$pB_by_year,
  `Surviving smolt` = res1surv$pB_by_year,
  `Post-smolt` = res2$pB_by_year,
  `Maturing post-smolt` = res2m$pB_by_year,
  `Non-maturing post-smolt` = res2nm$pB_by_year,
  check.names = FALSE
)
print(pB_df1, row.names = FALSE)
write.csv2(pB_df1, file.path("results",project,"ppc_length1.csv"))

## -----------------------------------------------------------------------------
## 7. Bayesian p-values for test 2
## -----------------------------------------------------------------------------
compute_pB <- function(T_obs, T_rep) {
  list(
    pB_by_year = colMeans(T_rep >= T_obs), 
    pB_pooled  = mean(rowSums(T_rep) >= rowSums(T_obs))
  )
}

# Smolt
ppc_R1 <- nf_ppc2_scale(model_nimble, MCMC_matrix2, 'Mu_L1', 'Sd_L1', data$S1, const$I1, 0.25)
ppc_C1 <- compileNimble(ppc_R1, project = compiled_model)

ppc_C1$run(MCMC_matrix2)
T_obs1 <- ppc_C1$getTobs()
T_rep1 <- ppc_C1$getTrep()
res1 <- compute_pB(T_obs1, T_rep1)

# Surviving smolt
ppc_R1surv <- nf_ppc2_scale(model_nimble, MCMC_matrix2, 'Mu_L1surv', 'Sd_L1surv', data$S1surv, const$I1surv, 0.25)
ppc_C1surv <- compileNimble(ppc_R1surv, project = compiled_model)

ppc_C1surv$run(MCMC_matrix2)
T_obs1surv <- ppc_C1surv$getTobs()
T_rep1surv <- ppc_C1surv$getTrep()
res1surv <- compute_pB(T_obs1surv, T_rep1surv)

# Post-smolt
ppc_R2 <- nf_ppc2_scale(model_nimble, MCMC_matrix2, 'Mu_L2', 'Sd_L2', data$S2, const$I2, 0.25)
ppc_C2 <- compileNimble(ppc_R2, project = compiled_model)

ppc_C2$run(MCMC_matrix2)
T_obs2 <- ppc_C2$getTobs()
T_rep2 <- ppc_C2$getTrep()
res2 <- compute_pB(T_obs2, T_rep2)

# Maturing post-smolt
ppc_R2m <- nf_ppc2_scale(model_nimble, MCMC_matrix2, 'Mu_L2m', 'Sd_L2m', data$S2m, const$I2m, 0.25)
ppc_C2m <- compileNimble(ppc_R2m, project = compiled_model)

ppc_C2m$run(MCMC_matrix2)
T_obs2m <- ppc_C2m$getTobs()
T_rep2m <- ppc_C2m$getTrep()
res2m <- compute_pB(T_obs2m, T_rep2m)

# Non-maturing post-smolt
ppc_R2nm <- nf_ppc2_scale(model_nimble, MCMC_matrix2, 'Mu_L2nm', 'Sd_L2nm', data$S2nm, const$I2nm, 0.25)
ppc_C2nm <- compileNimble(ppc_R2nm, project = compiled_model)

ppc_C2nm$run(MCMC_matrix2)
T_obs2nm <- ppc_C2nm$getTobs()
T_rep2nm <- ppc_C2nm$getTrep()
res2nm <- compute_pB(T_obs2nm, T_rep2nm)

# Dataframe
pB_df1 <- data.frame(
  year = 1996:2019,
  Smolt  = res1$pB_by_year,
  `Surviving smolt` = res1surv$pB_by_year,
  `Post-smolt` = res2$pB_by_year,
  `Maturing post-smolt` = res2m$pB_by_year,
  `Non-maturing post-smolt` = res2nm$pB_by_year,
  check.names = FALSE
)
print(pB_df1, row.names = FALSE)
write.csv2(pB_df1, file.path("results",project,"ppc_length2.csv"))

## -----------------------------------------------------------------------------
## 8. Bayesian p-values for test 2
## -----------------------------------------------------------------------------
compute_pB <- function(T_obs, T_rep) {
  list(
    pB_by_year = colMeans(T_rep >= T_obs), 
    pB_pooled  = mean(rowSums(T_rep) >= rowSums(T_obs))
  )
}

breaks_L1 <- seq(const$Min_L1, const$Max_L1, length.out = (const$L + 1))
midbins_L1 <- (breaks_L1[1:const$L] - breaks_L1[2:(const$L+1)])/2
breaks_L2 <- seq(const$Min_L2, const$Max_L2, length.out = (const$L + 1))
midbins_L2 <- (breaks_L2[1:const$L] - breaks_L2[2:(const$L+1)])/2

# Smolt
ppc_R1 <- nf_ppc3_scale(model_nimble, MCMC_matrix2, 'Mu_L1', 'Sd_L1', data$S1, const$I1, breaks_L1, midbins_L1)
ppc_C1 <- compileNimble(ppc_R1, project = compiled_model)

ppc_C1$run(MCMC_matrix2)
T_obs1 <- ppc_C1$getTobs()
T_rep1 <- ppc_C1$getTrep()
res1 <- compute_pB(T_obs1, T_rep1)

# Surviving smolt
ppc_R1surv <- nf_ppc3_scale(model_nimble, MCMC_matrix2, 'Mu_L1surv', 'Sd_L1surv', data$S1surv, const$I1surv, breaks_L1, midbins_L1)
ppc_C1surv <- compileNimble(ppc_R1surv, project = compiled_model)

ppc_C1surv$run(MCMC_matrix2)
T_obs1surv <- ppc_C1surv$getTobs()
T_rep1surv <- ppc_C1surv$getTrep()
res1surv <- compute_pB(T_obs1surv, T_rep1surv)

# Post-smolt
ppc_R2 <- nf_ppc3_scale(model_nimble, MCMC_matrix2, 'Mu_L2', 'Sd_L2', data$S2, const$I2, breaks_L2, midbins_L2)
ppc_C2 <- compileNimble(ppc_R2, project = compiled_model)

ppc_C2$run(MCMC_matrix2)
T_obs2 <- ppc_C2$getTobs()
T_rep2 <- ppc_C2$getTrep()
res2 <- compute_pB(T_obs2, T_rep2)

# Maturing post-smolt
ppc_R2m <- nf_ppc3_scale(model_nimble, MCMC_matrix2, 'Mu_L2m', 'Sd_L2m', data$S2m, const$I2m, breaks_L2, midbins_L2)
ppc_C2m <- compileNimble(ppc_R2m, project = compiled_model)

ppc_C2m$run(MCMC_matrix2)
T_obs2m <- ppc_C2m$getTobs()
T_rep2m <- ppc_C2m$getTrep()
res2m <- compute_pB(T_obs2m, T_rep2m)

# Non-maturing post-smolt
ppc_R2nm <- nf_ppc3_scale(model_nimble, MCMC_matrix2, 'Mu_L2nm', 'Sd_L2nm', data$S2nm, const$I2nm, breaks_L2, midbins_L2)
ppc_C2nm <- compileNimble(ppc_R2nm, project = compiled_model)

ppc_C2nm$run(MCMC_matrix2)
T_obs2nm <- ppc_C2nm$getTobs()
T_rep2nm <- ppc_C2nm$getTrep()
res2nm <- compute_pB(T_obs2nm, T_rep2nm)

# Dataframe
pB_df1 <- data.frame(
  year = 1996:2019,
  Smolt  = res1$pB_by_year,
  `Surviving smolt` = res1surv$pB_by_year,
  `Post-smolt` = res2$pB_by_year,
  `Maturing post-smolt` = res2m$pB_by_year,
  `Non-maturing post-smolt` = res2nm$pB_by_year,
  check.names = FALSE
)
print(pB_df1, row.names = FALSE)
write.csv2(pB_df1, file.path("results",project,"ppc_length3.csv"))
