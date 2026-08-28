### ============================================================================
### Creating table of Posterior Predictive Check on scale length data
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
MCMC_matrix2 <- MCMC_matrix[, theta_order]   # reordered/expanded to match NIMBLE's internal graph

## =============================================================================
## 3. NIMBLE function: compute chi-squared discrepancies T_obs and T_rep
##    for one data source (e.g. LogN1, LogN3, or LogN4), for every posterior draw
## =============================================================================

## =============================================================================
## PPC test 1 -- mean of scale lengths vs fitted mu (chi2)
## =============================================================================
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
    
    obsMat <- as.matrix(obsMat)   # ensure matrix, not data.frame
    
    obsMean <- numeric(nC)
    for (y in 1:nC) {
      n <- nIndivPerCohort[y]
      x <- obsMat[1:n, y]        # trim NA padding using the true sample size
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


## =============================================================================
## PPC test 2 -- 90% interquantile range of scale lengths vs fitted sd (chi2)
## =============================================================================
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

## =============================================================================
## PPC test 3 -- binned length-class proportions vs fitted Normal(mu, sd) (chi2)
## =============================================================================
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

## =============================================================================
## 4. Build, compile, and run for each stage: smolts (1), 1SW (3), 2SW (4)
## =============================================================================

## =============================================================================
## PPC test 1 -- 
## =============================================================================

## ---- s1 ----
ppc_R1 <- nf_ppc1_scale(model_nimble, MCMC_matrix2, 'Mu_L1', 'Sd_L1', data$S1, const$I1)
ppc_C1 <- compileNimble(ppc_R1, project = compiled_model)

ppc_C1$run(MCMC_matrix2)
T_obs1 <- ppc_C1$getTobs()
T_rep1 <- ppc_C1$getTrep()

## ---- s1surv ----
ppc_R1surv <- nf_ppc1_scale(model_nimble, MCMC_matrix2, 'Mu_L1surv', 'Sd_L1surv', data$S1surv, const$I1surv)
ppc_C1surv <- compileNimble(ppc_R1surv, project = compiled_model)

ppc_C1surv$run(MCMC_matrix2)
T_obs1surv <- ppc_C1surv$getTobs()
T_rep1surv <- ppc_C1surv$getTrep()

## ---- s2 ----
ppc_R2 <- nf_ppc1_scale(model_nimble, MCMC_matrix2, 'Mu_L2', 'Sd_L2', data$S2, const$I2)
ppc_C2 <- compileNimble(ppc_R2, project = compiled_model)

ppc_C2$run(MCMC_matrix2)
T_obs2 <- ppc_C2$getTobs()
T_rep2 <- ppc_C2$getTrep()

## ---- s2m ----
ppc_R2m <- nf_ppc1_scale(model_nimble, MCMC_matrix2, 'Mu_L2m', 'Sd_L2m', data$S2m, const$I2m)
ppc_C2m <- compileNimble(ppc_R2m, project = compiled_model)

ppc_C2m$run(MCMC_matrix2)
T_obs2m <- ppc_C2m$getTobs()
T_rep2m <- ppc_C2m$getTrep()

## ---- s2nm ----
ppc_R2nm <- nf_ppc1_scale(model_nimble, MCMC_matrix2, 'Mu_L2nm', 'Sd_L2nm', data$S2nm, const$I2nm)
ppc_C2nm <- compileNimble(ppc_R2nm, project = compiled_model)

ppc_C2nm$run(MCMC_matrix2)
T_obs2nm <- ppc_C2nm$getTobs()
T_rep2nm <- ppc_C2nm$getTrep()


compute_pB <- function(T_obs, T_rep) {
  list(
    pB_by_year = colMeans(T_rep >= T_obs),               # one p_B per year
    pB_pooled  = mean(rowSums(T_rep) >= rowSums(T_obs))   # one p_B, all years pooled
  )
}

res1 <- compute_pB(T_obs1, T_rep1)
res1surv <- compute_pB(T_obs1surv, T_rep1surv)
res2 <- compute_pB(T_obs2, T_rep2)
res2m <- compute_pB(T_obs2m, T_rep2m)
res2nm <- compute_pB(T_obs2nm, T_rep2nm)

pB_df1 <- data.frame(
  year = 1996:2019,     # adjust if c doesn't index years 1:1
  Smolt  = res1$pB_by_year,
  `Surviving smolt` = res1surv$pB_by_year,
  `Post-smolt` = res2$pB_by_year,
  `Maturing post-smolt` = res2m$pB_by_year,
  `Non-maturing post-smolt` = res2nm$pB_by_year,
  check.names = FALSE
)
print(pB_df1, row.names = FALSE)
write.csv2(pB_df1, file.path("results",project,"ppc_length1.csv"))

pB_df2 <- pB_df1 %>%
  pivot_longer(cols = -year, names_to = "stage", values_to = "pB") %>%
  mutate(
    stage = factor(stage, levels = c("Smolt", "Surviving smolt", "Post-smolt", "Maturing post-smolt", "Non-maturing post-smolt")),
    dev   = pB - 0.5
  )
ggplot(pB_df2) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.05, ymax = 0.95,
           alpha = 0.2, fill = "grey", colour = "black") +
  geom_hline(yintercept = 0.5, linewidth = 1.1) +
  geom_point(aes(x = year, y = pB), colour = "black", size = 5) +
  geom_point(aes(x = year, y = pB, colour = dev), size = 3.8) +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-0.5, 0.5),
                         name = "Deviation\nfrom 0.5") +
  facet_wrap(~ stage, ncol = 1) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Year of smolt migration", y = "Bayesian p-value") +
  theme_minimal()

## =============================================================================
## PPC test 2 -- 
## =============================================================================

## ---- s1 ----
ppc_R1 <- nf_ppc2_scale(model_nimble, MCMC_matrix2, 'Mu_L1', 'Sd_L1', data$S1, const$I1, 0.25)
ppc_C1 <- compileNimble(ppc_R1, project = compiled_model)

ppc_C1$run(MCMC_matrix2)
T_obs1 <- ppc_C1$getTobs()
T_rep1 <- ppc_C1$getTrep()

## ---- s1surv ----
ppc_R1surv <- nf_ppc2_scale(model_nimble, MCMC_matrix2, 'Mu_L1surv', 'Sd_L1surv', data$S1surv, const$I1surv, 0.25)
ppc_C1surv <- compileNimble(ppc_R1surv, project = compiled_model)

ppc_C1surv$run(MCMC_matrix2)
T_obs1surv <- ppc_C1surv$getTobs()
T_rep1surv <- ppc_C1surv$getTrep()

## ---- s2 ----
ppc_R2 <- nf_ppc2_scale(model_nimble, MCMC_matrix2, 'Mu_L2', 'Sd_L2', data$S2, const$I2, 0.25)
ppc_C2 <- compileNimble(ppc_R2, project = compiled_model)

ppc_C2$run(MCMC_matrix2)
T_obs2 <- ppc_C2$getTobs()
T_rep2 <- ppc_C2$getTrep()

## ---- s2m ----
ppc_R2m <- nf_ppc2_scale(model_nimble, MCMC_matrix2, 'Mu_L2m', 'Sd_L2m', data$S2m, const$I2m, 0.25)
ppc_C2m <- compileNimble(ppc_R2m, project = compiled_model)

ppc_C2m$run(MCMC_matrix2)
T_obs2m <- ppc_C2m$getTobs()
T_rep2m <- ppc_C2m$getTrep()

## ---- s2nm ----
ppc_R2nm <- nf_ppc2_scale(model_nimble, MCMC_matrix2, 'Mu_L2nm', 'Sd_L2nm', data$S2nm, const$I2nm, 0.25)
ppc_C2nm <- compileNimble(ppc_R2nm, project = compiled_model)

ppc_C2nm$run(MCMC_matrix2)
T_obs2nm <- ppc_C2nm$getTobs()
T_rep2nm <- ppc_C2nm$getTrep()

compute_pB <- function(T_obs, T_rep) {
  list(
    pB_by_year = colMeans(T_rep >= T_obs),               # one p_B per year
    pB_pooled  = mean(rowSums(T_rep) >= rowSums(T_obs))   # one p_B, all years pooled
  )
}

res1 <- compute_pB(T_obs1, T_rep1)
res1surv <- compute_pB(T_obs1surv, T_rep1surv)
res2 <- compute_pB(T_obs2, T_rep2)
res2m <- compute_pB(T_obs2m, T_rep2m)
res2nm <- compute_pB(T_obs2nm, T_rep2nm)

pB_df1 <- data.frame(
  year = 1996:2019,     # adjust if c doesn't index years 1:1
  Smolt  = res1$pB_by_year,
  `Surviving smolt` = res1surv$pB_by_year,
  `Post-smolt` = res2$pB_by_year,
  `Maturing post-smolt` = res2m$pB_by_year,
  `Non-maturing post-smolt` = res2nm$pB_by_year,
  check.names = FALSE
)
print(pB_df1, row.names = FALSE)
write.csv2(pB_df1, file.path("results",project,"ppc_length2_0.25.csv"))

pB_df2 <- pB_df1 %>%
  pivot_longer(cols = -year, names_to = "stage", values_to = "pB") %>%
  mutate(
    stage = factor(stage, levels = c("Smolt", "Surviving smolt", "Post-smolt", "Maturing post-smolt", "Non-maturing post-smolt")),
    dev   = pB - 0.5
  )
ggplot(pB_df2) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.05, ymax = 0.95,
           alpha = 0.2, fill = "grey", colour = "black") +
  geom_hline(yintercept = 0.5, linewidth = 1.1) +
  geom_point(aes(x = year, y = pB), colour = "black", size = 5) +
  geom_point(aes(x = year, y = pB, colour = dev), size = 3.8) +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-0.5, 0.5),
                         name = "Deviation\nfrom 0.5") +
  facet_wrap(~ stage, ncol = 1) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Year of smolt migration", y = "Bayesian p-value") +
  theme_minimal()

## =============================================================================
## PPC test 3 -- 
## =============================================================================
breaks_L1 <- seq(const$Min_L1, const$Max_L1, length.out = (const$L + 1))
midbins_L1 <- (breaks_L1[1:const$L] - breaks_L1[2:(const$L+1)])/2
breaks_L2 <- seq(const$Min_L2, const$Max_L2, length.out = (const$L + 1))
midbins_L2 <- (breaks_L2[1:const$L] - breaks_L2[2:(const$L+1)])/2

## ---- s1 ----
ppc_R1 <- nf_ppc3_scale(model_nimble, MCMC_matrix2, 'Mu_L1', 'Sd_L1', data$S1, const$I1, breaks_L1, midbins_L1)
ppc_C1 <- compileNimble(ppc_R1, project = compiled_model)

ppc_C1$run(MCMC_matrix2)
T_obs1 <- ppc_C1$getTobs()
T_rep1 <- ppc_C1$getTrep()

## ---- s1surv ----
ppc_R1surv <- nf_ppc3_scale(model_nimble, MCMC_matrix2, 'Mu_L1surv', 'Sd_L1surv', data$S1surv, const$I1surv, breaks_L1, midbins_L1)
ppc_C1surv <- compileNimble(ppc_R1surv, project = compiled_model)

ppc_C1surv$run(MCMC_matrix2)
T_obs1surv <- ppc_C1surv$getTobs()
T_rep1surv <- ppc_C1surv$getTrep()

## ---- s2 ----
ppc_R2 <- nf_ppc3_scale(model_nimble, MCMC_matrix2, 'Mu_L2', 'Sd_L2', data$S2, const$I2, breaks_L2, midbins_L2)
ppc_C2 <- compileNimble(ppc_R2, project = compiled_model)

ppc_C2$run(MCMC_matrix2)
T_obs2 <- ppc_C2$getTobs()
T_rep2 <- ppc_C2$getTrep()

## ---- s2m ----
ppc_R2m <- nf_ppc3_scale(model_nimble, MCMC_matrix2, 'Mu_L2m', 'Sd_L2m', data$S2m, const$I2m, breaks_L2, midbins_L2)
ppc_C2m <- compileNimble(ppc_R2m, project = compiled_model)

ppc_C2m$run(MCMC_matrix2)
T_obs2m <- ppc_C2m$getTobs()
T_rep2m <- ppc_C2m$getTrep()

## ---- s2nm ----
ppc_R2nm <- nf_ppc3_scale(model_nimble, MCMC_matrix2, 'Mu_L2nm', 'Sd_L2nm', data$S2nm, const$I2nm, breaks_L2, midbins_L2)
ppc_C2nm <- compileNimble(ppc_R2nm, project = compiled_model)

ppc_C2nm$run(MCMC_matrix2)
T_obs2nm <- ppc_C2nm$getTobs()
T_rep2nm <- ppc_C2nm$getTrep()

compute_pB <- function(T_obs, T_rep) {
  list(
    pB_by_year = colMeans(T_rep >= T_obs),               # one p_B per year
    pB_pooled  = mean(rowSums(T_rep) >= rowSums(T_obs))   # one p_B, all years pooled
  )
}

res1 <- compute_pB(T_obs1, T_rep1)
res1surv <- compute_pB(T_obs1surv, T_rep1surv)
res2 <- compute_pB(T_obs2, T_rep2)
res2m <- compute_pB(T_obs2m, T_rep2m)
res2nm <- compute_pB(T_obs2nm, T_rep2nm)

pB_df1 <- data.frame(
  year = 1996:2019,     # adjust if c doesn't index years 1:1
  Smolt  = res1$pB_by_year,
  `Surviving smolt` = res1surv$pB_by_year,
  `Post-smolt` = res2$pB_by_year,
  `Maturing post-smolt` = res2m$pB_by_year,
  `Non-maturing post-smolt` = res2nm$pB_by_year,
  check.names = FALSE
)
print(pB_df1, row.names = FALSE)
write.csv2(pB_df1, file.path("results",project,"ppc_length3.csv"))

pB_df2 <- pB_df1 %>%
  pivot_longer(cols = -year, names_to = "stage", values_to = "pB") %>%
  mutate(
    stage = factor(stage, levels = c("Smolt", "Surviving smolt", "Post-smolt", "Maturing post-smolt", "Non-maturing post-smolt")),
    dev   = pB - 0.5
  )
ggplot(pB_df2) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.05, ymax = 0.95,
           alpha = 0.2, fill = "grey", colour = "black") +
  geom_hline(yintercept = 0.5, linewidth = 1.1) +
  geom_point(aes(x = year, y = pB), colour = "black", size = 5) +
  geom_point(aes(x = year, y = pB, colour = dev), size = 3.8) +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-0.5, 0.5),
                         name = "Deviation\nfrom 0.5") +
  facet_wrap(~ stage, ncol = 1) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Year of smolt migration", y = "Bayesian p-value") +
  theme_minimal()