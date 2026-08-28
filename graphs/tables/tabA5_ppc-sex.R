### ============================================================================
### Creating table of Posterior Predictive Check on sexes data
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
nf_ppc_sex <- nimbleFunction(
  
  setup = function(model, samples, obsNodes, probNodes, sizeFixed) {
    
    theta <- colnames(samples)
    theta <- model$topologicallySortNodes(theta)
    theta <- model$expandNodeNames(theta, returnScalarComponents = TRUE)
    deps  <- model$getDependencies(theta, self = FALSE)
    
    # block names (one block per c, holding both categories together)
    obsNodesExpanded  <- model$expandNodeNames(obsNodes,  returnScalarComponents = FALSE)
    probNodesExpanded <- model$expandNodeNames(probNodes, returnScalarComponents = FALSE)
    n <- length(obsNodesExpanded)
    
    sizeFixed <- as.numeric(sizeFixed)   # Nt, fixed trial sizes, length n
    
    # fixed observed FEMALE counts, read once
    obsDataFem <- numeric(n)
    for (j in 1:n) {
      vals <- values(model, obsNodesExpanded[j])
      obsDataFem[j] <- vals[1]
    }
    
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
      
      for (j in 1:n) {
        
        probFitVals <- values(model, probNodesExpanded[c])
        gammaFem    <- probFitVals[1]                 # fitted female probability
        varFem      <- gammaFem * (1 - gammaFem)
        
        # --- T_obs: observed female proportion vs fitted Gamma_f ---
        propObs <- obsDataFem[j] / sizeFixed[j]
        T_obs_store[i, j] <<- (propObs - gammaFem) / varFem
      }
      
      # --- draw ONE multinomial replicate per class c, this iteration ---
      model$simulate(obsNodes, includeData = TRUE)
      
      for (j in 1:n) {
        
        probFitVals <- values(model, probNodesExpanded[c])
        gammaFem    <- probFitVals[1]
        varFem      <- gammaFem * (1 - gammaFem)
        
        repVals    <- values(model, obsNodesExpanded[j])
        propRep    <- repVals[1] / sizeFixed[j]        # Nt_obs used, per your spec
        T_rep_store[i, j] <<- (propRep - gammaFem) / varFem
      }
    }
  },
  
  # accessor methods -- call these after run() to retrieve each matrix separately
  methods = list(
    getTobs = function() {
      returnType(double(2));
      return(T_obs_store) },
    getTrep = function() {
      returnType(double(2));
      return(T_rep_store)
    }
  )
)

## =============================================================================
## 4. Build, compile, and run for each stage: smolts (1), 1SW (3), 2SW (4)
## =============================================================================

## ---- Ns1 ----
ppc_R1 <- nf_ppc_sex(model_nimble, MCMC_matrix2, 'Ns1', 'Gamma1', const$Nt1)
ppc_C1 <- compileNimble(ppc_R1, project = compiled_model)

ppc_C1$run(MCMC_matrix2)
T_obs1 <- ppc_C1$getTobs()
T_rep1 <- ppc_C1$getTrep()

## ---- Ns3 ----
ppc_R3 <- nf_ppc_sex(model_nimble, MCMC_matrix2, 'Ns3', 'Gamma3', const$Nt3)
ppc_C3 <- compileNimble(ppc_R3, project = compiled_model)

ppc_C3$run(MCMC_matrix2)
T_obs3 <- ppc_C3$getTobs()
T_rep3 <- ppc_C3$getTrep()

## ---- Ns4 ----
ppc_R4 <- nf_ppc_sex(model_nimble, MCMC_matrix2, 'Ns4', 'Gamma4', const$Nt4)
ppc_C4 <- compileNimble(ppc_R4, project = compiled_model)

ppc_C4$run(MCMC_matrix2)
T_obs4 <- ppc_C4$getTobs()
T_rep4 <- ppc_C4$getTrep()

## =============================================================================
## 5. Bayesian p-values
##
## p_B = P( T(y_rep, theta) >= T(y_obs, theta) ), estimated across posterior draws.
## Values outside [0.05, 0.95] indicate the model fails to replicate the
## observed discrepancy pattern (Gelman et al., 2021).
## =============================================================================
compute_pB <- function(T_obs, T_rep) {
  list(
    pB_by_year = colMeans(T_rep >= T_obs),               # one p_B per year
    pB_pooled  = mean(rowSums(T_rep) >= rowSums(T_obs))   # one p_B, all years pooled
  )
}

res1 <- compute_pB(T_obs1, T_rep1)   # smolts
res3 <- compute_pB(T_obs3, T_rep3)   # 1SW
res4 <- compute_pB(T_obs4, T_rep4)   # 2SW

pB_df1 <- data.frame(
  year = 1996:2019,     # adjust if c doesn't index years 1:1
  Smolt  = res1$pB_by_year,
  `1SW`  = res2$pB_by_year,
  `2SW`  = res3$pB_by_year,
  check.names = FALSE
)
print(pB_df1, row.names = FALSE)
write.csv2(pB_df1, file.path("results",project,"ppc_sex.csv"))


pB_df2 <- pB_df1 %>%
  pivot_longer(cols = -year, names_to = "stage", values_to = "pB") %>%
  mutate(
    stage = factor(stage, levels = c("Smolts", "1SW", "2SW")),
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
