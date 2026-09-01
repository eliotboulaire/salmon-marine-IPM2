### ============================================================================
### Disentangling length-dependent and length-independent variations
### in survival and maturation in Atlantic salmon
### ----------------------------------------------------------------------------
### File:    f_fillNA.R
### Purpose: Generate a set of initial values for the model's stochastic nodes,
###          drawn from vague distributions around plausible values so each MCMC
###          chain can start from a different, valid point. Called once per chain
###          in data_save.R with distinct seeds.
###          `n` is a length-16 vector giving how many values to draw for each
###          node group (typically C for per-cohort nodes, S for per-sex nodes,
###          1 for scalars), in the fixed order below.
### Author:  ©BOULAIRE Eliot, NEVOUX Marie & RIVOT Etienne
### Version: 01/09/2026
### ============================================================================

## -----------------------------------------------------------------------------
## Draw one set of initial values
##   n : integer vector of length 16 draw (typically C for per-cohort nodes,
##       S for per-sex nodes, 1 for scalars, in the fixed order below)
##   Returns a named list of inits.
## -----------------------------------------------------------------------------
f_geninit <- function(n) {
  out <- list(
    Mu_N1 = runif(n = n[1], min = 1200, max = 12000),
    Gamma1_F = rbeta(n = n[2], shape1 = 1, shape2 = 1),
    
    Mu_LogMu_L1 = rnorm(n = n[3], mean = 0, sd = 0.01),
    Sd_LogMu_L1 = rnorm(n = n[4], mean = 0.05, sd = 0.005),
    Mu_LogSd_L1 = rnorm(n = n[5], mean = -1.5, sd = 0.05),
    Sd_LogSd_L1 = rnorm(n = n[6], mean = 0.05, sd = 0.005),
    
    Beta1 = rnorm(n = n[7], mean = 0, sd = 5),
    Mu_Alpha1 = rnorm(n = n[8], mean = 2, sd = 0.1),
    Sd_Alpha1 = rnorm(n = n[9], mean = 0.2, sd = 0.02),
    
    Mu_LogMu_L2 = rnorm(n = n[10], mean = 0.6, sd = 0.02),
    Sd_LogMu_L2 = rnorm(n = n[11], mean = 0.05, sd = 0.005),
    Mu_LogSd_L2 = rnorm(n = n[12], mean = -1.5, sd = 0.05),
    Sd_LogSd_L2 = rnorm(n = n[13], mean = 0.05, sd = 0.005),
    
    Beta2 = rnorm(n = n[14], mean = 0, sd = 5),
    Mu_Alpha2 = rnorm(n = n[15], mean = 2, sd = 0.1),
    Sd_Alpha2 = rnorm(n = n[16], mean = 0.2, sd = 0.02)
  )
  return(out)
}