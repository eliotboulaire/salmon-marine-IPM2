### ============================================================================
### Disentangling length-dependent and length-independent variations
### in survival and maturation in Atlantic salmon
### ----------------------------------------------------------------------------
### File:    nf_pi.R
### Purpose: NIMBLE function returning the mean and sd of a population-level
###          scale-length distribution, given the number of individuals in each
###          length class and the mean length of each class.
###          Used to summarize the reconstructed length distributions of
###          survivors and of maturing / non-maturing post-smolts as an
###          abundance-weighted Normal (mu, sd), which then feeds the likelihood.
### Author:  ©BOULAIRE Eliot, NEVOUX Marie & RIVOT Etienne
### Version: 01/09/2026
### ============================================================================

## -----------------------------------------------------------------------------
## Abundance-weighted length moments
##    N : [length class x sex] matrix of individual numbers per class
##    L : vector of class mean lengths (from nf_l)
##    Returns c(mu, sd): mean and sd of the weighted length distribution.
## -----------------------------------------------------------------------------

nf_res <- nimbleFunction(
  run = function(N = double(2), L = double(1)) {

    nRows <- dim(N)[1]
    N_sum <- numeric(nRows)
    
    for (i in 1:nRows) {
      N_sum[i] <- sum(N[i, ])
    }
    total_N <- sum(N_sum)
    
    mu <- sum(L * N_sum) / total_N
    sd <- sqrt(sum(N_sum * (L - mu)^2) / total_N)
    omega <- nimC(mu, sd)
               
    returnType(double(1))
    return(omega)
  }
)