### ============================================================================
### Disentangling length-dependent and length-independent variations
### in survival and maturation in Atlantic salmon
### ----------------------------------------------------------------------------
### File:    nf_pi.R
### Purpose: NIMBLE function giving the proportion of individuals in each length
###          class when a Normal length distribution is discretised into
###          `nbclass` classes (Pi_{st,l,t} in Appendix S2).
###          Each class proportion is the Normal probability mass over the class
###          interval, obtained from differences of the cumulative distribution.
### Author:  ©BOULAIRE Eliot, NEVOUX Marie & RIVOT Etienne
### Version: 01/09/2026
### ============================================================================

## -----------------------------------------------------------------------------
## Class proportions Pi_{st,l,t}  (Appendix S2.2)
##    mu, sd    : mean and sd of the Normal length distribution
##    nbclass   : number of length classes (60 in the fitted model)
##    Lmin/Lmax : bounds of the length axis
##    Returns a vector of length `nbclass` giving each class's probability mass.
## -----------------------------------------------------------------------------
nf_pi <- nimbleFunction(
  run = function(mu = double(0), sd = double(0), nbclass = double(0), Lmin = double(0), Lmax = double(0)) {

    omega <- nimC(mu, sd)
    
    breaks <- seq(Lmin, Lmax, length.out = (nbclass + 1))
    
    prop_intervals <- pnorm(breaks, mean = omega[1], sd = omega[2])
    prop <- prop_intervals[2:(nbclass + 1)] - prop_intervals[1:(nbclass)]
    
    returnType(double(1))
    return(prop)
  }
)