### ============================================================================
### Disentangling length-dependent and length-independent variations
### in survival and maturation in Atlantic salmon
### ----------------------------------------------------------------------------
### File:    nf_l.R
### Purpose: NIMBLE functions giving the mean length within each length class
###          when a Normal length distribution is discretised into `nbclass`
###          classes (L_{st,l,t} in Appendix S2).
###          Each class mean is the density-weighted average of length
###          over the class interval.
### Author:  ©BOULAIRE Eliot, NEVOUX Marie & RIVOT Etienne
### Version: 01/09/2026
### ============================================================================

## -----------------------------------------------------------------------------
## 1. Numerator integrand x*f(x) (Lambda, Appendix S2.2).
## -----------------------------------------------------------------------------
nf1 <- nimbleFunction(
  run = function(x = double(1), omega = double(1)) {
    
    numerat <- x * dnorm(x, omega[1], omega[2])
    
    returnType(double(1))
    return(numerat)
  }
)

## -----------------------------------------------------------------------------
## 2. Denominator integrand f(x)  (Omega, Appendix S2.2).
## -----------------------------------------------------------------------------
nf2 <- nimbleFunction(
  run = function(x = double(1), omega = double(1)) {
    
    denomin <- dnorm(x, omega[1], omega[2])
    
    returnType(double(1))
    return(denomin)
  }
)

## -----------------------------------------------------------------------------
## 3. Class means L_{st,l,t}  (L, Appendix S2.2)
##    mu, sd    : mean and sd of the Normal length distribution
##    nbclass   : number of length classes (60 in the fitted model)
##    Lmin/Lmax : bounds of the length axis
##    Returns a vector of length `nbclass` giving each class's mean lengths.
## -----------------------------------------------------------------------------
nf_l <- nimbleFunction(
  run = function(mu = double(0), sd = double(0), nbclass = double(0), Lmin = double(0), Lmax = double(0)) {
    
    omega <- nimC(mu, sd)
    
    breaks <- seq(Lmin, Lmax, length.out = (nbclass + 1))
    
    int_num <- numeric(length = nbclass)
    int_denom <- numeric(length = nbclass)
    midbins <- numeric(length = nbclass)
    
    for (l in 1:nbclass) {
      int_num[l] <- nimIntegrate(nf1, lower = breaks[l], upper = breaks[l + 1], param = omega)[1]
      
      int_denom[l] <- nimIntegrate(nf2, lower = breaks[l], upper = breaks[l + 1], param = omega)[1]
    
      if (int_denom[l] == 0 | is.na(int_denom[l]) | is.nan(int_denom[l])) {
        midbins[l] <- (breaks[l] + breaks[l + 1])/2
      } else {
        midbins[l] <- int_num[l] / int_denom[l]
      }
    }
    
    returnType(double(1))
    return(midbins)
  }
)