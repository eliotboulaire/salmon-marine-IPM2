### ============================================================================
### Disentangling length-dependent and length-independent variations
### in survival and maturation in Atlantic salmon
### ----------------------------------------------------------------------------
### File:    f_round.R
### Purpose: R function to round a value down (floor) or up (ceiling)
###          to a chosen number of decimals.
###          Used to set the bounds of the length axis (Min_L, Max_L) so that 
###          all observed scale lengths fall inside the length range.
### Author:  ©BOULAIRE Eliot, NEVOUX Marie & RIVOT Etienne
### Version: 01/09/2026
### ============================================================================

## -----------------------------------------------------------------------------
## Directional rounding
##   x      : value to round
##   digit  : number of decimal places to keep
##   method : 1 = round down (floor), 2 = round up (ceiling)
## -----------------------------------------------------------------------------
f_round <- function(x, digit = 1, method = 1) {
  ratio <- 10^digit
  
  if (method == 1) {
    return(floor(x * ratio) / ratio)  # Arrondit vers le bas
  } else if (method == 2) {
    return(ceiling(x * ratio) / ratio)  # Arrondit vers le haut
  } else {
    stop("Invalid method. Use 1 for floor or 2 for ceiling.")  # Message d'erreur si méthode invalide
  }
}