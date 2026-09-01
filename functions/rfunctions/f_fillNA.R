### ============================================================================
### Disentangling length-dependent and length-independent variations
### in survival and maturation in Atlantic salmon
### ----------------------------------------------------------------------------
### File:    f_fillNA.R
### Purpose: Compact each column of a data frame with its non-missing values,
###          then equalize columns to the same length with trailing NAs.
###          Used in data_create.R to turn per-cohort scale length samples
###          (with unequal numbers of scales per year) into a rectangular
###          [scale x cohort] matrix.
### Author:  ©BOULAIRE Eliot, NEVOUX Marie & RIVOT Etienne
### Version: 01/09/2026
### ============================================================================

## -----------------------------------------------------------------------------
## Compact each column values and equalize columns to the same length
##   df : data frame with columns of unequal sample size
##   Returns a data frame where each column holds its non-NA values at the top,
##   followed by trailing NAs so all columns share the length of the longest.
## -----------------------------------------------------------------------------
f_fillNA <- function(df) {
  df_cleaned <- apply(df, 2, na.omit)
  
  max_length <- max(sapply(df_cleaned, length))
  
  for (i in seq_along(df_cleaned)) {
    diff_length <- max_length - length(df_cleaned[[i]])
    if (diff_length > 0) {
      df_cleaned[[i]] <- c(df_cleaned[[i]], rep(NA, diff_length))
    }
  }
  
  df_cleaned <- as.data.frame(df_cleaned)
  
  return(df_cleaned)
}