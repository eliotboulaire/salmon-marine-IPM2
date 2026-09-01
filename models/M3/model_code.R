### ============================================================================
### Disentangling length-dependent and length-independent variations
### in survival and maturation in Atlantic salmon
### ----------------------------------------------------------------------------
### File:    model_code.R (M3)
### Purpose: NIMBLE model code for the model (M3) of the IPM2.
###          Length-structured marine phase of the salmon life cycle.
###          Composed of 60-length-class, 24 years, 2-sex process linking three
###          data sources (abundance, sex ratio, scale length) to latent states.
###          - Survival (Stage 1->2) is a length-dependent logistic function,
###          - Maturation (Stage 2 -> 2m/2nm) is a length-dependent logistic function,
###          - Post-maturation survival (stage 2m/2nm- > 3/4) is fixed.
### Author:  ©BOULAIRE Eliot, NEVOUX Marie & RIVOT Etienne
### Version: 01/09/2026
### ============================================================================

model_code <- nimbleCode({
  ## ===========================================================================
  ##                  GLOSSARY OF PRINCIPAL VARIABLE NAMES
  ## ===========================================================================
  ## Indices   c = cohort/year of smolt migration (1:C, C=24);
  ##           l = length class (1:L, L=60);
  ##           s = sex (1 = female, 2 = male).
  ## Abund.    N1 = smolts;
  ##           N1surv Surviving smolts;
  ##           N2 post-smolts;
  ##           N2m/N2nm maturing/non-maturing post-smolts
  ##           N3/N4 returning 1SW/2SW adults.
  ## Length    Mu_L*/Sd_L* Normal mean/sd of scale length;
  ##           L* class means (nf_l);
  ##           Pi* class proportions (nf_pi).
  ## Rates     Theta1 survival (per length class and cohort);
  ##           Theta2 maturation (per length class, cohort and sex);
  ##           Theta3/Theta4 fixed post-maturation survival (constant).
  ##           Alpha* length-independent random intercept (per cohort & sex);
  ##           Beta* length effect (slope), constant.
  ## Obs.      LogN* (1, 3 & 4) CMR log-abundance;
  ##           Ns* (1, 3 & 4) genotyped (F,M) counts;
  ##           S* (1, 1surv, 2, 2m & 2nm) measured individual scale lengths.
  ## ===========================================================================
  
  ## ===========================================================================
  ##                          STAGE 1 -> STAGE 1surv
  ##   Smolts at sea entry -> Surviving smolt (length-dependent survival)
  ## ===========================================================================
  
  ## ---------------------------------------------------------------------------
  ##                                 PRIOR
  ## ---------------------------------------------------------------------------
  
  # Stage 1
  # ------------------------------------
  # Abundance parameters (Mu_N1)
  for (c in 1:C) {
    Mu_LogN1[c] <- log(Mu_N1[c]) - ((1/2) * Sd_LogN1[c])
    Mu_N1[c] ~ dlnorm(meanlog = 8, sdlog = 2)
  }
  
  # Length structure parameters (Mu_L1 & Sd_L1)
  # + Derived parameters (L1 & Pi1)
  Mu_LogMu_L1 ~ dnorm(mean = 1, sd = 1)
  Sd_LogMu_L1 ~ T(dt(0, 1/(1^2), df = 1), 0, )
  Mu_Mu_L1 <- exp(Mu_LogMu_L1+((Sd_LogMu_L1^2)/2))
  
  Mu_LogSd_L1 ~ dnorm(mean = 1, sd = 1)
  Sd_LogSd_L1 ~ T(dt(0, 1/(1^2), df = 1), 0, )
  Mu_Sd_L1 <- exp(Mu_LogSd_L1+((Sd_LogSd_L1^2)/2))
  
  for (c in 1:C) {
    Mu_L1[c] ~ dlnorm(meanlog = Mu_LogMu_L1, sdlog = Sd_LogMu_L1)
    Sd_L1[c] ~ dlnorm(meanlog = Mu_LogSd_L1, sdlog = Sd_LogSd_L1)
    
    L1[1:L, c] <- nf_l(mu = Mu_L1[c], sd = Sd_L1[c], nbclass = L, Lmin = Min_L1, Lmax = Max_L1)
    Pi1[1:L, c] <- nf_pi(mu = Mu_L1[c], sd = Sd_L1[c], nbclass = L, Lmin = Min_L1, Lmax = Max_L1)
  }
  
  # Sex ratio parameters (Gamma1)
  for (c in 1:C) {
    Gamma1_F[c] ~ dbeta(shape1 = 2, shape2 = 2)
    Gamma1[c, 1:2] <- nimC(Gamma1_F[c], (1 - Gamma1_F[c]))
  }
  
  # Length-dependent survival
  # ------------------------------------
  # Logistic parameters (Alpha1 & Beta1)
  Beta1 ~ dnorm(mean = 0, sd = 10)
  Mu_Alpha1 ~ dnorm(mean = 0, sd = 10)
  Sd_Alpha1 ~ T(dt(0, 1/(1^2), 1), 0, )
  
  # Survival rate (Theta 1)
  for (c in 1:C) {
    Alpha1[c] ~ dnorm(mean = Mu_Alpha1, sd = Sd_Alpha1)
    
    Logit_Theta1[1:L, c] <- Alpha1[c] + Beta1 * (L1[1:L, c] - Mean_L1)
    Theta1[1:L, c] <- (1 / (1 + exp(-Logit_Theta1[1:L, c])))
  }
  
  # Population dynamics
  # ------------------------------------
  # Application of the survival rate
  for (c in 1:C) {
    for (s in 1:S) {
      N1[1:L, c, s] <- Mu_N1[c] * Gamma1[c, s] * Pi1[1:L, c]
      N1surv[1:L, c, s] <- N1[1:L, c, s] * Theta1[1:L, c]
      
      Mu_N2[c, s] <- sum(N1surv[1:L, c, s])
    }
  }
  
  # Stage 1surv
  # ------------------------------------
  # Length structure parameters (Mu_L1surv & Sd_L1surv)
  for (c in 1:C) {
    Mu_L1surv[c] <- nf_res(N = N1surv[1:L, c, 1:S], L = L1[1:L, c])[1]
    Sd_L1surv[c] <- nf_res(N = N1surv[1:L, c, 1:S], L = L1[1:L, c])[2]
  }
  Mu_Mu_L1surv <- mean(Mu_L1surv[1:C])
  Mu_Sd_L1surv <- mean(Sd_L1surv[1:C])
  
  ## ---------------------------------------------------------------------------
  ##                                LIKELIHOOD
  ## ---------------------------------------------------------------------------
  
  # Length structure likelihood
  # -----------------------------
  for (c in 1:C) {
    for (i in 1:I1[c]) {
      S1[i, c] ~ dnorm(mean = Mu_L1[c], sd = Sd_L1[c])
    }
    
    for (i in 1:I1surv[c]) {
      S1surv[i, c] ~ dnorm(mean = Mu_L1surv[c], sd = Sd_L1surv[c])
    }
  }
  
  # Abundance likelihood
  # -----------------------------
  for (c in 1:C) {
    LogN1[c] ~ dlnorm(meanlog = Mu_LogN1[c], sdlog = Sd_LogN1[c])
  }
  
  # Sex ratio likelihood
  # -----------------------------
  for (c in 1:C) {
    Ns1[c, 1:2] ~ dmulti(prob = Gamma1[c, 1:2], size = Nt1[c])
  }
  
  
  ## ===========================================================================
  ##                      STAGE 2 -> STAGE 2M & STAGE 2NM
  ##             Post-smolt at the end of the first summer
  ##                -> Maturing post-smolt 
  ##                -> Non-maturing post-smolt
  ##             => Length-dependent maturation
  ## ===========================================================================
  
  ## ---------------------------------------------------------------------------
  ##                                 PRIOR
  ## ---------------------------------------------------------------------------
  
  # Stage 2
  # -----------------------------
  # Length structure parameters (Mu_L2 & Sd_L2)
  # + Derived parameters (L2 & Pi2)
  Mu_LogMu_L2 ~ dnorm(mean = 1, sd = 1)
  Sd_LogMu_L2 ~ T(dt(0, 1/(1^2), df = 1), 0, )
  Mu_Mu_L2 <- exp(Mu_LogMu_L2+((Sd_LogMu_L2^2)/2))
  
  Mu_LogSd_L2 ~ dnorm(mean = 1, sd = 1)
  Sd_LogSd_L2 ~ T(dt(0, 1/(1^2), df = 1), 0, )
  Mu_Sd_L2 <- exp(Mu_LogSd_L2+((Sd_LogSd_L2^2)/2))
  
  for (c in 1:C) {
    Mu_L2[c] ~ dlnorm(meanlog = Mu_LogMu_L2, sdlog = Sd_LogMu_L2)
    Sd_L2[c] ~ dlnorm(meanlog = Mu_LogSd_L2, sdlog = Sd_LogSd_L2)
    
    L2[1:L, c] <- nf_l(mu = Mu_L2[c], sd = Sd_L2[c], nbclass = L, Lmin = Min_L2, Lmax = Max_L2)
    Pi2[1:L, c] <- nf_pi(mu = Mu_L2[c], sd = Sd_L2[c], nbclass = L, Lmin = Min_L2, Lmax = Max_L2)
  }
  
  # Length-dependent maturation
  # ------------------------------------
  # Logistic parameters (Alpha2 & Beta2)
  Beta2 ~ dnorm(mean = 0, sd = 10)
  Mu_Alpha2 ~ dnorm(mean = 0, sd = 10)
  Sd_Alpha2 ~ T(dt(0, 1/(1^2), 1), 0, )
  
  # Maturation rate (Theta 2)
  for (c in 1:C) {
    Alpha2[c] ~ dnorm(mean = Mu_Alpha2, sd = Sd_Alpha2)
    
    Logit_Theta2[1:L, c] <- Alpha2[c] + Beta2 * (L1[1:L, c] - Mean_L2)
    Theta2[1:L, c] <- (1 / (1 + exp(-Logit_Theta2[1:L, c])))
  }
  
  # Population dynamics
  # ------------------------------------
  # Application of the maturation rate
  for (c in 1:C) {
    for (s in 1:S) {
      N2[1:L, c, s] <- Mu_N2[c, s] * Pi2[1:L, c]
      
      N2m[1:L, c, s] <- N2[1:L, c, s] * Theta2[1:L, c]
      N2nm[1:L, c, s] <- N2[1:L, c, s] * (1-Theta2[1:L, c])
      
      Mu_N2m[c, s] <- sum(N2m[1:L, c, s])
      Mu_N2nm[c, s] <- sum(N2nm[1:L, c, s]) 
    }
  }
  
  # Stage 2m
  # ------------------------------------
  # Length structure parameters (Mu_L2m & Sd_L2m)
  for (c in 1:C) {
    Mu_L2m[c] <- nf_res(N = N2m[1:L, c, 1:S], L = L2[1:L, c])[1]
    Sd_L2m[c] <- nf_res(N = N2m[1:L, c, 1:S], L = L2[1:L, c])[2]
  }
  Mu_Mu_L2m <- mean(Mu_L2m[1:C])
  Mu_Sd_L2m <- mean(Sd_L2m[1:C])
  
  # Stage 2nm
  # ------------------------------------
  # Length structure parameters (Mu_L2nm & Sd_L2nm)
  for (c in 1:C) {
    Mu_L2nm[c] <- nf_res(N = N2nm[1:L, c, 1:S], L = L2[1:L, c])[1]
    Sd_L2nm[c] <- nf_res(N = N2nm[1:L, c, 1:S], L = L2[1:L, c])[2]
  }
  Mu_Mu_L2nm <- mean(Mu_L2nm[1:C])
  Mu_Sd_L2nm <- mean(Sd_L2nm[1:C])
  
  ## ---------------------------------------------------------------------------
  ##                                LIKELIHOOD
  ## ---------------------------------------------------------------------------
  
  # Length structure likelihood
  # ------------------------------------
  for (c in 1:C) {
    for (i in 1:I2[c]) {
      S2[i,c] ~ dnorm(mean = Mu_L2[c], sd = Sd_L2[c])
    }
    
    for (i in 1:I2m[c]) {
      S2m[i,c] ~ dnorm(mean = Mu_L2m[c], sd = Sd_L2m[c])
    }
    
    for (i in 1:I2nm[c]) {
      S2nm[i,c] ~ dnorm(mean = Mu_L2nm[c], sd = Sd_L2nm[c])
    }
  }
  
  
  ## ===========================================================================
  ##                  STAGE 2M & STAGE 2NM -> STAGE 3 & STAGE 4
  ##             Maturing post-smolt -> 1SW adults
  ##             Non-maturing post-smolt -> 2SW adults
  ##             => Additionnal fixed survival
  ## ===========================================================================
  
  ## ---------------------------------------------------------------------------
  ##                                 PRIOR
  ## ---------------------------------------------------------------------------
  
  # Population dynamics
  # ------------------------------------
  # Application of the fixed survival rate
  for (s in 1:S) {
    for (c in 1:C) {
      N3[c, s] <- Mu_N2m[c, s] * Theta3
      N4[c, s] <- Mu_N2nm[c, s] * Theta4
    }
  }
  
  # Stage 3
  # ------------------------------------
  for (c in 1:C) {
  # Abundance parameters
    Mu_N3[c] <- N3[c, 1] + N3[c, 2]
    Mu_LogN3[c] <- log(Mu_N3[c]) - ((1/2) * Sd_LogN3[c])
    
  # Sex ratio parameters
    Gamma3[c, 1] <- N3[c, 1]/Mu_N3[c]
    Gamma3[c, 2] <- N3[c, 2]/Mu_N3[c]
  }
  
  # Stage 4
  # ------------------------------------
  for (c in 1:C) {
  # Abundance parameters
    Mu_N4[c] <- N4[c, 1] + N4[c, 2]
    Mu_LogN4[c] <- log(Mu_N4[c]) - ((1/2) * Sd_LogN4[c])
    
  # Sex ratio parameters
    Gamma4[c, 1] <- N4[c, 1]/Mu_N4[c]
    Gamma4[c, 2] <- N4[c, 2]/Mu_N4[c]
  }
  
  ## ---------------------------------------------------------------------------
  ##                                LIKELIHOOD
  ## ---------------------------------------------------------------------------
  
  # Abundance likelihood
  # ------------------------------------
  for (c in 1:C) {
    LogN3[c] ~ dlnorm(meanlog = Mu_LogN3[c], sdlog = Sd_LogN3[c])
    
    LogN4[c] ~ dlnorm(meanlog = Mu_LogN4[c], sdlog = Sd_LogN4[c])
  }
  
  # Sex ratio likelihood
  # ------------------------------------
  for (c in 1:C) {
    Ns3[c, 1:2] ~ dmulti(prob = Gamma3[c, 1:2], size = Nt3[c])
    
    Ns4[c, 1:2] ~ dmulti(prob = Gamma4[c, 1:2], size = Nt4[c])
  }
  
  
  ## ===========================================================================
  ##                            END OF THE MODEL
  ## ===========================================================================
})