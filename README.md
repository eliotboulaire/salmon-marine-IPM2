-   [Disentangling length-dependent and length-independent variations in
    survival and maturation in Atlantic
    salmon](#disentangling-length-dependent-and-length-independent-variations-in-survival-and-maturation-in-atlantic-salmon)
    -   [Introduction](#introduction)
        -   [What is this ?](#what-is-this)
    -   [Overview](#overview)
        -   [Abstract](#abstract)
        -   [Key words](#key-words)
        -   [Author(s)](#authors)
        -   [Corresponding author](#corresponding-author)
        -   [Acknowledgements](#acknowledgements)
    -   [Data availability &
        requirements](#data-availability-requirements)
        -   [Data availability](#data-availability)
        -   [Software and packages
            requirement](#software-and-packages-requirement)
    -   [File & folder structure](#file-folder-structure)
    -   [File descriptions](#file-descriptions)
        -   [R project](#r-project)
        -   [Manuscript](#manuscript)
        -   [Data](#data)
        -   [Functions](#functions)
        -   [Models](#models)
        -   [Run](#run)
        -   [Saves](#saves)
        -   [Graphs](#graphs)
        -   [Results](#results)
    -   [Execution instructions](#execution-instructions)
    -   [License](#license)

------------------------------------------------------------------------

# Disentangling length-dependent and length-independent variations in survival and maturation in Atlantic salmon

## Introduction

### What is this ?

Code and workflow accompanying the manuscript: “Disentangling
length-dependent and length-independent temporal variation in survival
and maturation in Atlantic salmon”.  
The project implements an Integrated Integral Projection Model (IPM²)
for the marine phase of the Scorff River Atlantic salmon (*Salmo salar*)
population, inferring length-dependent survival and maturation from
abundance, molecular-sexing and scale-length data.

------------------------------------------------------------------------

## Overview

### Abstract

Identifying the proximate mechanisms through which environmental changes
shape population dynamics is critical to analyse and forecast responses
under global change and support ecosystem-based management. Declines in
size-at-age have been widely reported across taxa, but understanding how
temporal variations in body size at key life stages propagate through
the life cycle and shape population dynamics remains challenging,
especially for species with partially observable life cycles. Atlantic
salmon (*Salmo salar*) have undergone concurrent declines in abundance,
length-at-age, and sea age at maturation, largely attributed to
deteriorating marine growth during the unobservable oceanic phase. The
marine phase of Atlantic salmon is shaped by a complex interplay between
survival and maturation schedules, yet the relative contributions of
body size to this interplay remain poorly understood. We developed a
Bayesian Integrated Integral Projection Model (IPM²) that integrates 24
years (1996-2019) of monitoring data from the Scorff River (France),
combining abundance estimates, molecular sexing, and length-at-stage
distributions reconstructed from archived fish scales. The model allowed
us to disentangle the relative contributions of length-dependent and
length-independent processes to temporal variation in survival and
maturation, separately for males and females. Length-at-stage explained
a substantial part of the pathways to variation in marine vital rates,
though considerable length-independent variation remained. Larger
individuals consistently experienced higher survival and higher
maturation than smaller conspecifics, with females maturing at lower
rates than males for the same length. However, temporal variation in the
smolt length distribution explained only 9% of the temporal variations
in marine survival. In contrast, temporal variation in the post-smolt
length distribution explained 45% of the temporal variation in marine
maturation for females and 21% for males. The IPM² framework is
transferable to other salmon and fish populations with partially
observable life cycles to understand how changing growth conditions
shape population dynamics under environmental change, and for building
climate-enhanced population dynamics and stock assessment models.

### Key words

body size; integral projection model; integrated population model;
life-history; *Salmo salar*; temporal variation; vital rates

### Author(s)

***Eliot BOULAIRE***\*\*  
UMR DECOD, Institut Agro, Rennes, France  
pôle MIAME,  
eliot.boulaire@institut-agro.fr

***Marie NEVOUX***\*  
UMR DECOD, INRAE, Rennes, France  
marie.nevoux@inrae.fr

***Etienne RIVOT***\*  
UMR DECOD, Institut Agro, Rennes, France  
etienne.rivot@institut-agro.fr

\*M. NEVOUX & E. RIVOT share senior co-authorship

### Corresponding author

\*\*Comments and requests should be addressed to Eliot BOULAIRE:
eliot.boulaire@institut-agro.fr or eliotboulaire@gmail.com. All material
is free of use, but I would appreciate being told, and this dataet and
the matching paper cited if appropriate.

### Acknowledgements

We are grateful to all those involved in the preparation and collection
of the data used in this study, including both current and former
members of the INRAE UE U3E, the UMR DECOD and the local angling
association. We especially thank Frédéric Marchand for overseeing scale
data collection as part of the COLISA program, Ludivine Lamireau for
conducting the scale readings, and Anne-Laure Besnard, Lisa Meslier, and
Thibaut Jousseaume for performing the genetic sexing. Also, a special
thanks to Cécile Tréhin for her earlier work on this topic, which laid
the groundwork for this publication. We thank Institut Agro, the
Brittany Region as well as OFB for funding this study.

------------------------------------------------------------------------

## Data availability & requirements

### Data availability

Raw data are publicly archived and must be placed under `data/rawdata/`
(see structure below) before running:

-   Trapping data (counts & biometry): GBIF archive
    (<https://doi.org/10.15468/yvcw8n>),

-   CMR abundance model (code & results): Zenodo archive
    (<https://doi.org/10.5281/zenodo.3275148>) & GitHub repository
    (<https://github.com/ORE-DiaPFC/Abundance>),

-   Scale samples & molecular sexing: COLISA collection
    (<https://doi.org/10.15454/D3ODJM>),

-   IPM² model (code & results; this repository): French data gouv
    archive (XXX) & GitHub repository
    (<https://github.com/eliotboulaire/salmon-marine-IPM2>).

### Software and packages requirement

R version 4.4.2 (or later). Rstudio 2024.12.0.467 (or later).  
MCMC via NIMBLE, which requires a working C++ toolchain (Rtools 4.4 or
later on Windows).

| Purpose | Packages |
|------------------------------------|------------------------------------|
| Data wrangling | `dplyr` (1.1.4); `tidyr` (1.3.1); `purrr` (1.0.4) |
| Running model (MCMC) | `nimble` (1.4.0); `parallel` (4.4.1); `coda` (0.19.4.1) |
| Saving models input & output | `qs` (0.27.3) |
| Model diagnostics & comparison | `ggmcmc` (1.5.1.1); `coda` (0.19.4.1) |
| Model comparison | `loo` (2.8.0) |
| Variance partitioning | `relaimpo` (2.2.7) |
| Figures | `ggplot2` (3.5.2); `ggridges` (3.5.2); `ggnewscale` (0.5.2) |

Each script auto-installs any missing packages via a small
`pkginstall()` helper at the top of the file.

------------------------------------------------------------------------

## File & folder structure

    ## .
    ## ├── data
    ## │   ├── data_create.R
    ## │   ├── data_save.R
    ## │   ├── rawdata
    ## │   │   ├── abundances
    ## │   │   │   ├── LogN1.csv
    ## │   │   │   ├── LogN3.csv
    ## │   │   │   └── LogN4.csv
    ## │   │   ├── data_Scorff_CMR.csv
    ## │   │   ├── data_Scorff_ORE.csv
    ## │   │   ├── scales
    ## │   │   │   ├── S1.csv
    ## │   │   │   ├── S1surv.csv
    ## │   │   │   ├── S2.csv
    ## │   │   │   ├── S2m.csv
    ## │   │   │   └── S2nm.csv
    ## │   │   └── sexes
    ## │   │       ├── Ns1.csv
    ## │   │       ├── Ns3.csv
    ## │   │       └── Ns4.csv
    ## │   └── realdata
    ## │       ├── const.qs
    ## │       ├── data.qs
    ## │       ├── M0
    ## │       │   ├── inits_1chain.qs
    ## │       │   ├── inits_chain1.qs
    ## │       │   ├── inits_chain2.qs
    ## │       │   ├── inits_chain3.qs
    ## │       │   ├── inits_chain4.qs
    ## │       │   ├── inits_chain5.qs
    ## │       │   ├── inits_chain6.qs
    ## │       │   ├── inits_nchains.qs
    ## │       │   ├── monitor1.qs
    ## │       │   └── monitor2.qs
    ## │       ├── M1
    ## │       │   ├── inits_1chain.qs
    ## │       │   ├── inits_chain1.qs
    ## │       │   ├── inits_chain2.qs
    ## │       │   ├── inits_chain3.qs
    ## │       │   ├── inits_chain4.qs
    ## │       │   ├── inits_chain5.qs
    ## │       │   ├── inits_chain6.qs
    ## │       │   ├── inits_nchains.qs
    ## │       │   ├── monitor1.qs
    ## │       │   └── monitor2.qs
    ## │       ├── M2
    ## │       │   ├── inits_1chain.qs
    ## │       │   ├── inits_chain1.qs
    ## │       │   ├── inits_chain2.qs
    ## │       │   ├── inits_chain3.qs
    ## │       │   ├── inits_chain4.qs
    ## │       │   ├── inits_chain5.qs
    ## │       │   ├── inits_chain6.qs
    ## │       │   ├── inits_nchains.qs
    ## │       │   ├── monitor1.qs
    ## │       │   └── monitor2.qs
    ## │       ├── M3
    ## │       │   ├── inits_1chain.qs
    ## │       │   ├── inits_chain1.qs
    ## │       │   ├── inits_chain2.qs
    ## │       │   ├── inits_chain3.qs
    ## │       │   ├── inits_chain4.qs
    ## │       │   ├── inits_chain5.qs
    ## │       │   ├── inits_chain6.qs
    ## │       │   ├── inits_nchains.qs
    ## │       │   ├── monitor1.qs
    ## │       │   └── monitor2.qs
    ## │       ├── M4
    ## │       │   ├── inits_1chain.qs
    ## │       │   ├── inits_chain1.qs
    ## │       │   ├── inits_chain2.qs
    ## │       │   ├── inits_chain3.qs
    ## │       │   ├── inits_chain4.qs
    ## │       │   ├── inits_chain5.qs
    ## │       │   ├── inits_chain6.qs
    ## │       │   ├── inits_nchains.qs
    ## │       │   ├── monitor1.qs
    ## │       │   └── monitor2.qs
    ## │       ├── M5
    ## │       │   ├── inits_1chain.qs
    ## │       │   ├── inits_chain1.qs
    ## │       │   ├── inits_chain2.qs
    ## │       │   ├── inits_chain3.qs
    ## │       │   ├── inits_chain4.qs
    ## │       │   ├── inits_chain5.qs
    ## │       │   ├── inits_chain6.qs
    ## │       │   ├── inits_nchains.qs
    ## │       │   ├── monitor1.qs
    ## │       │   └── monitor2.qs
    ## │       ├── M6
    ## │       │   ├── inits_1chain.qs
    ## │       │   ├── inits_chain1.qs
    ## │       │   ├── inits_chain2.qs
    ## │       │   ├── inits_chain3.qs
    ## │       │   ├── inits_chain4.qs
    ## │       │   ├── inits_chain5.qs
    ## │       │   ├── inits_chain6.qs
    ## │       │   ├── inits_nchains.qs
    ## │       │   ├── monitor1.qs
    ## │       │   └── monitor2.qs
    ## │       ├── M7
    ## │       │   ├── inits_1chain.qs
    ## │       │   ├── inits_chain1.qs
    ## │       │   ├── inits_chain2.qs
    ## │       │   ├── inits_chain3.qs
    ## │       │   ├── inits_chain4.qs
    ## │       │   ├── inits_chain5.qs
    ## │       │   ├── inits_chain6.qs
    ## │       │   ├── inits_nchains.qs
    ## │       │   ├── monitor1.qs
    ## │       │   └── monitor2.qs
    ## │       ├── M8
    ## │       │   ├── inits_1chain.qs
    ## │       │   ├── inits_chain1.qs
    ## │       │   ├── inits_chain2.qs
    ## │       │   ├── inits_chain3.qs
    ## │       │   ├── inits_chain4.qs
    ## │       │   ├── inits_chain5.qs
    ## │       │   ├── inits_chain6.qs
    ## │       │   ├── inits_M8_chain3.qs
    ## │       │   ├── inits_nchains.qs
    ## │       │   ├── monitor1.qs
    ## │       │   └── monitor2.qs
    ## │       └── M9
    ## │           ├── inits_1chain.qs
    ## │           ├── inits_chain1.qs
    ## │           ├── inits_chain2.qs
    ## │           ├── inits_chain3.qs
    ## │           ├── inits_chain4.qs
    ## │           ├── inits_chain5.qs
    ## │           ├── inits_chain6.qs
    ## │           ├── inits_nchains.qs
    ## │           ├── monitor1.qs
    ## │           └── monitor2.qs
    ## ├── functions
    ## │   ├── f_fillNA.R
    ## │   ├── f_geninit.R
    ## │   ├── f_round.R
    ## │   ├── nf_l.R
    ## │   ├── nf_pi.R
    ## │   └── nf_res.R
    ## ├── graphs
    ## │   ├── fig4_structure.R
    ## │   ├── fig5_size.R
    ## │   ├── fig6_effects.R
    ## │   ├── figA1_data.R
    ## │   ├── figA5_fit.R
    ## │   ├── figA7_cor.R
    ## │   ├── figA7_scatter.R
    ## │   ├── tab1_loo.R
    ## │   ├── tabA1_data.R
    ## │   ├── tabA4_conv.R
    ## │   ├── tabA4_sum.R
    ## │   ├── tabA5_ppc-abund.R
    ## │   ├── tabA5_ppc-alpha.R
    ## │   ├── tabA5_ppc-length.R
    ## │   └── tabA5_ppc-sex.R
    ## ├── LICENSE
    ## ├── manuscript
    ## ├── models
    ## │   ├── M0
    ## │   │   └── model_code.R
    ## │   ├── M1
    ## │   │   └── model_code.R
    ## │   ├── M2
    ## │   │   └── model_code.R
    ## │   ├── M3
    ## │   │   └── model_code.R
    ## │   ├── M4
    ## │   │   └── model_code.R
    ## │   ├── M5
    ## │   │   └── model_code.R
    ## │   ├── M6
    ## │   │   └── model_code.R
    ## │   ├── M7
    ## │   │   └── model_code.R
    ## │   ├── M8
    ## │   │   └── model_code.R
    ## │   └── M9
    ## │       └── model_code.R
    ## ├── README.md
    ## ├── README.Rmd
    ## ├── results
    ## │   ├── fig1.pdf
    ## │   ├── fig2.pdf
    ## │   ├── fig3.pdf
    ## │   ├── fig4.pdf
    ## │   ├── fig4_A.pdf
    ## │   ├── fig4_B.pdf
    ## │   ├── fig5.pdf
    ## │   ├── fig5_A.pdf
    ## │   ├── fig5_B.pdf
    ## │   ├── fig6.pdf
    ## │   ├── fig6_A.pdf
    ## │   ├── fig6_B.pdf
    ## │   ├── figA1.1.pdf
    ## │   ├── figA1.2.pdf
    ## │   ├── figA1.3.pdf
    ## │   ├── M0
    ## │   │   ├── tabA4.1.csv
    ## │   │   └── tabA4.2.csv
    ## │   ├── M1
    ## │   │   ├── tabA4.1.csv
    ## │   │   ├── tabA4.2.csv
    ## │   │   ├── tabA5.2.1.csv
    ## │   │   ├── tabA5.2.2.csv
    ## │   │   ├── tabA5.2.3.csv
    ## │   │   ├── tabA5.2.4.csv
    ## │   │   ├── tabA5.2.5.csv
    ## │   │   ├── tabA5.2.6A.csv
    ## │   │   └── tabA5.2.6B.csv
    ## │   ├── M2
    ## │   │   ├── tabA4.1.csv
    ## │   │   └── tabA4.2.csv
    ## │   ├── M3
    ## │   │   ├── tabA4.1.csv
    ## │   │   └── tabA4.2.csv
    ## │   ├── M4
    ## │   │   ├── tabA4.1.csv
    ## │   │   └── tabA4.2.csv
    ## │   ├── M5
    ## │   │   ├── tabA4.1.csv
    ## │   │   └── tabA4.2.csv
    ## │   ├── M6
    ## │   │   ├── tabA4.1.csv
    ## │   │   └── tabA4.2.csv
    ## │   ├── M7
    ## │   │   ├── tabA4.1.csv
    ## │   │   └── tabA4.2.csv
    ## │   ├── M8
    ## │   │   ├── tabA4.1.csv
    ## │   │   └── tabA4.2.csv
    ## │   ├── M9
    ## │   │   ├── tabA4.1.csv
    ## │   │   └── tabA4.2.csv
    ## │   ├── tab1.csv
    ## │   ├── tabA1.2.csv
    ## │   └── tabA1.3.csv
    ## ├── run
    ## │   ├── debug_model.R
    ## │   ├── model_save.R
    ## │   └── run_model.R
    ## ├── salmon-marine-IPM2.Rproj
    ## └── saves
    ##     ├── M0
    ##     │   └── MCMC.qs
    ##     ├── M1
    ##     │   └── MCMC.qs
    ##     ├── M2
    ##     │   └── MCMC.qs
    ##     ├── M3
    ##     │   └── MCMC.qs
    ##     ├── M4
    ##     │   └── MCMC.qs
    ##     ├── M5
    ##     │   └── MCMC.qs
    ##     ├── M6
    ##     │   └── MCMC.qs
    ##     ├── M7
    ##     │   └── MCMC.qs
    ##     ├── M8
    ##     │   └── MCMC.qs
    ##     └── M9
    ##         └── MCMC.qs

------------------------------------------------------------------------

## File descriptions

### R project

*MODESTA.Rproj* - R project that you need to simplify the running of all
the scripts.

### Manuscript

*Preprint_Manusrcipt.pdf* - PDF of the preprint, you can also check :

-   XXX (Link to the preprint),

-   XXX (Link to the final publication).

### Data

#### Rawdata

##### General

-   *data_Scorff_CMR.csv* - csv file of the original data file for
    abundance estimates.

| Columns | Description |
|-------------------|-----------------------------------------------------|
| stage | Life cycle stage, “Smolt” downstream migrating juvenile, “1SW” first year returning adult, “2SW” second year returning adult. |
| cohort | The year of downstream migrating juvenile (smolt) and returning adults (1SW & 2SW). dataets consist of years between 1996 and 2019. |
| mean | The stage & year mean abundance estimate from the lognormal CMR model. |
| sd | The stage & year standard-deviation abundance estimate from the lognormal CMR model. |

-   *data_Scorff_ORE.csv* - csv file of the original data file for scale
    length and scale molecular sexing.

| Columns | Description |
|-------------------|-----------------------------------------------------|
| index | Number for each distinct individual scale |
| stage | Life cycle stage: “Smolt” downstream migrating juvenile, “Adults” upstream migrating adults. |
| cohort | The year of downstream migrating juvenile (smolt) and returning adults (1SW & 2SW): dataets consist of years between 1996 and 2019. |
| sexe | Sex of the individual (not all individual where sexed): F for females, M for males and NA otherwise. |
| age | Age inside each stage: “1SW” or “2SW” for adults and “Smolt” for “Smolt” as all ages where pooled together. |
| migration | Scale length at the migration (when smolt migrates downstream to the sea): The mark was taken from adult scales retrospectively, for smolts it’s considered the edge of the scale. |
| end1sum | Scale length at the end of the first summer (when post-smolt are believed to “choose” between migrating = 1SW or waiting another year = 2SW : The mark was taken from adult scales retrospectively. |

-   *data_create.R* - R script reading both *data_Scorff* files to
    create dataets used in the model (*LogN*, *Ns*, and *S*)
-   *data_save.R* - R script using all dataets for the model (*LogN*,
    *Ns*, and *S*) to create model elements *data.qs* and *const.qs*

##### abundances/

-   *LogN1.csv* - csv of the estimated abundance of smolts in the Scorff
    river between cohort year 1996 and 2019 by the CMR model

-   *LogN2.csv* - csv of the estimated abundance of 1SW in the Scorff
    river between cohort year 1996 and 2019 by the CMR model

-   *LogN3.csv* - csv of the estimated abundance of 2SW in the Scorff
    river between cohort year 1996 and 2019 by the CMR model

##### sexes/

-   *Ns1.csv* - csv of the genotyped individual scales of smolts in the
    Scorff river between cohort year 1996 and 2019

-   *Ns2.csv* - csv of the genotyped individual scales of 1SW in the
    Scorff river between cohort year 1996 and 2019

-   *Ns3.csv* - csv of the genotyped individual scales of 2SW in the
    Scorff river between cohort year 1996 and 2019

##### scales/

-   *S1.csv* - csv of the individual scale length of the smolts (from
    smolt scales) in the Scorff river between cohort year 1996 and 2019

-   *S1surv.csv* - csv of the individual scale length of the surviving
    smolts (from adult scales) in the Scorff river between cohort year
    1996 and 2019

-   *S2.csv* - csv of the individual scale length of the post-smolts
    (from adult scales) in the Scorff river between cohort year 1996 and
    2019

-   *S2m.csv* - csv of the individual scale length of the maturing
    post-smolts (from 1SW adult scales) in the Scorff river between
    cohort year 1996 and 2019

-   *S2nm.csv* - csv of the individual scale length of the non-maturing
    post-smolts (from 2SW adult scales) in the Scorff river between
    cohort year 1996 and 2019

#### Realdata

Composed of models folders : M\[0-9\]/

##### General

-   *data.qs* - qs file of all the data used in the model

-   *const.qs* - qs file of all the constants used in the model

##### M\[0-9\]/

-   *monitor1.qs* - qs file of all monitored estimated parameters for
    results

-   *monitor2.qs* - qs file of all monitored data parameters for loo

-   *inits_chain\[1-6\].qs* - qs file of starting values for high level
    parameters to start all the MCMC chains

-   *inits_1chain.qs* - qs file of starting values for all parameters to
    start one MCMC chain

-   *inits_nchains.qs* - qs file of starting values for all parameters
    to start all the MCMC chains

### Functions

#### R functions

-   *f_fillNA.R* - R function that fill the scale length dataets with
    *NA* as all years don’t have the same sample size

-   *f_geninit.R* - R function that generates initial values for high
    level parameters to start MCMC chains.

-   *f_round.R* - R function that uses ceiling or floor functions
    depending to calculate the min and max values of scale length
    structure and ensure to include extreme values

#### Nimble functions

-   *nf_l.R* - nimbleFunction that calculate the mean length within each
    length class based on the mean and the sd of an assumed gaussian
    distribution of length at the population level.

-   *nf_pi.R* - nimbleFunction that calculate the proportion of
    individuals within each length class based on the mean and the sd of
    an assumed gaussian distribution of length at the population level.

-   *nf_res.R* - nimbleFunction that back calculate mean and sd of an
    assumed gaussian distribution of length at the population level
    based on the number of individuals within each length class and the
    mean length within each length class.

### Models

Composed of models folders : M\[0-9\]/

-   *model_code.R* - R code of the different models

### Run

-   *model_save.R* - R script used before running the model to set up
    all the variables before running the model. It construct the
    *monitor.qs* and *inits.qs* files

-   ***run_model.R*** **- R script to run the different models**

-   *debug_model.R* - R script to debug any problem when running a model

### Saves

Composed of models folders : M\[0-9\]/

-   *MCMC.qs* - qs file saved at the end of the MCMC procedure, it
    includes both monitored posteriors **(not included in the repository
    as too heavy)**

-   *run_info.txt* - txt about the information of the running model

### Graphs

#### Figure

-   *fig4_structure.R* - R script to recreate Figure 4: Annual variation
    in length distribution

-   *fig5_size.R* - R script to recreate Figure 5: Annual
    length-dependent vital rates

-   *fig6_effects.R* - R script to recreate Figure 6: Annual variation
    in population-level vital rates

-   *figA1_data.R* - R script to recreate figures of Annexe 1 (A1.1;
    A1.2 & A1.3): Annual data values

-   *figA5_fit.R* - R script to recreate figures of Annexe 5.1 (A5.1.1;
    A5.1.2 & A5.1.3): Annual estimated posterior distributions against
    data values

-   *figA7_cor.R* - R script to recreate figures of Annexe 7.1 (A7.1.1 &
    A7.1.2): Distribution of Pearson correlations between the
    sex-specific parameters of maturation

-   *figA7_scatter.R* - R script to recreate figures of Annexe 7.2
    (A7.2.1 & A7.2.2): Annual variation in length-dependent effect on
    population-level vital rates

#### Table

-   *tab1_loo.R* - R script to recreate Table 1: PSIS-LOO model
    comparison for the influence of sex, length-independent, and
    length-dependent effects on vital rates modelling

-   *tabA1_data.R* - R script to recreate tables of Annexe 1 (A1.2 &
    A1.3): Number of scale samples analysed to characterise sex-ratios
    and length distributions

-   *tabA4_conv.R* - R script to recreate table of Annexe 4 (A4.1):
    Convergence diagnosis of each model (Gelman-Rubin-Brooks diagnostic,
    Geweke diagnostic & ESS)

-   *tabA4_sum.R* - R script to recreate table of Annexe 4 (A4.2):
    Estimated parameters summaries of each model

-   *tabA5_ppc-abund.R* - R script to recreate table of Annexe 5.2
    (A5.2.1): Posterior Predictive Check (PPC) of abundance estimates

-   *tabA5_ppc-sex.R* - R script to recreate table of Annexe 5.2
    (A5.2.2): Posterior Predictive Check (PPC) of sex-ratio data

-   *tabA5_ppc-length.R* - R script to recreate table of Annexe 5.2
    (A5.2.3; A5.2.4; A5.2.5): Posterior Predictive Check (PPC) of
    length-distribution data. 3) Proportions; 4) Mean length; 5) 95% IQR

-   *tabA5_ppc-alpha.R* - R script to recreate table of Annexe 5.2
    (A5.2.6): Posterior Predictive Check (PPC) of vital rates alpha
    parameters estimates

### Results

Composed of models folders : M\[0-9\]/

#### General

-   *tab1.csv* - Table 1: PSIS-LOO model comparison for the influence of
    sex, length-independent, and length-dependent effects on vital rates
    modelling

#### M\[0-9\]/

-   *figA5.1.1.pdf* - Figure A5.1.1: Annual estimated abundance
    posterior distributions against data values

-   *figA5.1.2.pdf* - Figure A5.1.2: Annual estimated sex-ratio
    posterior distributions against data values

-   *figA5.1.3.pdf* - Figure A5.1.3: Annual estimated
    length-distribution posterior distributions against data values

<!-- -->

-   *tabA4.1.csv* - Table A4.1: Convergence diagnosis of each model
    (Gelman-Rubin-Brooks diagnostic, Geweke diagnostic & ESS)

-   tabA4.2.csv *-* Table A4.2: Estimated parameters summaries of each
    model

#### M1/

-   *tabA5.2.1.pdf* - Table A5.2.1: Posterior Predictive Check (PPC) of
    abundance estimates

-   *tabA5.2.2.pdf* - Table A5.2.2: Posterior Predictive Check (PPC) of
    sex-ratio data

-   *tabA5.2.3.pdf* - Table A5.2.3: Posterior Predictive Check (PPC) of
    length-distribution data using proportions test

-   *tabA5.2.4.pdf* - Table A5.2.4: Posterior Predictive Check (PPC) of
    length-distribution data using mean length test

-   *tabA5.2.5.pdf* - Table A5.2.5: Posterior Predictive Check (PPC) of
    length-distribution data using 95% IQR test

-   *tabA5.2.6A.pdf* - Table A5.2.6A: Posterior Predictive Check (PPC)
    of survival alpha parameter estimates

-   *tabA5.2.6B.pdf* - Table A5.2.6B: Posterior Predictive Check (PPC)
    of maturation alpha parameters estimates

------------------------------------------------------------------------

## Execution instructions

1.  You first need to run *data_create.R*, *data_save.R* (at this step
    you can also check data using *figA1_data.R* & *tabA1_data.R*).

2.  You will also need to run the *model_save.R*.

3.  Then you can use the *run_model.R* script (if any issues occur,
    please verify previous script and use the *debug_model.R*).

4.  Finally, when *MCMC.qs* and *run_info.txt* output are created, you
    can use the graphs script. We recommend the following order:

    1.  *tabA4_conv.R* & *tabA4_sum.R*

    2.  *tab1_loo.R*

    3.  *figA5_fit.R* & *tabA5_ppc-abund.R*; *tabA5_ppc-sex.R*;
        *tabA5_ppc-length.R*; *tabA5_ppc-alpha.R*

    4.  *fig4_structure.R*; *fig5_size.R* & *fig6_effects.R*

    5.  *figA7_cor.R* & *figA7_scatter.R*

------------------------------------------------------------------------

## License

The content of this project itself is licensed under the [GNU GPL v3.0
license](https://www.gnu.org/licenses/gpl-3.0.en.html), except the
dataets which are under the [Creative Commons Attribution 4.0
license](https://creativecommons.org/licenses/by/4.0/).
