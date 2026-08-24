-   [Introduction](#introduction)
    -   [What is this ?](#what-is-this)
-   [Overview](#overview)
    -   [Abstract](#abstract)
    -   [Key words](#key-words)
    -   [Author(s)](#authors)
    -   [Corresponding author](#corresponding-author)
    -   [Acknowledgements](#acknowledgements)
-   [Data availability & requirements](#data-availability-requirements)
    -   [Data availability](#data-availability)
    -   [Software and packages
        requirement](#software-and-packages-requirement)
-   [File & folder structure](#file-folder-structure)
-   [File descriptions](#file-descriptions)
    -   [R project](#r-project)
    -   [Manuscript](#manuscript)
    -   [Data](#data)
        -   [Rawdata](#rawdata)
        -   [Realdata](#realdata)
        -   [R functions](#r-functions)
        -   [Nimble functions](#nimble-functions)
    -   [Models](#models)
    -   [Run](#run)
    -   [Saves](#saves)
    -   [Graphs](#graphs)
        -   [Figure](#figure)
        -   [Table](#table)
    -   [Results](#results)
        -   [General](#general-2)
        -   [M\[0-9\]/](#m0-9-1)
        -   [M1/](#m1)
-   [Execution instructions](#execution-instructions)
-   [License](#license)

------------------------------------------------------------------------

# Introduction

## What is this ?

Code and workflow accompanying the manuscript. The project implements an
Integrated Integral Projection Model (IPM²) for the marine phase of the
Scorff River Atlantic salmon (*Salmo salar*) population, inferring
length-dependent survival and maturation from abundance,
molecular-sexing and scale-length data.

------------------------------------------------------------------------

# Overview

## Abstract

1.  Identifying the proximate mechanisms through which environmental
    changes shape population dynamics is critical to forecast responses
    under global change and support ecosystem-based management.

2.  Global change is driving widespread declines in size-at-age across
    taxa. Therefore, understanding how variations in the size
    distribution at key life stages propagate through the life cycle and
    shape population dynamics is critical. However, this remains
    particularly challenging for species with partially observable life
    cycles.

3.  Atlantic salmon (*Salmo salar*) have undergone concurrent declines
    in abundance, length-at-age, and sea age at maturation, largely
    attributed to deteriorating marine growth during an unobservable
    oceanic phase. In this paper, we investigate how temporal variation
    in the distribution of body length at seaward migration (smolts) and
    after the first summer at sea (post-smolts) drive temporal variation
    in marine survival and maturation.

4.  We developed a Bayesian Integrated Integral Projection Model (IPM²)
    that integrates 24 years (1996-2019) of monitoring data from the
    Scorff River (France), combining abundance estimates, molecular
    sexing, and length-at-stage distributions reconstructed from
    archived fish scales. The framework inferred latent length-dependent
    survival and maturation and quantified the relative contribution of
    the variation in length distribution to temporal variation in those
    vital rates.

5.  Results showed that length-at-stage explained a substantial part of
    the pathways to variations in vital rates, yet considerable
    length-independent variation remained. Larger individuals
    consistently experienced higher survival and higher maturation than
    smaller conspecifics, with females maturing at lower rates than
    males for the same length. Temporal variations in the smolt length
    distribution explained only 9% of the temporal variations in marine
    survival. In contrast, the variation in the post-smolt length
    distribution explained 45% of the temporal variation in maturation
    rates of females and 21% in males.

6.  The IPM² framework is transferable to other salmon and fish
    populations and offers a mechanistic basis for climate-enhanced
    population dynamics and stock assessment in length-structured
    species.

## Key words

body length; global change; integral projection model; integrated
population model; life-history; *Salmo salar*; temporal variations;
vital rates

## Author(s)

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

**\***M. NEVOUX & E. RIVOT share senior co-authorship

## Corresponding author

\*\*Comments and requests should be addressed to Eliot BOULAIRE:
eliot.boulaire@institut-agro.fr or eliotboulaire@gmail.com. All material
is free of use, but I would appreciate being told, and this dataset and
the matching paper cited if appropriate.

## Acknowledgements

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

# Data availability & requirements

## Data availability

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
    archive (XXX) & GitHub repository (XXX).

## Software and packages requirement

R version 4.4.2 (or later).  
MCMC via NIMBLE, which requires a working C++ toolchain (Rtools on
Windows; Xcode command-line tools on macOS; `build-essential` on Linux).

| Purpose                        | Packages                        |
|--------------------------------|---------------------------------|
| Data wrangling                 | `dplyr` (), `tidyr` ()          |
| Running model (MCMC)           | `nimble` (1.3.0), `parallel` () |
| Saving models output           | `qs` ()                         |
| Model diagnostics & comparison | `ggmcmc` (1.5.1.1), `coda` ()   |
| Mode comparison                | `loo` (2.8.0)                   |
| Variance partitioning          | `relaimpo` (2.2-7)              |
| Figures                        | `ggplot2`                       |

Each script auto-installs any missing packages via a small
`pkginstall()` helper at the top of the file.

------------------------------------------------------------------------

# File & folder structure

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
    ## ├── README.Rmd
    ## ├── results
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

# File descriptions

## R project

*MODESTA.Rproj* - R project that you need to simplify the running of all
the scripts.

## Manuscript

*Preprint_Manusrcipt.pdf* - PDF of the preprint, you can also check :

-   <https://www.biorxiv.org/content/10.64898/2025.12.30.696970v1> (Link
    to the preprint),

-   XXX (Link to the final publication).

## Data

### Rawdata

#### General

-   *data_Scorff_CMR.csv* - csv file of the original data file for
    abundance estimates.

| Columns | Description |
|---------------|---------------------------------------------------------|
| stage | Life cycle stage, “Smolt” downstream migrating juvenile, “1SW” first year returning adult, “2SW” second year returning adult. |
| cohort | The year of downstream migrating juvenile (smolt) and returning adults (1SW & 2SW). Datasets consist of years between 1995 and 2019. |
| mean | The stage & year mean abundance estimate from the lognormal CMR model. |
| sd | The stage & year standard-deviation abundance estimate from the lognormal CMR model. |

-   *data_Scorff_ORE.csv* - csv file of the original data file for scale
    length and scale molecular sexing.

| Columns | Description |
|--------------|----------------------------------------------------------|
| index | Number for each distinct individual scale |
| stage | Life cycle stage: “Smolt” downstream migrating juvenile, “Adults” upstream migrating adults. |
| cohort | The year of downstream migrating juvenile (smolt) and returning adults (1SW & 2SW): Datasets consist of years between 1995 and 2019. |
| sexe | Sex of the individual (not all individual where sexed): F for females, M for males and NA otherwise. |
| age | Age inside each stage: “1SW” or “2SW” for adults and “Smolt” for “Smolt” as all ages where pooled together. |
| migration | Scale length at the migration (when smolt migrates downstream to the sea): The mark was taken from adult scales retrospectively, for smolts it’s considered the edge of the scale. |
| end1sum | Scale length at the end of the first summer (when post-smolt are believed to “choose” between migrating = 1SW or waiting another year = 2SW : The mark was taken from adult scales retrospectively. |

-   *data_create.R* - R script reading both *data_Scorff* files to
    create datasets used in the model (*LogN*, *Ns*, and *S*)

<!-- -->

-   *data_save.R* - R script using all datasets for the model (*LogN*,
    *Ns*, and *S*) to create model elements *data.qs* and *const.qs*

#### abundances/

-   *LogN1.csv* - csv of the estimated abundance of smolts in the Scorff
    river between cohort year 1995 and 2019 by the CMR model

-   *LogN2.csv* - csv of the estimated abundance of 1SW in the Scorff
    river between cohort year 1995 and 2019 by the CMR model

-   *LogN3.csv* - csv of the estimated abundance of 2SW in the Scorff
    river between cohort year 1995 and 2019 by the CMR model

#### sexes/

-   *Ns1.csv* - csv of the genotyped individual scales of smolts in the
    Scorff river between cohort year 1995 and 2019

-   *Ns2.csv* - csv of the genotyped individual scales of 1SW in the
    Scorff river between cohort year 1995 and 2019

-   *Ns3.csv* - csv of the genotyped individual scales of 2SW in the
    Scorff river between cohort year 1995 and 2019

#### scales/

-   *S1.csv* - csv of the individual scale length of the smolts (from
    smolt scales) in the Scorff river between cohort year 1995 and 2019

-   *S1surv.csv* - csv of the individual scale length of the surviving
    smolts (from adult scales) in the Scorff river between cohort year
    1995 and 2019

-   *S2.csv* - csv of the individual scale length of the post-smolts
    (from adult scales) in the Scorff river between cohort year 1995 and
    2019

-   *S2m.csv* - csv of the individual scale length of the maturing
    post-smolts (from 1SW adult scales) in the Scorff river between
    cohort year 1995 and 2019

-   *S2nm.csv* - csv of the individual scale length of the non-maturing
    post-smolts (from 2SW adult scales) in the Scorff river between
    cohort year 1995 and 2019

### Realdata

Composed of models folders : M\[0-9\]/

#### General

-   *data.qs* - qs file of all the datas used in the model

-   *const.qs* - qs file of all the constants used in the model

#### M\[0-9\]/

-   *monitor1.qs* - qs file of all monitored estimated parameters for
    results

-   *monitor2.qs* - qs file of all monitored data parameters for loo

-   *inits_chain\[1-6\].qs* - qs file of starting values for high level
    parameters to start all the MCMC chains

-   *inits_1chain.qs* - qs file of starting values for all parameters to
    start one MCMC chain

-   *inits_nchains.qs* - qs file of starting values for all parameters
    to start all the MCMC chains

    ## Functions

### R functions

-   *f_fillNA.R* - R function that fill the scale length datasets with
    *NA* as all years don’t have the same sample size

-   *f_geninit.R* - R function that generates initial values for high
    level parameters to start MCMC chains.

-   *f_round.R* - R function that uses ceiling or floor functions
    depending to calculate the min and max values of scale length
    structure and ensure to include extreme values

### Nimble functions

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

## Models

Composed of models folders : M\[0-9\]/

-   *model_code.R* - R script of the different models

## Run

-   *model_save.R* - R script used before running the model to set up
    all the variables before running the model. It construct the
    *monitor.qs* and *inits.qs* files

-   ***run_model.R*** **- R script to run the different models**

-   *debug_model.R* - R script to debug any problem when running a model

## Saves

Composed of models folders : M\[0-9\]/

-   *MCMC.qs* - qs file saved at the end of the MCMC procedure, it
    includes both monitored posteriors **(not included in the repository
    as too heavy)**

-   *run_info.txt* - txt about the information of the running model

## Graphs

### Figure

-   *fig4_structure.R* -

-   *fig5_size.R* -

-   *fig6_effects.R* -

-   *figA1_data.R* -

-   *figA5_fit.R* -

-   *figA7_cor.R* -

-   *figA7_scatter.R* -

### Table

-   *tab1_loo.R* -

-   *tabA1_data.R* -

-   *tabA4_conv.R* -

-   *tabA4_sum.R* -

-   *tabA5_ppc-abund.R* -

-   *tabA5_ppc-sex.R* -

-   *tabA5_ppc-length.R* -

-   *tabA5_ppc-alpha.R* -

## Results

Composed of models folders : M\[0-9\]/

### General

-   *tab1.csv* -

### M\[0-9\]/

-   *figA5.1.1.pdf* -

-   *figA5.1.2.pdf* -

-   *figA5.1.3.pdf* -

<!-- -->

-   *tabA4.1.csv* -

-   tabA4.2.csv *-*

### M1/

-   *tabA5.2.1.pdf* -

-   *tabA5.2.2.pdf* -

-   *tabA5.2.3.pdf* -

-   *tabA5.2.4.pdf* -

-   *tabA5.2.5.pdf* -

-   *tabA5.2.6A.pdf* -

-   *tabA5.2.6B.pdf* -

------------------------------------------------------------------------

# Execution instructions

You need to run

------------------------------------------------------------------------

# License

The content of this project itself is licensed under the [GNU GPL v3.0
license](https://www.gnu.org/licenses/gpl-3.0.en.html), except the
datasets which are under the [Creative Commons Attribution 4.0
license](https://creativecommons.org/licenses/by/4.0/).
