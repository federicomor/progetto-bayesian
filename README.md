# Project A2

Goal: Clustering weekly data of one year of PM10 (plus covariates - see AGRIMONIA project)   
Models & R packages:  
- drpm (Page, Quintana, Dahl (2022) "Dependent Modeling of Temporal Sequences of Random Partitions", JCGS, R-package on Github)    
- ppmSuite (various models implemented, also PPMs) on https://cran.r-project.org/web/packages/ppmSuite/index.html   
Data: old https://zenodo.org/records/7563265 new https://zenodo.org/records/7956006   


# Models order
**Important**: for models fitting let's use `df_weekly_scaled_centered` loaded from `include.R`, which has (or will have, i will update it soon):

- numerical covariates scaled
- spatial coordinates also scaled (fitting appears to be better)
- PM10 just centered, not scaled

And for models evaluation see file `Metrics.Rmd` (i will create it soon).

**Models recap**:

- DRPM (Federico): with time and space
- Gaussian PPMx (Federica): just covariates
- SPPM (Oswaldo): just space
- Curve ppmx (Aby): time and covariates

| Model             | Time     | Space    | Covariates |
|-------------------|----------|----------|------------|
| DRPM (Federico)          |  ✓      | ✓        |  X    |
| Gaussian PPMx (Federica) |  X        |      X  |    ✓ |
| SPPM (Oswaldo)           |  X        |       ✓   | X    |
| Curve PPMx (Aby)         |   ✓     | X      |    ✓   |


The models without time means that time is not inside the model. Then we can "force" time to be included by looping or doing some other tricks, but time will never be inside the model.

The models without space could consider to use space as a covariate? like the locations longitude and latitude as numerical covariates.
