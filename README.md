# Fitting individual-based models of spatial population dynamics to long-term monitoring data


This is supplementary material to the publication: 

*"Fitting individual-based models of spatial population dynamics to long-term monitoring data"*

by Anne-Kathleen Malchow, Guillermo Fandos, Urs G. Kormann, Martin U. Grüebler, Marc Kéry, Florian Hartig, Damaris Zurell.  

It contains the R scripts and data needed to reproduce all results. The code uses a developmental version of the R package *RangeShiftR*, which can be found in the RangeShiftR repository under the development tag [v.1.1-beta.1](https://github.com/RangeShifter/RangeShiftR-package/releases/tag/v.1.1-beta.1).

This development version implements extended functionality of the interface to replace all file in- and outputs. Inputs such as habitat maps can now be given as matrices in R and the output abundance maps are (optionally) returned as raster objects. 
(Note: the output type that is stored in these return rasters is currently hard-coded and is here set to the number of adult individuals in each cell. Adults are all stages with strictly positive fecundity.) 
Avoiding the file connection saves runtime and facilitates thread safety as needed for parallel code execution.

**Funding**

AM and DZ were supported by Deutsche Forschungsgemeinschaft (DFG) under grant agreement No. ZU 361/1-1.


**Data provision**

The study uses data from the Swiss breeding bird atlas 1993-96 ([Schmid et al., 1998](#1)) and the Swiss breeding bird index ([Knaus et al., 2022](#2)) provided by the Swiss ornithological institute, Sempach.
The repository contains processed versions of this data.

<a id="1"></a>
Schmid, H., Luder, R., Naef-Daenzer, B., Graf, R., & Zbinden, N. (1998). Schweizer Brutvogelatlas 1993–1996. Schweizerische Vogelwarte.

<a id="2"></a>
Knaus, P., Strebel, N., & Sattler, T. (2022). The State of Birds in Switzerland 2022. Swiss Ornithological Institute. http://www.vogelwarte.ch/state


---

## Folder *data*

All input data required to run the RangeShiftR model and the calibration: 
- **habitatmaps**: ASCII rasters of yearly habitat maps (for the standard scenario as well as for the two sensitivity scenarios, labelled "plus"and "mnus")
- **init**: the initial distribution model
- **abund**: spatially aggregated counts from the Swiss breeding bird survey
- **spatial_blocking**: the spatial blocks of the gridded landscape over which to aggregate the counts / number of individuals for comparison in the calibration
- **hillshade**: Switzerland hillshading for plotting

## Folder *scripts*

#### *scripts/calibration*

The scripts to run the model calibration. The procedural calibration routine is contained in *RangeShiftR_calibration.R*, which loads functions defined in the other scripts.

#### *scripts/analysis*

The scripts to process and analyse the MCMC chains produced in the calibration. Reads and combines the independent MCMCs for each fold; plots diagnostics, marginal distributions, and correlations; performs the spatial-block cross-validation and produces prior/posterior predictions. Creates all plots. All routines are contained in "analysis_main.R*, which loads functions defined in the other scripts.

## Folder *model*

Contains the RangeShifter folder structure that is required by RangeShifter even though no file in- or output takes place.

## Folder *results*

Contains the folders to store all calibration and analysis results. Intermediate results of prediction and validation are provided in the sub-folder *predict*.


---

## Session info

The calculations were run under the following setup:

**R version 4.0.4 (2021-02-15)**  
Platform: x86_64-pc-linux-gnu (64-bit)  
Running under: Debian GNU/Linux 11 (bullseye)  

**Matrix products:**  
default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

**Locale:**  
 LC_CTYPE=en_US.UTF-8; LC_NUMERIC=C; LC_TIME=C; LC_COLLATE=en_US.UTF-8; LC_MONETARY=C; LC_MESSAGES=en_US.UTF-8; LC_PAPER=C; LC_NAME=C; LC_ADDRESS=C; LC_TELEPHONE=C; LC_MEASUREMENT=C; LC_IDENTIFICATION=C  

**Attached base packages:**  
 stats; graphics; grDevices; utils; datasets; methods; base  

**Other attached packages:**  
 RangeShiftR_1.0.3; raster_3.5-21; sp_1.5-0  

**Loaded via a namespace (and not attached):**  
 compiler_4.0.4; rgdal_1.5-27; tools_4.0.4; Rcpp_1.0.7; codetools_0.2-18; grid_4.0.4; rbibutils_2.2.8; Rdpack_2.3.1; lattice_0.20-41; terra_1.5-34 
 
