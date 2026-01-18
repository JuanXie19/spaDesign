# spaDesign: Simulation-based Experimental Design for Sequencing-based Spatial Transcriptomics Studies

[![R-CMD-check](https://github.com/JuanXie19/spaDesign/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/JuanXie19/spaDesign/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

`spaDesign` is a statistical framework for optimizing sequencing-based spatial transcriptomics experiments.

It helps researchers answer the critical question: **"How much sequencing depth do I actually need?"** by quantifying how sequencing depth and signal strength (e.g., domain effect size) impact domain recovery performance (NMI/ARI). It provides tools to detect the **saturation depth**—the point where additional sequencing yields diminishing returns.

### Workflow
The framework follows these steps:
1. **feature selection:** Select domain-informative genes
1.  **Model:** Learn spatial parameters from a pilot or public dataset (NNGP-based).
2.  **Simulate:** Generate new count matrices under varying sequencing depths and effect sizes/spatial disturbance.
3.  **Evaluate:** Measure domain recovery (NMI) to identify the saturation point.

---

## Installation

You can install the development version of `spaDesign` from GitHub:

```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("JuanXie19/spaDesign")
```
## Quick Start
library(spaDesign)

### 1. create a `spaDesign` object
```r
DATA <- createDesignObject(count_matrix = count_matrix, loc = loc_file)
```

### 2. Select domain-informative features
```r
DATA <- featureSelection(object = DATA, logfc_cutoff = 0.7,
                              mean_in_cutoff = 1.8,
                              max_num_gene = 10,
                              n_cores = 4)
```
### 3. Parameter estimation
```r
# learning spatial distribution pattern 
DATA <- estimation_FGEM(DATA, iter_max = 1000,
                        M_candidates = 2:6,
                        tol = 0.1,
                        n_cores = 4, verbose = F)
# learning spatial expression pattern
DATA <- estimation_NNGP(DATA, n_neighbors = 10,
                        order = 'AMMD',
                        X = NULL, verbose = FALSE)
```
### 4. Simulate sequencing depth gradients and calculate performance
```r
# simulate across defined sequencing depth
res <- powerAnalysisEffectSize(
    DATA, es_range = 1,
    seq_depth_range = 1:5,
    n_rep = 3,
    n_cores = 4
)
```
### 5. Detect saturation and plot results
```r
saturation_point <- saturationDetection(res)
plotSaturation(saturation_point)
```
## Vignette
