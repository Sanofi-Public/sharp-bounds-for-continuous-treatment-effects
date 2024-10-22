# Sharp Bounds for Continuous-Valued Treatment Effects with Unobserved Confounders

## Abstract

In causal inference, treatment effects are typically estimated under the ignorability, or unconfoundedness, assumption, which is often unrealistic in observational data. By relaxing this assumption and conducting a sensitivity analysis, we introduce novel bounds and derive confidence intervals for the Average Potential Outcome (APO) - a standard metric for evaluating continuous-valued treatment or exposure effects. We demonstrate that these bounds are sharp under a continuous sensitivity model, in the sense that they give the smallest possible interval under this model, and propose a doubly robust version of our estimators. In a comparative analysis with the method of Jesson et al. (2022), using both simulated and real datasets, we show that our approach not only yields sharper bounds but also achieves good coverage of the true APO, with significantly reduced computation times.

## Description
This repository is related to the paper *Sharp Bounds for Continuous-Valued Treatment Effects with Unobserved Confounders*, Baitairian et al. (2024). It is organized as follows:

- **sharp_bounds_cont_val_treatment.Rmd**: a R Markdown document to reproduce the figures from the paper;
- **utils.R**: functions used by the proposed algorithm and the one from Jesson et al. (2022);
- **cont_qb_fun.R**: functions used exclusively by the proposed algorithm;
- **jesson_fun.R**: functions used exclusively by the algorithm from Jesson et al. (2022);
- **simulated_data.fun**: functions used to create the simulated dataset;
- **sens_analys_MC_cmr_and_simu.R**: file to perform sensitivity analyses on Monte-Carlo simulated sample and on real data;
- **MC_simu_gamma_estim.R**: file to perform the exploratory analysis of the evolution of the sensitivity parameter $\Gamma$ with respect to the dataset generation parameters;
- **/images**: a folder to save the generated plots;
- **/params**: a folder to contain optimal neural network parameters obtained after the fine-tuning step;
- **/data**: a folder that contains the real data from the U.S. Environmental Protection Agency (Wyatt et al., 2020).

*Note:* all functions from Jesson et al. (2022) were reimplemented in R, as the original code is in Python. See https://github.com/oatml/overcast.

### References
Jesson, Andrew et al. (2022). “Scalable sensitivity and uncertainty analyses for causal-effect estimates of continuous-valued interventions”. In: Advances in Neural Information Processing Systems 35, pp. 13892–13907.

Wyatt, Lauren H et al. (2020). “Annual PM2. 5 and cardiovascular mortality rate data: Trends modified by county socioeconomic status in 2,132 US counties”. In: Data in brief 30, p. 105318

## Requirements and Licenses
The following libraries are used in this repository with R version 4.3.2.

|Library|Version|License|
|---|---|---|
|`foreach`|1.5.2|Apache License (== 2.0)|
|`ggplot2`|3.4.4|MIT|
|`latex2exp`|0.9.6|MIT|
|`tictoc`|1.2.1|Apache License (== 2.0)|
|`torch`|0.12.0|MIT|
