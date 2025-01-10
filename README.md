# Sharp Bounds for Continuous-Valued Treatment Effects with Unobserved Confounders

## Abstract [(Link to the article)](https://arxiv.org/abs/2411.02231)

In causal inference, treatment effects are typically estimated under the ignorability, or unconfoundedness, assumption, which is often unrealistic in observational data. By relaxing this assumption and conducting a sensitivity analysis, we introduce novel bounds and derive confidence intervals for the Average Potential Outcome (APO) - a standard metric for evaluating continuous-valued treatment or exposure effects. We demonstrate that these bounds are sharp under a continuous sensitivity model, in the sense that they give the smallest possible interval under this model, and propose a doubly robust version of our estimators. In a comparative analysis with the method of Jesson et al. (2022), using both simulated and real datasets, we show that our approach not only yields sharper bounds but also achieves good coverage of the true APO, with significantly reduced computation times.

## Description
This repository is related to the paper *Sharp Bounds for Continuous-Valued Treatment Effects with Unobserved Confounders*, Baitairian et al. (2024). It is organized as follows:

- **sharp_bounds_cont_val_treatment.Rmd**: a R Markdown document to reproduce the figures and numerical results from the paper. Users only need to execute code from this file;
- **utils.R**: functions used by the proposed algorithm and the one from Jesson et al. (2022);
- **cont_qb_fun.R**: functions used exclusively by the proposed algorithm;
- **jesson_fun.R**: functions used exclusively by the algorithm from Jesson et al. (2022);
- **simulated_data.fun**: functions used to create the simulated dataset;
- **sens_analys_MC_cmr_and_simu.R**: file to perform sensitivity analyses on Monte-Carlo simulated sample and on real data;
- **MC_simu_gamma_estim.R**: file to perform the exploratory analysis of the evolution of the sensitivity parameter $\Gamma$ with respect to the dataset generation parameters;
- **sens_param_estim_inf_bench.R**: file to estimate plausible values of the sensitivity parameter $\Gamma$ on the real data via informal benchmarking;
- **/images**: a folder to save the generated plots;
- **/params**: a folder to contain optimal neural network parameters obtained after the fine-tuning step;
- **/data**: a folder that contains the real data from the U.S. Environmental Protection Agency (Wyatt et al., 2020).

*Note:* all functions from Jesson et al. (2022) were reimplemented in R, as the original code is in Python. See https://github.com/oatml/overcast.

### References
Jesson, Andrew et al. (2022). “Scalable sensitivity and uncertainty analyses for causal-effect estimates of continuous-valued interventions”. In: Advances in Neural Information Processing Systems 35, pp. 13892–13907.

Wyatt, Lauren H et al. (2020). “Annual PM2. 5 and cardiovascular mortality rate data: Trends modified by county socioeconomic status in 2,132 US counties”. In: Data in brief 30, p. 105318

## Requirements and Licenses
This repository is under a non commercial license (see **LICENSE.txt** file). The following libraries are used in this repository with R version 4.3.2.

|Library|Version|License|
|---|---|---|
|`foreach`|1.5.2|Apache License (== 2.0)|
|`ggplot2`|3.4.4|MIT|
|`latex2exp`|0.9.6|MIT|
|`tictoc`|1.2.1|Apache License (== 2.0)|
|`torch`|0.12.0|MIT|

## Output of the `sessionInfo()` command

```
R version 4.3.2 (2023-10-31)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Rocky Linux 8.7 (Green Obsidian)

Matrix products: default
BLAS/LAPACK: /usr/lib64/libopenblas-r0.3.15.so;  LAPACK version 3.9.0

locale:
 [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8        LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
 [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C           LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   

time zone: UTC
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] tictoc_1.2.1    latex2exp_0.9.6 torch_0.12.0    foreach_1.5.2   ggplot2_3.4.4  

loaded via a namespace (and not attached):
 [1] bit_4.0.5         gtable_0.3.4      dplyr_1.1.4       compiler_4.3.2    tidyselect_1.2.0  Rcpp_1.0.11       stringr_1.5.1     parallel_4.3.2   
 [9] callr_3.7.3       scales_1.3.0      yaml_2.3.8        fastmap_1.1.1     R6_2.5.1          doSNOW_1.0.20     generics_0.1.3    knitr_1.45       
[17] iterators_1.0.14  tibble_3.2.1      snow_0.4-4        munsell_0.5.0     pillar_1.9.0      rlang_1.1.2       utf8_1.2.4        stringi_1.8.3    
[25] xfun_0.41         bit64_4.0.5       cli_3.6.2         withr_2.5.2       magrittr_2.0.3    ps_1.7.5          digest_0.6.34     grid_4.3.2       
[33] processx_3.8.3    rstudioapi_0.15.0 lifecycle_1.0.4   coro_1.0.3        vctrs_0.6.5       evaluate_0.23     glue_1.6.2        codetools_0.2-19 
[41] fansi_1.0.6       colorspace_2.1-0  rmarkdown_2.25    tools_4.3.2       pkgconfig_2.0.3   htmltools_0.5.7
```
