# ADGLM: Ancestry Double Generalized Linear Model

`ADGLM` tests for population variance structure and performs genome-wide association studies (GWAS) and variance quantitative trait loci (vQTL) tests that are corrected for population structure. This is code from the paper: [_Existence and implications of population variance structure_](https://www.biorxiv.org/content/early/2018/10/11/439661).

# Usage
```sh
# Basic usage
Rscript adglm_ancestry.R ancestry_parameters.txt &> logfile.log
```

## Dependencies
```R
install.packages("glmx")
install.packages("lmtest")
install.packages("dglm")
install.packages("plyr")
install.packages("dplyr")
install.packages("yaml")
```

## Options
These are read from a parameter file that is in yaml format (e.g. ancestry_parameters.txt).

```sh
# Name of phenotype to test
pheno_name = params$pheno_name         

# Boolean specifying test to run. T or F. if T, run logistic regresion on binary phenotype
is_binary = params$is_binary           

# Phenotype transform to apply. Options: 'orig' | 'quantnorm' | 'truncate'
pheno_transform = params$pheno_transform  

# Theta transform to apply. Options: 'orig' | 'quantnorm' | 'truncate'
theta_transform = params$theta_transform      

# White-space delimited file of phenotype and covariates including ancestry
pheno_file = params$pheno_file         

# Degrees of freedom of test
dof = as.numeric(params$dof)

# String: printed out
model_name = params$model_name

# Covariates of mean model
mean_covar = params$mean_covar

# Covariates of variance model
var_covar = params$var_covar

# Mean terms for LRT
mean_test = params$mean_test

# Variance terms for LRT
var_test = params$var_test

# Parameter estimates, standard errors, Wald p-values (do not use for testing)
fit_file = params$fit_file

# LRT p-values: use these for testing
lrt_file = params$lrt_file

# Log file with errors and warnings
log_file = params$log_file             

# Float: points more than this many standard deviations from the mean are outliers
sd_thresh = as.numeric(params$sd_thresh)
```

## Input
Phenotype file is a white-space delimited file of phenotype and covariates including ancestry. This simulated example file has the same variable names as the parameter file included:
```sh
IID asthma age sex afr nam eur PC1 PC2 PC3 PC4 PC5
1 0 17.5 M 0.303 0.122 0.574 -0.043 0.053 -0.147 -0.245 -0.088
2 0 20.2 M 0.028 0.59 0.383 -0.021 -0.364 -0.53 -0.3 -0.231
3 0 18.7 F 0.018 0.799 0.183 -0.331 0.184 0.146 0.084 -0.226
4 0 13 F 0.251 0.136 0.613 -0.357 -0.286 -0.585 -0.177 -0.161
5 0 20.6 M 0.334 0.092 0.574 0.214 -0.519 -0.159 0.105 -0.386
6 0 11.3 M 0.214 0.132 0.654 -0.09 0.034 -0.059 -0.255 0.219
7 0 18 F 0.038 0.498 0.465 -0.287 -0.236 0.703 -0.257 0.202
8 0 8.6 M 0.179 0.494 0.328 -0.512 0.496 -0.267 -0.14 0.331
```

## Output
Three files will be generated. The "_fit.txt" file contains all parameter estimates, standard errors, and Wald test p-values. The Wald test p-values are not to be used for testing because they do not allow for tests that are more than one degree of freedom, and are less powerful than likelihood ratio test (LRT). The "_lrt.txt" contains LRT results; use these p-values. The "_log.txt" is a log file. Samples of these files are provided.


## GWAS and vQLT
This code takes in as input plink map and ped files that have been recoded with the `--recode 12` option. The phenotype file must have an "IID" column that matches that of the plink files.

```sh
# GWAS
param_file=fev1_betag_params.txt
outdir=gwas_outdir
map_file=test_recode12.map
ped_file=test_recode12.ped
Rscript adglm_genetic.R ${param_file} ${outdir} ${ped_file} ${map_file}

# vQTL
param_file=fev1_sigmag_params.txt
outdir=vqlt_outdir
map_file=test_recode12.map
ped_file=test_recode12.ped
Rscript adglm_genetic.R ${param_file} ${outdir} ${ped_file} ${map_file}
```

## Questions?

Please [file an issue here](https://github.com/shailam/adglm/issues) on GitHub or email me at shaila [dot] musharoff [at] gmail [dot] com.
