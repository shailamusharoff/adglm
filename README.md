# adglm

ADGLM code

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
asthma age sex afr eur PC1 PC2 PC3 PC4 PC5
0 18.2 F 0.41 0.59 0.247 0.232 -0.414 0.091 0.111
0 13.9 F 0.51 0.49 -0.04 -0.484 0.265 0.084 -0.342
0 19.4 F 0.42 0.58 0.236 -0.207 0.099 0.035 0.356
0 17.8 F 0.63 0.37 0.014 -0.052 -0.107 -0.313 0.222
0 18.9 M 0.27 0.73 -0.533 -0.536 0.073 -0.127 -0.508
1 19.3 M 0.65 0.35 -0.178 -0.148 0.247 -0.017 -0.064
0 11.9 F 0.71 0.29 0.069 -0.33 0.326 0.151 -0.348
1 20.2 F 0.42 0.58 -0.15 -0.514 -0.176 -0.403 -0.622
0 9.5 F 0.59 0.41 -0.51 -0.234 -0.367 -0.167 -0.332
```

## Output
Three files will be generated. The "_fit.txt" file contains all parameter estimates, standard errors, and Wald test p-values. The Wald test p-values are not to be used for testing because they do not allow for tests that are more than one degree of freedom, and are less powerful than likelihood ratio test (LRT). The "_lrt.txt" contains LRT results; use these p-values. The "_log.txt" is a log file. Samples of these files are provided.

## Questions?

Please [file an issue here](https://github.com/shailam/adglm/issues) on GitHub or email me at shaila [dot] musharoff [at] gmail [dot] com.
