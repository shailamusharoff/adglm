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

# Data dictionary csv of covariates including pheno, age, ancestry, etc
pheno_file = params$pheno_file         

# File with genetic PCs
pc_file = params$pc_file               

# Individuals to analyze. one line per individual with two columns, no header
idv_file = params$idv_file             

# Degrees of freedom of test
dof = as.numeric(params$dof)

# String. printed out in table
model_name = params$model_name

# Covariates of mean model
mean_covar = params$mean_covar

# Covariates of variance model
var_covar = params$var_covar

# Mean terms for LRT
mean_test = params$mean_test

# Variance terms for LRT
var_test = params$var_test

# Parameter estimates, standard errors, Wald p-values (do not use)
fit_file = params$fit_file

# LRT p-values: use these for testing
lrt_file = params$lrt_file

# Log file with errors and warnings
log_file = params$log_file             

# Integer. points more than this many standard deviations from the mean are outliers
sd_thresh = as.numeric(params$sd_thresh)
```

## Questions?

Please [file an issue here](https://github.com/shailam/adglm/issues) on GitHub or email me at shaila [dot] musharoff [at] gmail [dot] com.
