## Runs null tests for a batch of populations and phenotypes
#  Output: for each population, ancestry variance estimates, along with ordinary linear or logistic regression estimates in "null mean" rows of dataframe
# Note: the dglm package only runs on a matrix with no NAs. the hetglm package can run with NAs.


# load libraries
library(glmx)
library(lmtest)
library(dglm)
library(plyr)
library(dplyr)
library(yaml)

source('adglm_ancestry_functions.R')

## Argument parsing -----------------------------------------

NARGS=1                 # number of required arguments to program
NPARAMS = 14            # number of required parameters in file

terse = T               # if T, writes one line per phenotype
args = commandArgs(trailingOnly=TRUE)
if((length(args) < NARGS)) {
  print(paste('Got', length(args), 'arguments\n'))
  stop(sprintf("At least %d argument must be supplied.n", NARGS), call.=F)
}
param_file = args[1]
params = yaml.load_file(param_file)
if((length(params) < NPARAMS)) {
  print(paste('Parameter file has', length(param), 'entries'))
  stop(sprintf("At least %d argument must be supplied.n", NPARAMS), call.=F)
}

pheno_name = params$pheno_name         # name of phenotype to test
is_binary = params$is_binary           # boolean specifying test to run. T or F. if T, run logistic regresion on binary phenotype
pheno_transform = params$pheno_transform  # string: orig | quantnorm | truncate
theta_transform = params$theta_transform      # string: orig | quantnorm | truncate

pheno_file = params$pheno_file         # data dictionary csv of covariates including pheno, age, ancestry, etc
pc_file = params$pc_file               # file with genetic PCs
idv_file = params$idv_file             # individuals to analyze. one line per individual with two columns, no header

dof = as.numeric(params$dof)
model_name = params$model_name         # string. printed out in table
mean_covar = params$mean_covar
var_covar = params$var_covar
mean_test = params$mean_test
var_test = params$var_test

fit_file = params$fit_file             # output
lrt_file = params$lrt_file             # output
log_file = params$log_file             # log file with errors and warnings
sd_thresh = as.numeric(params$sd_thresh)  # integer. points more than this many standard deviations from the mean are outliers

# process arguments
if(idv_file == 'NA') {
  idv_file=NA    
}
valid_transformations = c('orig', 'quantnorm', 'truncate')
if(! pheno_transform %in% valid_transformations) {
  stop(sprintf('pheno_transform: got value %s. Must be one of: %s', pheno_transform, toString(valid_transformations)))
}
if(is_binary & (pheno_transform != 'orig')) {
  next
}
if(! theta_transform %in% valid_transformations) {
  stop(sprintf('theta_transform: got value %s. Must be one of: %s', theta_transform, toString(valid_transformations)))
}

cat(param_file,
    '\n\nArguments:', 
    '\n    pheno_name: ', pheno_name, 
    '\n    is_binary: ', is_binary, 
    '\n    pheno_file: ', pheno_file, 
    '\n    fit_file: ', fit_file, 
    '\n    lrt_file: ', lrt_file, 
    '\n    log_file: ', log_file, 
    '\n    pc_file: ', pc_file, 
    '\n    idv_file: ', idv_file,
    '\n    pheno_transform: ', pheno_transform,
    '\n    theta_transform: ', theta_transform, 
    '\n    dof: ', dof, 
    '\n    sd_thresh: ', sd_thresh,
    '\n\n', file=log_file, append=F)

## Main
# 0. write headers
if(terse) {
  lrt_colnames = c('pheno_name', 'model_name', 'pheno_trans', 'theta_trans', 'll_null', 'll_alt', 'lrt_stat', 'df', 'lrt_pvalue')
  fit_colnames = c('pheno_name', 'model_name', 'pheno_trans', 'theta_trans', 'hypothesis', 'param_type', 'param_name', 'estimate', 'se', 'stat', 'pvalue')
}

cat(paste(fit_colnames, sep=' '), '\n', file=fit_file, append=F, sep=' ')
cat(paste(lrt_colnames, sep=' '), '\n', file=lrt_file, append=F, sep=' ')

# 1. read in a file of all phenotypes and covariates
covar = make_covar_df(pheno_name, pheno_file, pc_file, idv_file, is_binary, pheno_transform, theta_transform, sd_thresh, log_file)

# 2. exclude individuals with NAs because dglm cannot handle any NAs
model_terms = unique(c(mean_covar, mean_test, var_covar, var_test))
covar = subset_covar_df(covar, model_terms)
retval = make_models(mean_covar, var_covar, mean_test, var_test, dof, log_file)
mean_null_model = retval$mean_null_model
mean_alt_model = retval$mean_alt_model
var_null_model = retval$var_null_model
var_alt_model = retval$var_alt_model

# 3. check if number of phenotype values is compatible with test. will re-print lines to logfile
num_unique_values = categorize_pheno(covar$pheno, log_file, do_print=F)
if(is_binary & (num_unique_values != 2)) {
  stop('Quitting: refusing to run binary test on phenotype with %d values.\n', num_unique_values)
}
if(! is_binary & (num_unique_values == 2)) {
  cat('Warning: running continuous test on binary phenotype with %d values.\n', num_unique_values, file=log_file, append=T)
}

# 4. run models
if(is_binary) {
  retval = run_models_discrete(pheno_name, model_name, mean_null_model, mean_alt_model, var_null_model, var_alt_model, dof, covar, fit_file, lrt_file, log_file, terse)
} else {
  retval = run_models_continuous(pheno_name, model_name, mean_null_model, mean_alt_model, var_null_model, var_alt_model, dof, covar, fit_file, lrt_file, log_file, terse)
}
lrt_df = retval$lrt_df
write.table(lrt_df, lrt_file, row.names=F, col.names=F, quote=F, append=T)
write.table(retval$fit_df, fit_file, row.names=F, col.names=F, quote=F, append=T)

# print all warnings
warning_msg = toString(warnings())
cat('\n\nWarnings:\n', warning_msg, '\n', file=log_file, append=T)
