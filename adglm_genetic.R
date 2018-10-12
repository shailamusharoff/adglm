# Note: the dglm package only runs on a matrix with no NAs. the hetglm package can run with NAs.

library(glmx)
library(lmtest)
library(dglm)
library(plyr)
library(dplyr)
library(yaml)
source('adglm_genetic_functions.R')


## Argument parsing -----------------------------------------

NARGS=4                 # number of required arguments to program
NPARAMS = 14            # number of required parameters in file

args = commandArgs(trailingOnly=TRUE)
if((length(args) != NARGS)) {
  print(paste('Got', length(args), 'arguments\n'))
  stop(sprintf("Exactly %d argument(s) must be supplied.n", NARGS), call.=F)
}
param_file = args[1]
outdir = args[2]
ped_file = args[3]              # per-chrom genotype files
map_file = args[4]              # per-chrom genotype files

params = yaml.load_file(param_file)
if((length(params) != NPARAMS)) {
  stop(sprintf("Exactly %d parameters must be supplied. Got: %d.n", NPARAMS, length(params)), call.=F)
}

pheno_name = params$pheno_name         # name of phenotype to test
is_binary = params$is_binary           # boolean specifying test to run. T or F. if T, run logistic regresion on binary phenotype
outbase = params$outbase               # prefix of output file names
pheno_file = params$pheno_file         # data dictionary csv of covariates including pheno, age, ancestry, etc
model_name = params$model_name         # string. printed out in table
pheno_transform = params$pheno_transform      # string: orig | quantnorm | truncate
theta_transform = params$theta_transform      # string: orig | quantnorm | truncate
sd_thresh = as.numeric(params$sd_thresh)      # integer. points more than this many standard deviations from the mean are outliers
dof = params$dof
mean_covar = params$mean_covar
var_covar = params$var_covar
mean_test = params$mean_test
var_test = params$var_test
lrt_name = params$lrt_name


# Process arguments ---------------------
fit_file = paste0(outdir, '/', outbase, '_fit.txt')
lrt_file = paste0(outdir, '/', outbase, '_lrt.txt')
log_file = paste0(outdir, '/', outbase, '.log')
lm_file = paste0(outdir, '/', outbase, '_lm.txt')

valid_transformations = c('orig', 'quantnorm', 'truncate')
if(! pheno_transform %in% valid_transformations) {
  stop(sprintf('pheno_transform: got value %s. Must be one of: %s', pheno_transform, toString(valid_transformations)))
}
if(! theta_transform %in% valid_transformations) {
  stop(sprintf('theta_transform: got value %s. Must be one of: %s', theta_transform, toString(valid_transformations)))
}
cat(param_file,
    '\n\nArguments:', 
    '\n    param_file: ', param_file, 
    '\n    outdir: ', outdir, 
    '\n    fit_file: ', fit_file, 
    '\n    lrt_file: ', lrt_file, 
    '\n    log_file: ', log_file,
    '\n    ped_file: ', ped_file,
    '\n    map_file: ', map_file,    
    '\n    lm_file: ', lm_file,    
    '\n\n', file=log_file, append=F)

## Main from null tables ---------------------------------

# 0. write headers
fit_colnames = c('snp_name', 'chr', 'pos', 'model_name', 'param_type', 'param_name', 'estimate', 'se', 'stat', 'pvalue')
lrt_colnames = c('snp_name', 'chr', 'pos', 'model_name', 'param_type', 'param_name', 'loglik', 'df', 'chisq', 'pvalue')
cat(paste(fit_colnames, sep=' '), '\n', file=fit_file, append=F, sep=' ')
cat(paste(lrt_colnames, sep=' '), '\n', file=lrt_file, append=F, sep=' ')
cat(paste(fit_colnames, sep=' '), '\n', file=lm_file, append=F, sep=' ')

# 1. read in a file of all phenotypes and covariates
#    exclude individuals with NAs because dglm cannot handle any NAs
model_terms = unique(c(mean_covar, mean_test, var_covar, var_test))
model_terms = setdiff(model_terms, c('genotype', '1'))   # genotype not in this dataframe. 1 is for constant variance and not a covariate.
covar = make_covar_df(pheno_name, pheno_file, is_binary, pheno_transform, theta_transform, sd_thresh, model_terms, log_file)

# 2. make dataframe of single phenotype and covariates of interest
num_covar_cols = ncol(covar)             # col for IID, col for pheno, and one col per other covariates

# 3. check if number of phenotype values is compatible with test. will re-print lines to logfile
num_unique_values = categorize_pheno(covar$pheno, log_file)
if(is_binary & (num_unique_values != 2)) {
  stop('Quitting: refusing to run binary test on phenotype with %d values.\n', num_unique_values)
}
if(! is_binary & (num_unique_values == 2)) {
  cat('Warning: running continuous test on binary phenotype with %d values.\n', num_unique_values, file=log_file, append=T)
}

# 1. read genotypes. first col is IID, individual ID
map = read_map(map_file)
genos = read_genos(ped_file, map_file)
num_snps = ncol(genos) - 1

# 3. make dataframe of covar, phenos, genos. Done here for speed: joining each geno in loop below is slow.
all_df = left_join(covar, genos, by='IID')
if(nrow(all_df) != nrow(covar)) {
  print(cat('Warning: after join, number of individuals is', nrow(all_df), ', which is', nrow(covar)-nrow(all_df), 'less than in covar'))
}

for(idx in 1:num_snps) { 
  # make matrix with geno
  geno_idx = idx + num_covar_cols                    # all_df index
  snp_name = colnames(all_df)[geno_idx]
  chr = map %>% filter(ID==snp_name) %>% select(CHR) %>% unlist
  pos = map %>% filter(ID==snp_name) %>% select(BP) %>% unlist
  
  snp_df = all_df[, c(1:num_covar_cols, geno_idx)]   # includes genotypes
  colnames(snp_df)[num_covar_cols + 1] = 'genotype'  # original name is snp name
  na_idx = which(!complete.cases(snp_df))
  snp_df = snp_df[complete.cases(snp_df),]

  # assocations tests
  for(curr_model_name in c('lm', model_name)) {
    if(is_binary) {
      gwas_binary(snp_df, curr_model_name, lrt_name, dof, fit_file, lrt_file, lm_file, snp_name, chr, pos)
    } else {
      gwas_continuous(snp_df, curr_model_name, lrt_name, dof, fit_file, lrt_file, lm_file, snp_name, chr, pos)
    }
  }
}

# print all warnings
warning_msg = toString(warnings())
cat('\nWarnings:\n', warning_msg, '\n', file=log_file, append=T)
