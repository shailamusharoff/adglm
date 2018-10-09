library(dglm)
library(glmx)
library(lmtest)
library(plyr)
library(dplyr)

# For discrete or continuous phenotypes

## Functions -----------------------------------------
fit_colnames = c('test', 'model', 'param_type', 'param_name', 'estimate', 'se', 'stat', 'pvalue')
header = 'test model param_type param_name estimate se stat pvalue'

arg_to_boolean <- function(x) {
  ifelse( (!is.na(x) && (x == 'T')), T, F) 
}

is.error <- function(x) {
  inherits(x, "try-error")
}

quant_norm <- function(m) {
  # returns values quantile-normalized to a standard normal distribution
  probs = ppoints(sort(m))        # sequence of probability points
  normQuant = qnorm(probs)        # quantiles of a standard normal distribution
  m_norm = m
  m_norm[order(m_norm)] = normQuant
  if(! all.equal(order(m), order(m_norm))) {  # check if min and max indices match for original and transformed values
    stop('Error: quantile normalization did not work')
  }
  return(m_norm)
}

write_output <- function(fit_df, outfile) {
  colnames(fit_df) = fit_colnames
  write.table(fit_df, outfile, row.names=F, col.names=F, quote=F, append=T)
}

run_test_continuous <- function(covar_df, model_name, outfile, PC_name=NA) {
  if(! model_name %in% c('lm', 'adglm_pc', 'adglm_sex', 'adglm_age')) {
    stop('Invalid adglm model name')
  }

  if(model_name == 'lm') {
    fit = try(lm(pheno ~ cov_ASSESS_CENTER + cov_GENO_ARRAY + cov_SEX + cov_AGE + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covar_df))
  } else {
    fit_null = try(dglm(pheno ~ cov_ASSESS_CENTER + cov_GENO_ARRAY + cov_SEX + cov_AGE + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covar_df, family='gaussian'))
    if(model_name == 'adglm_pc') {
      if(is.na(PC_name))
        stop('Missing PC_name')
      covar_df = covar_df %>% mutate_(var = as.name(PC_name))
    } else if(model_name == 'adglm_sex') {
      covar_df = covar_df %>% mutate(var = cov_SEX)
    } else if(model_name == 'adglm_age') {
      covar_df = covar_df %>% mutate(var = cov_AGE)
    }
    fit_alt = try(dglm(pheno ~ cov_ASSESS_CENTER + cov_GENO_ARRAY + cov_SEX + cov_AGE + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, ~var, data = covar_df, family='gaussian'))
  }
  if(model_name == 'lm') {
    if(! is.error(fit)) {
      mean_fit = summary(fit)$coefficients
      fit_df = cbind.data.frame(model_name, 'null', 'mean', rownames(mean_fit), mean_fit)
      write_output(fit_df, outfile)
    }
  } else {
    if((! is.error(fit_null)) & (! is.error(fit_alt))) {
      if(model_name == 'adglm_sex') {     # print out once because same for all models. should be identical to lm output
        null_fit = summary(fit_null)$coefficients
        null_df = cbind.data.frame('adglm', 'null', 'mean', rownames(null_fit), null_fit)
        write_output(null_df, outfile)
      }
      if(model_name == 'adglm_pc') {
        model_name = paste('adglm', PC_name, sep='_')
      }
      chiStat = fit_null$m2loglik - fit_alt$m2loglik
      pval = pchisq(chiStat, df=1, lower.tail=F)
      # pval = exp(pchisq(chiStat, df=1, lower.tail=F, log.p=T))      # alternate for small stat
      res = c(-2*fit_alt$m2loglik, 1, chiStat, pval)     # NLL_alt df lambda_stat pval
      mean_fit = summary(fit_alt)$coefficients
      mean_df = cbind.data.frame(model_name, 'alt', 'mean', rownames(mean_fit), mean_fit)
      var_fit = summary(fit_alt)$dispersion.summary$coefficients
      var_df = cbind.data.frame(model_name, 'alt', 'var', rownames(var_fit), var_fit)
      lrt_df = data.frame(t(c(model_name, 'LRT', 'var', rownames(var_fit), res)))
      
      write_output(lrt_df, outfile)
      write_output(var_df, outfile)
      write_output(mean_df, outfile)
    }
  }
}


run_test_discrete <- function(covar_df, model_name, outfile, PC_name=NA) {
  if(! model_name %in% c('lm', 'adglm_pc', 'adglm_sex', 'adglm_age')) {
    stop('Invalid adglm model name')
  }
  
  if(model_name == 'lm') {
    fit = try(glm(pheno ~ cov_ASSESS_CENTER + cov_GENO_ARRAY + cov_SEX + cov_AGE + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = covar_df, family = binomial(link = "probit")))  # works
  } else {
    fit_null = try(hetglm(pheno ~ cov_ASSESS_CENTER + cov_GENO_ARRAY + cov_SEX + cov_AGE + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | 1, data = covar_df))  # same as fit above
    if(model_name == 'adglm_pc') {
      covar_df = covar_df %>% mutate_(var = as.name(PC_name))
    } else if(model_name == 'adglm_sex') {
      covar_df = covar_df %>% mutate(var = cov_SEX)
    } else if(model_name == 'adglm_age') {
      covar_df = covar_df %>% mutate(var = cov_AGE)
    }
    fit_alt = try(hetglm(pheno ~ cov_ASSESS_CENTER + cov_GENO_ARRAY + cov_SEX + cov_AGE + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 | var, data = covar_df))
  }
  
  if(model_name == 'lm') {
    if(! is.error(fit)) {
      mean_fit = summary(fit)$coefficients   # for glm call only. hetglm call needs $mean
      fit_df = cbind.data.frame(model_name, 'null', 'mean', rownames(mean_fit), mean_fit)
      write_output(fit_df, outfile)
    }
  } else {
    if((! is.error(fit_null)) & (! is.error(fit_alt))) {    # hetglm
      if(model_name == 'adglm_sex') {     # print out once because same for all models. should be identical to glm output
        null_fit = summary(fit_null)$coef$mean
        null_df = cbind.data.frame('adglm', 'null', 'mean', rownames(null_fit), null_fit)
        write_output(null_df, outfile)
      }
      if(model_name == 'adglm_pc') {
        model_name = paste('adglm', PC_name, sep='_')
      }
      res = lrtest(fit_null, fit_alt)
      mean_fit = summary(fit_alt)$coef$mean    # need this: fit under alt model
      mean_df = cbind.data.frame(model_name, 'alt', 'mean', rownames(mean_fit), mean_fit)
      var_fit = summary(fit_alt)$coef$scale
      var_df = cbind.data.frame(model_name, 'alt', 'var', rownames(var_fit), var_fit)
      lrt_df = data.frame(t(c(model_name, 'LRT', 'var', rownames(var_fit), res[2,2], res[2,3], res[2,4], res[2,5])))
      
      write_output(lrt_df, outfile)
      write_output(var_df, outfile)
      write_output(mean_df, outfile)
    }
  }
}


## Argument parsing -----------------------------------------
args = commandArgs(trailingOnly=TRUE)
if((length(args) != 7)) {
  print(paste('Got', length(args), 'arguments'))
  stop("Exactly 7 arguments must be supplied.n", call.=F)
}

pheno_name = args[1]            # name of phenotype to test
quantnorm_pheno_arg = args[2]   # T or F. Must be F for discrete pheno
pc_file = args[3]               # file with covariates including PCs, age, sex
pheno_file = args[4]            # file with phenotypes
discrete_arg = args[5]          # T or F. if T, run boolean tests.
outdir = args[6]                # output directory
diagnostic_arg = args[7]        # T or F. if T, estimates correlations and makes plots.

quantnorm_pheno = arg_to_boolean(quantnorm_pheno_arg)
discrete = arg_to_boolean(discrete_arg)
diagnostic = arg_to_boolean(diagnostic_arg)

cat('\nArguments:\n', pheno_name, quantnorm_pheno, pc_file, pheno_file, discrete, outdir, diagnostic, '\n\n')


## Read input files --------------------------------------------
pc_df = read.table(pc_file, header=T, stringsAsFactors=F)
pheno_df = read.table(pheno_file, header=T, stringsAsFactors=F)
all_df = merge(pc_df, pheno_df, by=c('FID',	'IID'))

# subset to include phenotype and covariates of interest
colnames(all_df)[which(colnames(all_df) == pheno_name)] = 'pheno'
covar_df = all_df[, c('FID', 'IID', 'pheno', 'cov_ASSESS_CENTER', 'cov_GENO_ARRAY', 'cov_SEX', 'cov_AGE', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')]

if(ncol(covar_df) != 17) {
  stop('Error: did not pull out single phenotype')
}

# exclude any individuals with NA measurements
covar_df = covar_df[complete.cases(covar_df), ]
if(quantnorm_pheno) {
  covar_df$pheno = quant_norm(covar_df$pheno)
} 

if(length(unique(covar_df$pheno)) == 2) {   # there should not be any NA's here
  print('Discrete pheno: has 2 values')
  if(quantnorm_pheno) {
    cat('Error: do not quantnorm discrete phenotype\n')
  }
}

# open file and write header
outfile = paste(outdir, '/', pheno_name, '_qnorm', quantnorm_pheno_arg, '.txt', sep='')  
cat(header, file=outfile, append=F, sep='\n')

# run all tests and append output lines of each test to a growing file

if(discrete)  {
  run_test_discrete(covar_df, 'lm', outfile)
  run_test_discrete(covar_df, 'adglm_sex', outfile)
  run_test_discrete(covar_df, 'adglm_age', outfile)
  for(PC_name in c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')) {
    run_test_discrete(covar_df, 'adglm_pc', outfile, PC_name)
  }
} else {  # continuous
  run_test_continuous(covar_df, 'lm', outfile)
  run_test_continuous(covar_df, 'adglm_sex', outfile)
  run_test_continuous(covar_df, 'adglm_age', outfile)
  for(PC_name in c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')) {
    run_test_continuous(covar_df, 'adglm_pc', outfile, PC_name)
  }
}

if(diagnostic) {
  # correlation of phenotype and covariate
  outfile = paste(outdir, '/', pheno_name, '_qnorm', quantnorm_pheno_arg, '.out', sep='')  
  covars = c('cov_SEX', 'cov_AGE', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10')
  cor_line = data.frame(t(c('pheno_cor', sapply(covars, function(x) cor(covar_df[,'pheno'], covar_df[, x])))))
  write.table(cor_line, outfile, row.names=F, quote=F, sep=' ')
  
  # single pdf of plots
  outfile = paste(outdir, '/', pheno_name, '_qnorm', quantnorm_pheno_arg, '.pdf', sep='')
  pdf(outfile)
  hist(covar_df[,'pheno'], xlab='y', main=pheno_name, col='gray')
  for(covar in covars) {
    hist(covar_df[,covar], xlab='y', main=covar, col='gray')
  }
  for(covar in covars) {
    plot(covar_df[,'pheno'], covar_df[,covar], xlab='y', ylab=covar, main=pheno_name)
  }
  dev.off()
}