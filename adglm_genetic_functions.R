
## Phenotype, covariate, and file processing --------------------------------------

is.error <- function(x) {
  inherits(x, "try-error")
}

quant_norm <- function(m) {
  # returns values quantile-normalized to a standard normal distribution
  
  probs = ppoints(sort(m))        # sequence of probability points
  normQuant = qnorm(probs)        # quantiles of a standard normal distribution
  m_norm = m
  m_norm[order(m_norm)] = normQuant
  if(! all.equal(order(m), order(m_norm))) {  # check that orders of original and transformed values match
    stop('Error: quantile normalization did not work')
  }
  return(m_norm)
}

qnorm_norm_with_NA <- function(x, debug=F) {
  # x: numerical vector that can have NAs
  # typical usage: x is a phenotype or covariate
  
  num_unique_values = length(unique(x[!is.na(x)]))
  if(num_unique_values < 5) {
    stop('Quitting: refusing to quantile normalize a vector with less than 5 unique values\n')
  } 
  num_na_before = sum(is.na(x))
  non_na_idx = which(!is.na(x))
  x[non_na_idx] = quant_norm(x[non_na_idx])
  num_na_after = sum(is.na(x))
  
  if(num_na_before != num_na_after) {
    cat(sprintf('Error: quantnorm vector x with NA did not work. Before: %d NAs. After: %d NAs\n', num_na_before, num_na_after), file=log_file, append=T)
  }
  return(x)
}

remove_outliers <- function(x, m, debug=F) {
  # sets values of x more than m standard deviations away from the mean to NA
  # typical usage: x is a phenotype or covariate
  # x: numerical vector that can have NAs
  # m: float
  
  num_unique_values = length(unique(x[!is.na(x)]))
  if(num_unique_values < 5) {
    stop('Quitting: refusing to remove outliers from a vector with less than 5 unique values\n')
  } 
  
  x_bar = mean(x, na.rm=T)
  sd_x = sd(x, na.rm=T)
  outlier_idx = which( abs(x - x_bar) > (m * sd_x) )
  x[outlier_idx] = NA
  
  if(debug) cat(sprintf('Removed %d outliers\n', length(outlier_idx)))
  return(x)
}

categorize_pheno <- function(pheno, log_file) {
  # returns number of unique phenotype values  
  # is_binary: determines whether logistic or continuous test is run
  
  num_unique_values = length(unique(pheno[!is.na(pheno)]))
  if(num_unique_values == 2) {
    cat('Binary phenotype: has 2 values\n', file=log_file, append=T)
  } else if(num_unique_values <= 10) {
    cat(sprintf('Categorical phenotype or binary phenotype with unidentified null characters: has %d values\n', num_unique_values), file=log_file, append=T)
  } else {
    cat(sprintf('Continuous phenotype or categorical phenotype with more than ten values: has %d values\n', num_unique_values), file=log_file, append=T)
  }
  return(num_unique_values)
}

make_covar_df <- function(pheno_name, pheno_file, is_binary, pheno_transform, theta_transform, sd_thresh, model_terms, log_file) {
  # Read in covar file and re-name tested pheno as "pheno" exclude individuals with NAs in tested phenotype, mean model terms, or variance model terms
  # Input:  pheno_name is a column in file
  # Output: dataframe with a column named pheno
  
  # read file of covariates
  all_df = read.table(pheno_file, header=T, stringsAsFactors=F, na.strings=c('NA'))
  
  # get phenotype and transform if indicated
  pheno = all_df %>% select_(as.name(pheno_name))
  colnames(pheno) = 'pheno'
  if(ncol(pheno) != 1) {
    stop(cat('Did not pull out one column based on pheno_name. Got ', ncol(pheno), ' columns.'))
  }
  num_unique_values = categorize_pheno(pheno, log_file)
  if(pheno_transform == 'quantnorm') {
    pheno[,1] = qnorm_norm_with_NA(pheno[,1])    
  } else if(pheno_transform == 'truncate') {
    pheno[,1] = remove_outliers(pheno[,1], sd_thresh)
  } 
  
  # transform theta if indicated
  if(theta_transform == 'quantnorm') {
    all_df$afr = qnorm_norm_with_NA(all_df$afr)
    all_df$eur = qnorm_norm_with_NA(all_df$eur)
    all_df$nam = qnorm_norm_with_NA(all_df$nam)
  } else if(theta_transform == 'truncate') {
    all_df$afr = remove_outliers(all_df$afr, sd_thresh)
    all_df$eur = remove_outliers(all_df$eur, sd_thresh)
    all_df$nam = remove_outliers(all_df$nam, sd_thresh)
  }
  
  covar = cbind.data.frame(all_df, pheno)     # add column named pheno
  cat('Covariate dataframe has', nrow(covar), 'individuals.\n', file=log_file, append=T)
  covar = covar[,c('IID', 'pheno', model_terms)]
  num_incomplete = sum(!complete.cases(covar))
  if(num_incomplete > 0) {
    covar = covar[complete.cases(covar),]
    cat('Removed', num_incomplete, 'individuals with incomplete phenotype and/or covariate measurements.\n', file=log_file, append=T)
    cat('Final dataframe with no NAs has', nrow(covar), 'individuals.\n', file=log_file, append=T)
  }
  return(covar)   
}


## Genotype functions ----------------------------------

read_map <- function(map_file) {
  map = read.table(map_file, header=F, stringsAsFactors=F)
  colnames(map) = c('CHR', 'ID', 'CM', 'BP')
  return(map)
}

read_genos <- function(ped_file, map_file) {
  # reads in a ped file and converts to a genotype matrix with correct rsids from map file
  # cols 1-6 are front matter; haplo1 for snp1 is in col7 and haplo2 for snp1 is in col8
  ped = read.table(ped_file, header=F, stringsAsFactors=F)
  map = read.table(map_file, header=F, stringsAsFactors=F)
  colnames(map) = c('CHR', 'ID', 'CM', 'BP')
  
  h1 = ped[, seq(7, (ncol(ped)-1), 2) ]
  h2 = ped[, seq(8, ncol(ped), 2) ]
  if(ncol(ped) == 8) {   # edge case: one SNP in file
    h1 = as.data.frame(h1)
    h2 = as.data.frame(h2)
  }
  if(ncol(h1) + ncol(h2) + 6 != ncol(ped)) print('error in dims')
  genos = h1 + h2 - 2
  genos[genos==-2] = NA
  genos = cbind.data.frame(ped[,1], genos, stringsAsFactors=F)
  colnames(genos) = c('IID', map$ID)    
  return(genos)
}

## Model fitting and testing ------------------------------------------------------

variance_test_binary <- function(pheno_name, model_name, mean_terms, var_terms, covar, fit_file, lrt_file, log_file) {
  mean_model = paste(mean_terms, collapse=' + ')
  var_model = paste(var_terms, collapse=' + ')
  
  mean_formula = as.formula(paste('pheno~', mean_model))
  het_formula = as.formula(paste('pheno~', mean_model, ' | ' ,var_model))   # hetglm package syntax
  
  fit_null = glm(mean_formula, data = covar, family = binomial(link = "probit"))
  fit_alt = hetglm(het_formula, data = covar)
  
  # TODO use hetglm convergence code below
  if(! fit_alt$converged) {
    cat(sprintf('Warning: hetglm did not converge for phenotype %s on model %s with mean parameters %s and variance parameters %s.\n', pheno_name, model_name, toString(mean_terms), toString(var_terms)), file=log_file, append=T)
  } else {
      cat(sprintf('\nSuccess: hetglm converged for phenotype %s on model %s with mean parameters %s and variance parameters %s.\n', pheno_name, model_name, toString(mean_terms), toString(var_terms)), file=log_file, append=T)
  }
  summarize_fits_binary(pheno_name, model_name, fit_null, fit_alt, fit_file)
  summarize_lrt_binary(pheno_name, model_name, fit_null, fit_alt, lrt_file)
}


variance_test_continuous <- function(pheno_name, model_name, mean_terms, var_terms, covar, fit_file, lrt_file, log_file) {
  # TODO figure out how to check for dglm convergence
  
  if(any(is.na(covar))) {
    cat('Error: covariate matrix has NAs and package dglm will fail. Run function subset_covar_df\n', file=log_file, append=T)
  }
  
  mean_model = paste(mean_terms, collapse=' + ')
  var_model = paste(var_terms, collapse=' + ')
  
  mean_formula = as.formula(paste('pheno~', mean_model))
  var_formula = as.formula(paste('~', var_model))   # dglm package syntax
  
  fit_null = lm(mean_formula, data = covar)
  fit_alt = dglm(mean_formula, var_formula, data = covar)
  
  num_iter_ran = fit_alt$iter
  max_iter = dglm.control()$maxit
  if(num_iter_ran >= max_iter) {
    cat(sprintf('Warning: dglm ran for %d iterations and reached maximum number of iterations %d for phenotype %s on model %s with mean parameters %s and variance parameters %s.\n', num_iter_ran, max_iter, pheno_name, model_name, toString(mean_terms), toString(var_terms)), file=log_file, append=T)
  } else {
    cat(sprintf('Ok: dglm ran for %d iterations which is less than maximum number of iterations %d for phenotype %s on model %s with mean parameters %s and variance parameters %s.\n', num_iter_ran, max_iter, pheno_name, model_name, toString(mean_terms), toString(var_terms)), file=log_file, append=T)
  }
  
  summarize_fits_continuous(pheno_name, model_name, fit_null, fit_alt, fit_file)
  summarize_lrt_continuous(pheno_name, model_name, fit_null, fit_alt, lrt_file)
}



## Writing output: variance tests -------------------------------------------------

write_output <- function(fit_df, colnames_df, outfile) {
  # append to a growing file
  
  colnames(fit_df) = colnames_df   # not written below. just to suppress warning?
  write.table(fit_df, outfile, row.names=F, col.names=F, quote=F, append=T)
}

# variance table
summarize_fits_binary <- function(pheno_name, model_name, fit_null, fit_alt, fit_file) {
  # write output
  # 1. table of model fits
  # 2. TODO R^2, adj R^2, AIC, BIC, likelihood
  
  tmp_df = cbind.data.frame(pheno_name, model_name)       # avoids warnings from cbind.data.frame
  rownames(tmp_df) = c()
  
  null_fit = summary(fit_null)$coef      # output of function glm
  null_df = cbind.data.frame(tmp_df, 'null', 'mean', rownames(null_fit), null_fit)
  
  mean_fit = summary(fit_alt)$coef$mean  # output of function hetglm
  mean_df = cbind.data.frame(tmp_df, 'alt', 'mean', rownames(mean_fit), mean_fit)
  
  var_fit = summary(fit_alt)$coef$scale  # output of function hetglm
  var_df = cbind.data.frame(tmp_df, 'alt', 'var', rownames(var_fit), var_fit)
  
  write_output(null_df, fit_colnames, fit_file)
  write_output(mean_df, fit_colnames, fit_file)
  write_output(var_df, fit_colnames, fit_file)
}

summarize_lrt_binary <- function(pheno_name, model_name, fit_null, fit_alt, lrt_file) {
  # > summary(fit_alt)$lrtest
  # LR         Df    p-value
  # 9.93457039 5.00000000 0.07711054
  # > -2*(summary(fit_alt)$loglik.null - summary(fit_alt)$loglik)
  # [1] 9.93457
  # tmp_df = cbind.data.frame(pheno_name, model_name)       # avoids warnings from cbind.data.frame
  # rownames(tmp_df) = c()
  
  null_row = c(pheno_name, model_name, 'null', attr(logLik(fit_null), 'nobs'), NA, logLik(fit_null), NA, attr(logLik(fit_null), 'df'), NA)         # glm
  alt_row = c(pheno_name, model_name, 'alt', fit_alt$nobs, summary(fit_alt)$loglik.null, summary(fit_alt)$loglik, summary(fit_alt)$lrtest)     # hetglm: c('ll_null', 'll_alt', 'LR_test_stat', 'df', 'pvalue')
  out_df = rbind.data.frame(null_row, alt_row)
  write_output(out_df, lrt_colnames, lrt_file)
}

  
summarize_fits_continuous <- function(pheno_name, model_name, fit_null, fit_alt, fit_file) {
  # write output
  # 1. table of model fits
  # 2. TODO R^2, adj R^2, AIC, BIC, likelihood
  
  tmp_df = cbind.data.frame(pheno_name, model_name)       # avoids warnings from cbind.data.frame
  rownames(tmp_df) = c()

  null_fit = summary(fit_null)$coef      # output of function lm
  null_df = cbind.data.frame(tmp_df, 'null', 'mean', rownames(null_fit), null_fit)
  
  mean_fit = summary(fit_alt)$coefficients  # output of function dglm
  mean_df = cbind.data.frame(tmp_df, 'alt', 'mean', rownames(mean_fit), mean_fit)
  
  var_fit = summary(fit_alt)$dispersion.summary$coefficients  # output of function dglm
  var_df = cbind.data.frame(tmp_df, 'alt', 'var', rownames(var_fit), var_fit)
  
  write_output(null_df, fit_colnames, fit_file)
  write_output(mean_df, fit_colnames, fit_file)
  write_output(var_df, fit_colnames, fit_file)
}

summarize_lrt_continuous <- function(pheno_name, model_name, fit_null, fit_alt, lrt_file) {
  # TODO figure out how to get other fields from dglm: nobs, ll of null model, lrt, etc. see binary code
  # glm
  null_row = c(pheno_name, model_name, 'null', attr(logLik(fit_null), 'nobs'), NA, logLik(fit_null), NA, attr(logLik(fit_null), 'df'), NA)
  
  # dglm
  ll = - fit_alt$m2loglik / 2 
  alt_row = c(pheno_name, model_name, 'alt', NA, NA, ll, NA, NA, NA)  # c('ll_null', 'll_alt', 'LR_test_stat', 'df', 'pvalue')
  
  out_df = rbind.data.frame(null_row, alt_row)
  write_output(out_df, lrt_file)
}


make_models <- function(mean_covar, var_covar, mean_test, var_test, dof, log_file) {
  # Input: lists of mean and variance covariates and mean and variance terms to test
  # Output: additive model with + separating terms
  
  dof_estimated = length(c(mean_test, var_test))
  if(dof != dof_estimated) {
    cat(sprintf('Error: argument dof %d does not match estimated dof %d based on number of parameters\n', dof, dof_estimated), file=log_file, append=T)
  }
  mean_null_model = paste(mean_covar, collapse=' + ')
  mean_alt_model = paste(c(mean_covar, mean_test), collapse=' + ')    
  var_null_model = paste(var_covar, collapse=' + ')
  var_alt_model = paste(c(var_covar, var_test), collapse=' + ')    
  retval = list(mean_null_model, mean_alt_model, var_null_model, var_alt_model)
  names(retval) = c('mean_null_model', 'mean_alt_model', 'var_null_model', 'var_alt_model')
  return(retval)
}


gwas_binary <- function(snp_df, model_name, lrt_name, dof, fit_file, lrt_file, lm_file, snp_name, chr, pos, log_file) {
  retval = make_models(mean_covar, var_covar, mean_test, var_test, dof, log_file)
  mean_null_model = retval$mean_null_model
  mean_alt_model = retval$mean_alt_model
  var_null_model = retval$var_null_model
  var_alt_model = retval$var_alt_model
  
  # run test
  if(model_name == 'lm') {
    mean_formula = as.formula(paste('pheno~', mean_alt_model))
    # specify family and link function to match hetglm defaults    
    fit_null = try(glm(mean_formula, data = snp_df, family = binomial(link = "probit")))
  } else {
    # hetglm package syntax
    null_formula = as.formula(paste('pheno~', mean_null_model, ' | ', var_null_model))  
    alt_formula = as.formula(paste('pheno~', mean_alt_model, ' | ', var_alt_model))   
    fit_null = try(hetglm(null_formula, data = snp_df))
    fit_alt = try(hetglm(alt_formula, data = snp_df))
  }    
  
  # write output
  if(model_name == 'lm') {
    if(! is.error(fit_null)) {
      null_fit = summary(fit_null)$coef      # output of function glm
      tmp_df = cbind.data.frame(snp_name, chr, pos, model_name)       # avoids warnings from cbind.data.frame below
      rownames(tmp_df) = c()
      null_df = cbind.data.frame(tmp_df, 'mean', rownames(null_fit), null_fit)
      write_output(null_df, fit_colnames, lm_file)
    } else {
      cat('Failed lm convergence:', snp_name, model_name, '\n')
    }
  } else {
    if((! is.error(fit_null)) & (! is.error(fit_alt)) & fit_null$converged & fit_alt$converged ) {
      # do LRT
      res = lrtest(fit_null, fit_alt)
      # TODO change to doing LRT as below
      
      tmp_df = cbind.data.frame(snp_name, chr, pos, model_name)       # avoids warnings from cbind.data.frame below
      rownames(tmp_df) = c()
      mean_fit = summary(fit_alt)$coef$mean  # output of function hetglm
      mean_df = cbind.data.frame(tmp_df, 'mean', rownames(mean_fit), mean_fit)
      lrt_df = data.frame(t(c(snp_name, chr, pos, model_name, 'LRT', lrt_name, res[2,2], res[2,3], res[2,4], res[2,5])))
      write_output(lrt_df, lrt_colnames, lrt_file)
      write_output(mean_df, fit_colnames, fit_file)
      if(var_alt_model != "1") {                 # constant dispersion model. hetglm does not fit an intercept so var_fit below will be an empty dataframe
        var_fit = summary(fit_alt)$coef$scale    # specific
        var_df = cbind.data.frame(tmp_df, 'var', rownames(var_fit), var_fit)
        write_output(var_df, fit_colnames, fit_file)
      }
    } else {
      errmsg = cat('Failed adglm convergence:', snp_name, model_name)
      if(! fit_null$converged) {
        errmsg = cat(errmsg, '(null model)')
      }
      if(! fit_alt$converged) {
        errmsg = cat(errmsg, '(alt model)')
      }
      cat(errmsg, '\n')
    }
  }
}

gwas_continuous <- function(snp_df, model_name, lrt_name, dof, fit_file, lrt_file, lm_file, snp_name, chr, pos) {
  retval = make_models(mean_covar, var_covar, mean_test, var_test, dof, log_file)
  mean_null_model = retval$mean_null_model
  mean_alt_model = retval$mean_alt_model
  var_null_model = retval$var_null_model
  var_alt_model = retval$var_alt_model
  
  
  if(model_name == 'lm') {
    mean_formula = as.formula(paste('pheno~', mean_alt_model))
    fit_null = try(lm(mean_formula, data = snp_df))
  } else {
    # dglm package syntax
    mean_formula = as.formula(paste('pheno~', mean_null_model))
    var_formula = as.formula(paste('~', var_null_model))
    fit_null = try(dglm(mean_formula, var_formula, data = snp_df))
    
    mean_formula = as.formula(paste('pheno~', mean_alt_model))
    var_formula = as.formula(paste('~', var_alt_model))
    fit_alt = try(dglm(mean_formula, var_formula, data = snp_df))
  }

  # write output
  if(model_name == 'lm') {
    if(! is.error(fit_null)) {
      null_fit = summary(fit_null)$coef      # output of function dglm
      tmp_df = cbind.data.frame(snp_name, chr, pos, model_name)       # avoids warnings from cbind.data.frame below
      rownames(tmp_df) = c()
      null_df = cbind.data.frame(tmp_df, 'mean', rownames(null_fit), null_fit)
      write_output(null_df, fit_colnames, lm_file)
    } else {
      cat('Failed lm convergence:', snp_name, model_name, '\n')
    }
  } else {
    if((! is.error(fit_null)) & (! is.error(fit_alt))) {
      # do LRT
      chiStat = fit_null$m2loglik - fit_alt$m2loglik    # m2loglik is minus twice the log likelihood = 2 * NLL
      pval = exp(pchisq(chiStat, df=dof, lower.tail=F, log.p=T))      # alternate for small stat
      res = c(-2*fit_alt$m2loglik, dof, chiStat, pval)     # NLL_alt df lambda_stat pval  # TODO WRONG FOR SOME TESTS
      
      tmp_df = cbind.data.frame(snp_name, chr, pos, model_name)       # avoids warnings from cbind.data.frame below
      rownames(tmp_df) = c()
      mean_fit = summary(fit_alt)$coefficients
      mean_df = cbind.data.frame(tmp_df, 'mean', rownames(mean_fit), mean_fit)
      var_fit = summary(fit_alt)$dispersion.summary$coefficients
      var_df = cbind.data.frame(tmp_df, 'var', rownames(var_fit), var_fit)
      lrt_df = data.frame(t(c(snp_name, chr, pos, model_name, 'LRT', sub('adglm_', '', model_name), res)))
      
      write_output(lrt_df, lrt_colnames, lrt_file)
      write_output(mean_df, fit_colnames, fit_file)
      write_output(var_df, fit_colnames, fit_file)
    } else {
      cat('Failed adglm convergence:', snp_name, model_name, '\n')
    }
  }
}

