

## File processing ----------------------------------------------------------------

qnorm_norm_with_NA <- function(x, debug=F) {
  # x: numerical vector that can have NAs
  # typical usage: x is a phenotype or covariate
  
  num_unique_values = length(unique(x[!is.na(x)]))
  if(num_unique_values < 5) {
    stop('Quitting: refusing to quantile normalize a vector with less than 5 unique values\n')
  } 
  if(debug) cat('Before qnorm summary and sd: ', summary(x), sd(x, na.rm=T), '\n')
  num_na_before = sum(is.na(x))
  
  non_na_idx = which(!is.na(x))
  x[non_na_idx] = quant_norm(x[non_na_idx])
  
  if(debug) cat('After qnorm summary and sd:  ', summary(x), sd(x, na.rm=T), '\n')
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

categorize_pheno <- function(pheno, log_file, do_print=T) {
  # returns number of unique phenotype values  
  # is_binary: determines whether logistic or continuous test is run
  
  num_unique_values = length(unique(pheno[!is.na(pheno)]))
  if(do_print) {
    if(num_unique_values == 2) {
      cat('Binary phenotype: has 2 values\n', file=log_file, append=T)
    } else if(num_unique_values <= 10) {
      cat(sprintf('Categorical phenotype or binary phenotype with unidentified null characters: has %d values\n', num_unique_values), file=log_file, append=T)
    } else {
      cat(sprintf('Continuous phenotype or categorical phenotype with more than ten values: has %d values\n', num_unique_values), file=log_file, append=T)
    }
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
  
  # TODO the same variable is in dataframe twice if pheno
  covar = cbind.data.frame(all_df, pheno)     # add column named pheno
  cat('Covariate dataframe has', nrow(covar), 'individuals.\n', file=log_file, append=T)
  covar = covar[,c('pheno', model_terms)]
  num_incomplete = sum(!complete.cases(covar))
  if(num_incomplete > 0) {
    covar = covar[complete.cases(covar),]
    cat('Removed', num_incomplete, 'individuals with incomplete phenotype and/or covariate measurements.\n', file=log_file, append=T)
    cat('Final dataframe with no NAs has', nrow(covar), 'individuals.\n', file=log_file, append=T)
  }
  return(covar)   
}


## Model fitting and testing ------------------------------------------------------

make_models <- function(mean_covar, var_covar, mean_test, var_test, dof, log_file) {
  # Input: lists of mean and variance covariates and mean and variance terms to test
  # Output: additive model with + separating terms
  
  mean_null_model = paste(mean_covar, collapse=' + ')
  mean_alt_model = paste(c(mean_covar, mean_test), collapse=' + ')    
  var_null_model = paste(var_covar, collapse=' + ')
  var_alt_model = paste(c(var_covar, var_test), collapse=' + ')    
  retval = list(mean_null_model, mean_alt_model, var_null_model, var_alt_model)
  names(retval) = c('mean_null_model', 'mean_alt_model', 'var_null_model', 'var_alt_model')
  return(retval)
}

run_models_discrete <- function(pheno_name, model_name, mean_null_model, mean_alt_model, var_null_model, var_alt_model, dof, covar, fit_file, lrt_file, log_file, terse=F) {

  mean_formula = as.formula(paste('pheno~', mean_alt_model))
  # specify family and link function to match hetglm defaults    
  fit_lm = try(glm(mean_formula, data = covar, family = binomial(link = "probit")))
  
  null_formula = as.formula(paste('pheno~', mean_null_model, ' | ', var_null_model))  
  alt_formula = as.formula(paste('pheno~', mean_alt_model, ' | ', var_alt_model))   
  fit_null = try(hetglm(null_formula, data = covar))
  fit_alt = try(hetglm(alt_formula, data = covar))
  
  if(! fit_alt$converged) {
    cat(sprintf('Warning: hetglm did not converge'), file=log_file, append=T)
                # for phenotype %s on model %s with mean parameters %s and variance parameters %s.\n', pheno_name, model_name, toString(mean_terms), toString(var_terms)), file=log_file, append=T)
  } else {
      cat(sprintf('Success: hetglm converged'), file=log_file, append=T)
                  # for phenotype %s on model %s with mean parameters %s and variance parameters %s.\n', pheno_name, model_name, toString(mean_terms), toString(var_terms)), file=log_file, append=T)
  }
  
  if(terse) {
    tmp_df = cbind.data.frame(pheno_name, model_name, pheno_transform, theta_transform)       # avoids warnings from cbind.data.frame
    rownames(tmp_df) = c()
    null_fit = summary(fit_lm)$coef      # output of function glm
    null_df = cbind.data.frame(tmp_df, 'null', 'mean', rownames(null_fit), null_fit)
    mean_fit = summary(fit_alt)$coef$mean  # output of function hetglm
    mean_df = cbind.data.frame(tmp_df, 'alt', 'mean', rownames(mean_fit), mean_fit)
    var_fit = summary(fit_alt)$coef$scale  # output of function hetglm
    var_df = cbind.data.frame(tmp_df, 'alt', 'var', rownames(var_fit), var_fit)
    colnames(null_df) = fit_colnames;     colnames(mean_df) = fit_colnames;     colnames(var_df) = fit_colnames
    fit_df = rbind.data.frame(mean_df, var_df)    
    
    ll_null = summary(fit_alt)$loglik.null    
    ll_alt = summary(fit_alt)$loglik
    test_stat = - 2 * (ll_null - ll_alt)
    lrt_pval = exp(pchisq(test_stat, df=dof, lower.tail=F, log.p=T))
    lrt_df = cbind.data.frame(tmp_df, ll_null, ll_alt, test_stat, dof, lrt_pval)
    retval = list(lrt_df, fit_df)
    names(retval) = c('lrt_df', 'fit_df')
  } else {
    summarize_fits_discrete(pheno_name, model_name, fit_null, fit_alt, fit_file)
    summarize_lrt_discrete(pheno_name, model_name, fit_null, fit_alt, lrt_file)
    retval = NA
  }
  invisible(retval)
}

run_models_continuous <- function(pheno_name, model_name, mean_null_model, mean_alt_model, var_null_model, var_alt_model, dof, covar, fit_file, lrt_file, log_file, terse=F) {
  if(any(is.na(covar))) {
    cat('Error: covariate matrix has NAs and package dglm will fail. Run function subset_covar_df\n', file=log_file, append=T)
  }
  
  mean_formula = as.formula(paste('pheno~', mean_alt_model))
  fit_lm = try(lm(mean_formula, data = covar))   # TODO not used below

  mean_formula = as.formula(paste('pheno~', mean_null_model))
  var_formula = as.formula(paste('~', var_null_model))
  fit_null = try(dglm(mean_formula, var_formula, data = covar))
  
  mean_formula = as.formula(paste('pheno~', mean_alt_model))
  var_formula = as.formula(paste('~', var_alt_model))
  fit_alt = try(dglm(mean_formula, var_formula, data = covar))
  
  num_iter_ran = fit_alt$iter
  max_iter = dglm.control()$maxit
  if(num_iter_ran >= max_iter) {
    cat(sprintf('Warning: dglm ran for %d iterations and reached maximum number of iterations %d for phenotype %s on model %s\n', num_iter_ran, max_iter, pheno_name, model_name), file=log_file, append=T)
                # with mean parameters %s and variance parameters %s.\n', num_iter_ran, max_iter, pheno_name, model_name, toString(mean_terms), toString(var_terms)), file=log_file, append=T)
  } else {
    cat(sprintf('Ok: dglm ran for %d iterations which is less than maximum number of iterations %d for phenotype %s on model %s\n', num_iter_ran, max_iter, pheno_name, model_name), file=log_file, append=T)
                # with mean parameters %s and variance parameters %s.\n', num_iter_ran, max_iter, pheno_name, model_name, toString(mean_terms), toString(var_terms)), file=log_file, append=T)
  }

  if(terse) {
    tmp_df = cbind.data.frame(pheno_name, model_name, pheno_transform, theta_transform)       # avoids warnings from cbind.data.frame
    rownames(tmp_df) = c()
    null_fit = summary(fit_lm)$coef      # output of function lm
    null_df = cbind.data.frame(tmp_df, 'null', 'mean', rownames(null_fit), null_fit)
    mean_fit = summary(fit_alt)$coefficients  # output of function dglm
    mean_df = cbind.data.frame(tmp_df, 'alt', 'mean', rownames(mean_fit), mean_fit)
    var_fit = summary(fit_alt)$dispersion.summary$coefficients  # output of function dglm
    var_df = cbind.data.frame(tmp_df, 'alt', 'var', rownames(var_fit), var_fit)
    colnames(null_df) = fit_colnames;     colnames(mean_df) = fit_colnames;     colnames(var_df) = fit_colnames
    fit_df = rbind.data.frame(mean_df, var_df)
    
    ll_null = logLik(fit_null)
    ll_alt = - fit_alt$m2loglik / 2
    test_stat = - 2 * (ll_null - ll_alt)
    lrt_pval = exp(pchisq(test_stat, df=dof, lower.tail=F, log.p=T))
    lrt_df = cbind.data.frame(tmp_df, ll_null, ll_alt, test_stat, dof, lrt_pval)
    retval = list(lrt_df, fit_df)
    names(retval) = c('lrt_df', 'fit_df')
  } else {
    summarize_fits_continuous(pheno_name, model_name, fit_null, fit_alt, fit_file)
    summarize_lrt_continuous(pheno_name, model_name, fit_null, fit_alt, lrt_file, dof)
    retval = NA
  }
  invisible(retval)
}



## Writing output -----------------------------------------------------------------

write_output <- function(fit_df, outfile) {
  colnames(fit_df) = fit_colnames
  write.table(fit_df, outfile, row.names=F, col.names=F, quote=F, append=T)
}

summarize_fits_discrete <- function(pheno_name, model_name, fit_null, fit_alt, fit_file) {
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
  
  write_output(null_df, fit_file)
  write_output(mean_df, fit_file)
  write_output(var_df, fit_file)
}

summarize_lrt_discrete <- function(pheno_name, model_name, fit_null, fit_alt, lrt_file) {
  # > summary(fit_alt)$lrtest
  # LR         Df    p-value
  # 9.93457039 5.00000000 0.07711054
  # > -2*(summary(fit_alt)$loglik.null - summary(fit_alt)$loglik)
  # [1] 9.93457

  null_row = c(pheno_name, model_name, 'null', attr(logLik(fit_null), 'nobs'), NA, logLik(fit_null), NA, attr(logLik(fit_null), 'df'), NA)         # glm
  alt_row = c(pheno_name, model_name, 'alt', fit_alt$nobs, summary(fit_alt)$loglik.null, summary(fit_alt)$loglik, summary(fit_alt)$lrtest)     # hetglm: c('ll_null', 'll_alt', 'LR_test_stat', 'df', 'pvalue')
  out_df = rbind.data.frame(null_row, alt_row)
  write_output(out_df, lrt_file)
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
  
  write_output(null_df, fit_file)
  write_output(mean_df, fit_file)
  write_output(var_df, fit_file)
}

summarize_lrt_continuous <- function(pheno_name, model_name, fit_null, fit_alt, lrt_file, dof) {
  # TODO figure out how to get other fields from dglm: nobs, ll of null model, lrt, etc. see discrete code
  # glm
  # m2loglikminus is twice the log-likelihood or adjusted profile likelihood of the fitted model

  null_row = c(pheno_name, model_name, 'null', attr(logLik(fit_null), 'nobs'), NA, logLik(fit_null), NA, attr(logLik(fit_null), 'df'), NA)
  
  # dglm
  ll_null = logLik(fit_null)
  ll_alt = - fit_alt$m2loglik / 2
  test_stat = - 2 * (ll_null - ll_alt)
  lrt_pval = exp(pchisq(test_stat, df=dof, lower.tail=F, log.p=T))
  alt_row = c(pheno_name, model_name, 'alt', NA, NA, ll_alt, test_stat, dof, lrt_pval)  # c('ll_null', 'll_alt', 'LR_test_stat', 'df', 'pvalue')
  out_df = rbind.data.frame(null_row, alt_row)
  colnames(out_df) = lrt_colnames
  write.table(out_df, lrt_file, row.names=F, col.names=T, quote=F)
}

