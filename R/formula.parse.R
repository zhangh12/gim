

formula.parse <- function(formula, model, data){
  
  form0 <- as.formula(formula)
  #outcome <- all.vars(form0)[1] # this does not work for log(y) ~ x. It will return 'y' rather than 'log(y)'
  mf0 <- model.frame(form0, data = data)
  outcome <- colnames(mf0)[1]
  mat0 <- model.matrix(form0, data = data)
  
  miss.var <- NULL
  nform <- length(model)
  for(i in 1:nform){
    f <- as.formula(model[[i]]$form)
    mat <- model.matrix(f, data = data)
    new.var <- setdiff(colnames(mat), colnames(mat0))
    if(length(new.var) > 0){
      mat0 <- cbind(mat0, mat[, new.var, drop = FALSE])
    }
    
    miss.var <- c(miss.var, setdiff(model[[i]]$info$var, colnames(mat)))
    covar <- setdiff(colnames(mat), model[[i]]$info$var)
    if(length(covar) == 0){
      covar <- NULL
    }
    
    model[[i]] <- list(f, covar, model[[i]]$info)
  }
  
  if(length(miss.var) > 0){
    msg <- paste0('The following variable(s) specified in model$info does not match those in model$formula: \n"', paste0(miss.var, collapse = '", "'), '"\nTry to find proper form(s) for these missing variable(s) from the following names and modify your model input: \n"', paste0(colnames(mat0), collapse = '", "'), '"')
    stop(msg)
  }
  
  mat0 <- cbind(mf0[, outcome], mat0)
  colnames(mat0)[1] <- outcome
  
  ori.var <- setdiff(colnames(data), colnames(mat0))
  if(length(ori.var) > 0){
    mat0 <- cbind(mat0, data[, ori.var, drop = FALSE])
  }
  
  mat0 <- as.data.frame(mat0)
  
  list(model = model, data = mat0, outcome = outcome)
  
}
