

formula.parse <- function(formula, model, data, ref = NULL){
  
  if(is.null(ref)){
    ref <- data
  }
  int.id <- 1:nrow(data)
  ref <- ref[, colnames(data)]
  data <- rbind(data, ref)
  
  form0 <- as.formula(formula)
  #outcome <- all.vars(form0)[1] # this does not work for log(y) ~ x. It will return 'y' rather than 'log(y)'
  mf0 <- model.frame(form0, data = data)
  outcome <- colnames(mf0)[1]
  mat0 <- model.matrix(form0, data = data)
  # need the following otherwise a formula like y~. will result in an error
  #form0 <- paste0(outcome, ' ~ ', paste0(setdiff(colnames(mat0), '(Intercept)'), collapse = ' + '))
  
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
  
  data <- mat0[int.id, , drop = FALSE]
  ref <- mat0[-int.id, , drop = FALSE]
  
  list(model = model, data = data, ref = ref, outcome = outcome, formula = form0)
  
}
