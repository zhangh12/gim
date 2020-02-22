
## formula: full model for internal data
## model: working models for external data
## data: internal data
## ref: reference panel


#####################################################################################
## list of formula supported in this version
## (1) y ~ x - 1 # real example in paper (bet1, bet2 given). Feature: -1
## (2) log(y) ~ x # example in help doc. Feature: log in outcome
## (3) y ~ I(x==0) -1 # real example in paper (bet1 given)
#####################################################################################
formula.parse <- function(formula, family, data, model, ref = NULL){
  
  form0 <- as.formula(formula)
  if(family == 'gaussian'){
    fit0 <- lm(form0, data)
  }else{
    fit0 <- glm(form0, family = 'binomial', data = data)
  }
  
  if(1){
    # this does not work for log(y) ~ x. It will return 'y' rather than 'log(y)'
    # use the following instead
    ori.outcome <- all.vars(form0)[1]
  }
  
  # if ref is not specified, it assumes that internal data could be used as reference
  # which is used to construct constrain equations
  type <- 'others'
  if(is.null(ref)){
    ref <- data
  }else{
    if(family == 'case-control'){
      type <- 'cc-ref' # case-control data with specific external reference
    }
  }
  
  if(ori.outcome %in% colnames(ref)){
    # do not discard a line in reference if it only misses its outcome
    if(any(is.na(ref[, ori.outcome]))){
      id <- sample(1:nrow(data), nrow(ref), TRUE)
      # it doesn't matter how to impute NA in outcome in ref
      # because gim will not use its value
      # I did this because I do not want model.matrix to delete incomplete lines in ref
      # outcome is usually missing in ref, so I have to impute it.
      ref[, ori.outcome] <- data[id, ori.outcome]
    }
  }
  
  # I am going to parsing formula and extract/create necessary variables (transformed according to formula)
  # so need to bind data and ref together
  # after those variables are ready, data and ref will be split again
  # int.id lets me know which rows are internal data, and which rows are ref
  int.id <- 1:nrow(data)
  data0 <- data
  
  # ref may miss some columns in data (e.g. outcome)
  # I will add those missed columns to ref
  # and rbind ref and data for parsing
  # values in added columns are copied from data, which do not have to make sense
  # those values could not be NA, otherwise model.frame or model.matrix will delete uncomplete rows
  # I will remove those added columns before returning (see the end of this function)
  # so their imputed values do not matter
  var1 <- colnames(data)
  var2 <- colnames(ref)
  var3 <- intersect(var1, var2)
  if(length(var3) == 0){
    stop('data and ref do not share any variable')
  }
  
  ref <- ref[, var3, drop = FALSE]
  var4 <- setdiff(var1, var2)
  if(length(var4) > 0){
    id <- sample(1:nrow(data), nrow(ref), TRUE)
    tmp <- data[id, var4, drop = FALSE]
    rownames(tmp) <- NULL
    # randomly impute some value so that model.matrix below can work
    # I use values in data to impute
    # it is neccessary as it does not introduce new factor level
    # if imputed to be 0, then 0 might be new level of a factor variable which affect parsing
    # I will set it to NA before returning
    ref[, var4] <- tmp
    rm(tmp)
  }
  ref <- ref[, colnames(data)]
  data <- rbind(data, ref)
  data$'(Intercept)' <- 1
  
  
  # extract name of transformed outcome 
  # if log(y) ~ ..., then outcome should be 'log(y)', not 'y'
  # we aim to creat a data frame with a column named log(y)
  if(1){
    mf0 <- model.frame(form0, data = data)
    outcome <- colnames(mf0)[1]
  }
  
  # create design matrix
  # expand factor to be dummy variables
  # expand interaction
  # expand transformation, e.g. I(), log, etc
  # add intercept term (all 1, named '(Intercept)')
  mat0 <- model.matrix(form0, data = data)
  
  miss.var <- NULL
  nform <- length(model)
  # parse each of working models
  # find requested covariates and put it to mat0
  for(i in 1:nform){
    f <- as.formula(model[[i]]$form)
    mat <- model.matrix(f, data = data) # desing matrix for working model
    var.in.form <- colnames(mat)
    add.var <- setdiff(var.in.form, model[[i]]$info$var)
    if(length(add.var) > 0){
      tmp1 <- mat[, add.var, drop = FALSE]
      tmp2 <- mat[, model[[i]]$info$var, drop = FALSE]
      rm.col <- NULL
      for(j in 1:ncol(tmp1)){
        for(k in 1:ncol(tmp2)){
          if(all(tmp1[, j] == tmp2[, k])){
            rm.col <- c(rm.col, j)
          }
        }
      }
      if(!is.null(rm.col)){
        rm.col <- colnames(tmp1)[rm.col]
        mat <- mat[, setdiff(var.in.form, rm.col), drop = FALSE]
      }
    }
    
    new.var <- setdiff(colnames(mat), colnames(mat0))
    if(length(new.var) > 0){
      mat0 <- cbind(mat0, mat[, new.var, drop = FALSE])
    }
    
    miss.var <- c(miss.var, setdiff(model[[i]]$info$var, colnames(mat)))
    covar <- setdiff(colnames(mat), model[[i]]$info$var)
    if(length(covar) == 0){
      covar <- NULL
    }else{
      covar <- as.character(covar)
    }
    
    if(family == 'gaussian'){
      fit <- lm(f, data = data0)
    }else{
      fit <- glm(f, family = 'binomial', data = data0)
    }
    # covar is nuisance variables specified in working model but no summary data available
    model[[i]] <- list(f, covar, model[[i]]$info, fit)
    rm(fit)
  }
  
  if(length(miss.var) > 0){
    msg <- paste0('The following variable(s) specified in model$info does not match those in model$formula: \n"', paste0(miss.var, collapse = '", "'), '"\nTry to find proper form(s) for these missing variable(s) from the following names and modify your model input: \n"', paste0(colnames(mat0), collapse = '", "'), '"')
    stop(msg)
  }
  
  mat0 <- cbind(mf0[, outcome, drop = FALSE], mat0)
  
  ori.var <- setdiff(colnames(data), colnames(mat0))
  if(length(ori.var) > 0){
    mat0 <- cbind(mat0, data[, ori.var, drop = FALSE])
  }
  
  mat0 <- as.data.frame(mat0)
  
  data <- mat0[int.id, , drop = FALSE]
  ref <- mat0[-int.id, , drop = FALSE]
  
  if(length(var4) > 0){
    ref[, var4] <- NA
  }
  
  list(model = model, data = data, ref = ref, fit0 = fit0, 
       outcome = outcome, formula = form0, type = type)
  
}
