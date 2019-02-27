
clean.data <- function(dat){
  
  dat <- as.data.frame(dat)
  n.miss <- sum(!complete.cases(dat))
  dat <- dat[complete.cases(dat), ]
  if(n.miss > 0){
    msg <- paste0(n.miss, ' samples with missing entries are excluded from analysis. If you think this number is incorrect, remove useless variables from \'data\' and try again')
    message(msg)
  }
  
  if(nrow(dat) == 0){
    msg <- paste0('No data available for analysis after data cleaning')
    stop(msg)
  }
  
  var.class <- sapply(dat, class)
  id.num <- which(var.class %in% c('numeric', 'integer'))
  id.fac <- which(var.class %in% 'factor')
  for(i in id.fac){
    dat[, i] <- as.character(dat[, i])
  }
  id.cha <- which(var.class %in% c('character'))
  
  rm.var <- NULL
  for(i in id.num){
    if(sd(dat[, i]) == 0){
      rm.var <- c(rm.var, i)
    }
  }
  
  for(i in id.cha){
    if(length(unique(dat[, i])) == 1){
      rm.var <- c(rm.var, i)
    }
  }
  
  if(!is.null(rm.var)){
    rm.var <- colnames(dat)[rm.var]
    msg <- paste('Variable(s)', paste(rm.var, collapse = ', '), 'are constant. Remove them from \'data\', modify \'formula\' and \'model\' accordingly, then try again')
    stop(msg)
  }
  
  dat
  
}
