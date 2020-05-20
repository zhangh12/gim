
# lam: Lagrange multiplier
# the: parameters in full model (internal data)
# alp: nuisance parameter in working model (external data)
# bet: parameters in working model (with summary data)
# para := c(lam, the, alp, bet)
# bet0 := summary statistics to be used in quadratic form


init.lo <- function(fit0, data, model, nsample){
  
  #message('Initializing integration analysis...')
  
  if(is.null(nsample)){
    msg <- 'nsample could not be NULL'
    stop(msg)
  }
  
  nsample <- as.matrix(nsample)
  
  #fit0 <- glm(formula, data = data, family = 'binomial')
  the <- coef(fit0)
  
  n <- nrow(data)
  nmodel <- length(model)
  
  alp <- NULL
  bet <- NULL
  bet0 <- NULL
  
  map <- list()
  map$alp <- list(0:0)
  map$bet <- list(0:0)
  
  no.alp <- NULL
  for(i in 1:nmodel){
    
    form <- model[[i]][[1]]
    #fit <- glm(form, data = data, family = 'binomial')
    fit <- model[[i]][[4]]
    alp.var <- as.character(model[[i]][[2]])
    bet.var <- as.character(model[[i]][[3]]$var) # critical for meta-analysis of bet below
    N <- diag(nsample)[i]
    
    if(length(alp.var) > 0){
      alp0 <- coef(fit)[alp.var]
    }else{
      alp0 <- NULL
    }
    
    meta.bet <- (n * coef(fit)[bet.var] + N * model[[i]][[3]]$bet) / (n + N)
    #meta.bet <- model[[i]][[3]]$bet
    names(meta.bet) <- model[[i]][[3]]$var
    alp <- c(alp, alp0)
    bet <- c(bet, meta.bet)
    bet0 <- c(bet0, model[[i]][[3]]$bet)
    
    if(!is.null(alp0)){
      map$alp[[i + 1]] <- max(map$alp[[i]]) + 1:length(alp0)
    }else{
      map$alp[[i + 1]] <- map$alp[[i]]
      no.alp <- c(no.alp, i)
    }
    
    map$bet[[i + 1]] <- max(map$bet[[i]]) + 1:length(meta.bet)
    
    rm(form, fit, alp.var, bet.var, N, alp0, meta.bet)
  }
  
  map$alp[[1]] <- NULL
  map$bet[[1]] <- NULL
  
  nlam <- length(alp) + length(bet)
  lam <- rep(0.01, nlam)
  names(lam) <- paste0('lam', 1:nlam)
  
  para <- c(lam, the, alp, bet)
  
  nthe <- length(the)
  nalp <- length(alp)
  
  map$lam <- 1:nlam
  map$the <- (nlam + 1) : (nlam + nthe)
  all.alp <- NULL
  all.bet <- NULL
  for(i in 1:nmodel){
    if(i %in% no.alp){
      map$alp[[i]] <- NA
    }else{
      map$alp[[i]] <- map$alp[[i]] + nlam + nthe
      all.alp <- c(all.alp, map$alp[[i]])
    }
    map$bet[[i]] <- map$bet[[i]] + nlam + nthe + nalp
    all.bet <- c(all.bet, map$bet[[i]])
  }
  
  if(!is.null(all.alp)){
    map$all.alp <- sort(unique(all.alp))
  }
  
  map$all.bet <- sort(unique(all.bet))
  
  list(para = para, map = map, bet0 = bet0, sample.info = nsample)
  
}

