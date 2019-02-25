
# lam: Lagrange multiplier
# the: parameters in full model (internal data)
# alp: nuisance parameter in working model (external data)
# bet: parameters in working model (with summary data)
# para := c(lam, the, alp, bet)
# bet0 := summary statistics to be used in quadratic form

## formula: full model (internal data)
## data: internal data
## model: working model and summary data
## nsample: overlapped sample size

init.lm <- function(formula, data, model, nsample){
  
  #message('Initializing integration analysis...')
  
  if(is.null(nsample)){
    msg <- 'nsample could not be NULL'
    stop(msg)
  }
  
  nsample <- as.matrix(nsample)
  
  fit0 <- glm(formula, data = data, family = 'gaussian')
  the <- c(mean(fit0$residuals^2), coef(fit0))
  names(the)[1] <- 'sigma'
  
  n <- nrow(data)
  nmodel <- length(model)
  
  alp <- NULL
  bet <- NULL
  bet0 <- NULL
  
  map <- list()
  map$alp <- list(0:0)
  map$bet <- list(0:0)
  
  for(i in 1:nmodel){
    
    form <- model[[i]][[1]]
    fit <- glm(form, data = data, family = 'gaussian')
    alp.var <- as.character(model[[i]][[2]]) # nuisance variables in working model
    bet.var <- as.character(model[[i]][[3]]$var)
    N <- diag(nsample)[i] # sample size for working model
    # we assume estimate of error parameter always not provided
    # so alp0 has at least one entry
    alp0 <- c(mean(fit$residuals^2), coef(fit)[alp.var])
    names(alp0)[1] <- paste0('tau', i) # always be the variance parameter in the first entry
    meta.bet <- (n * coef(fit)[bet.var] + N * model[[i]][[3]]$bet) / (n + N)
    #meta.bet <- model[[i]][[3]]$bet
    names(meta.bet) <- model[[i]][[3]]$var
    alp <- c(alp, alp0)
    bet <- c(bet, meta.bet)
    bet0 <- c(bet0, model[[i]][[3]]$bet)
    
    map$alp[[i + 1]] <- max(map$alp[[i]]) + 1:length(alp0)
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
    map$alp[[i]] <- map$alp[[i]] + nlam + nthe
    map$bet[[i]] <- map$bet[[i]] + nlam + nthe + nalp
    all.alp <- c(all.alp, map$alp[[i]])
    all.bet <- c(all.bet, map$bet[[i]])
  }
  
  map$all.alp <- all.alp
  map$all.bet <- all.bet
  
  list(para = para, map = map, bet0 = bet0, sample.info = nsample)
  
}

