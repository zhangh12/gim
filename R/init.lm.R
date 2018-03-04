

init.lm <- function(formula, data, model, nsample){
  
  message('Initializing integration analysis...')
  
  fit0 <- glm(formula, data = data, family = 'gaussian')
  the <- c(mean(fit0$residuals^2), coef(fit0))
  names(the)[1] <- 'sigma'
  
  n <- nrow(data)
  nmodel <- length(model)
  
  alp <- NULL
  bet <- NULL
  bet0 <- NULL
  id.alp <- data.frame(start = 0, end = 0)
  id.bet <- data.frame(start = 0, end = 0)
  
  for(i in 1:nmodel){
    
    form <- model[[i]][[1]]
    fit <- glm(form, data = data, family = 'gaussian')
    alp.var <- as.character(model[[i]][[2]])
    bet.var <- as.character(model[[i]][[3]]$var)
    N <- diag(nsample)[i]
    alp0 <- c(mean(fit$residuals^2), coef(fit)[alp.var])
    names(alp0)[1] <- paste0('tau', i) # always be the variance parameter in the first entry
    meta.bet <- (n * coef(fit)[bet.var] + N * model[[i]][[3]]$bet) / (n + N)
    #meta.bet <- model[[i]][[3]]$bet
    names(meta.bet) <- model[[i]][[3]]$var
    alp <- c(alp, alp0)
    bet <- c(bet, meta.bet)
    bet0 <- c(bet0, model[[i]][[3]]$bet)
    
    tmp <- data.frame(start = id.alp$end[i] + 1, end = id.alp$end[i] + length(alp0))
    id.alp <- rbind(id.alp, tmp)
    tmp <- data.frame(start = id.bet$end[i] + 1, end = id.bet$end[i] + length(meta.bet))
    id.bet <- rbind(id.bet, tmp)
    rm(form, fit, alp.var, bet.var, N, alp0, meta.bet, tmp)
  }
  
  id.alp <- id.alp[-1, ]
  id.bet <- id.bet[-1, ]
  
  nlam <- length(alp) + length(bet)
  lam <- rep(0.01, nlam)
  names(lam) <- paste0('lam', 1:nlam)
  
  para <- c(lam, the, alp, bet)
  
  nthe <- length(the)
  nalp <- length(alp)
  
  id.lam <- data.frame(start = 1, end = nlam)
  id.the <- data.frame(start = nlam + 1, end = nlam + nthe)
  id.alp$start <- id.alp$start + nlam + nthe
  id.alp$end <- id.alp$end + nlam + nthe
  id.bet$start <- id.bet$start + nlam + nthe + nalp
  id.bet$end <- id.bet$end + nlam + nthe + nalp
  
  para.id <- list(id.lam = id.lam, id.the = id.the, id.alp = id.alp, id.bet = id.bet)
  
  list(para = para, para.id = para.id, bet0 = bet0)
  
}

