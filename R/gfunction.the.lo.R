

gfunction.the.lo <- function(para, map, data){
  
  nmodel <- length(map$bet)
  
  the <- para[map$the]
  fx <- as.matrix(data[, names(the), drop = FALSE])
  
  g.the <- list()
  
  nthe <- length(the)
  n <- nrow(data)
  
  nlam <- max(map$lam)
  offset <- max(map$the)
  
  for(j in 1:nthe){
    gt <- matrix(0, nrow = n, ncol = nlam)
    fx0 <- fx[, names(the)[j]]
    for(i in 1:nmodel){
      
      id.a <- alp.index.lo(map, i)
      alp.exist <- !is.null(id.a)
      if(alp.exist){
        alp <- para[id.a]
      }else{
        alp <- NULL
      }
      
      id.b <- map$bet[[i]]
      bet <- para[id.b]
      gam <- c(alp, bet)
      
      rx <- as.matrix(data[, names(gam), drop = FALSE])
      
      tmp <- as.vector(exp(fx %*% the))
      yh <- tmp/(1+tmp)
      
      id <- c(id.a, id.b)
      gt[, id - offset] <- rx[, names(para)[id], drop = FALSE] * (fx0 * yh * (1-yh))
      
      rm(id.a, id.b, alp.exist)
      
    }
    g.the[[j]] <- gt
    rm(gt)
  }
  
  g.the
  
}


