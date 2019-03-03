
# tau is not in alp
gfunction.alp.bet.lm <- function(para, map, ref, lam){
  
  nmodel <- length(map$bet)
  
  g.alp.bet <- list()
  
  foo <- function(j, l){
    paste0(j,'-',l)
  }
  
  offset <- max(map$the)
  for(i in 1:nmodel){
    id.tau <- map$alp[[i]][1]
    id.a <- alp.index.lm(map, i)
    if(is.null(id.a)){
      next
    }
    
    id.b <- map$bet[[i]]
    
    for(j in id.b){
      fxj <- ref[, names(para)[j]]
      for(l in id.a){
        fxl <- ref[, names(para)[l]]
        g.alp.bet[[foo(j,l)]] <- 2 * fxj * fxl * lam[id.tau - offset]
      }
    }
  }
  
  if(length(g.alp.bet) == 0){
    return(NULL)
  }
  
  g.alp.bet
  
}

