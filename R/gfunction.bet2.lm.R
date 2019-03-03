
gfunction.bet2.lm <- function(para, map, ref, lam){
  
  nmodel <- length(map$bet)
  
  g.bet2 <- list()
  
  foo <- function(j, l){
    paste0(j,'-',l)
  }
  
  offset <- max(map$the)
  for(i in 1:nmodel){
    id.tau <- map$alp[[i]][1]
    id.b <- map$bet[[i]]
    for(j in id.b){
      fxj <- ref[, names(para)[j]]
      for(l in id.b){
        if(l < j){
          next
        }
        fxl <- ref[, names(para)[l]]
        g.bet2[[foo(j,l)]] <- 2 * fxj * fxl * lam[id.tau - offset]
      }
    }
  }
  
  g.bet2
  
}

