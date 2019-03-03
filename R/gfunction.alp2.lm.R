
# tau is not in alp
gfunction.alp2.lm <- function(para, map, ref, lam){
  
  nmodel <- length(map$bet)
  
  g.alp2 <- list()
  
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
    
    for(j in id.a){
      fxj <- ref[, names(para)[j]]
      for(l in id.a){
        if(l < j){
          next
        }
        fxl <- ref[, names(para)[l]]
        g.alp2[[foo(j,l)]] <- 2 * fxj * fxl * lam[id.tau - offset]
      }
    }
  }
  
  g.alp2
  
}

