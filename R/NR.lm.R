
# Newton-Raphson algorithm

NR.lm <- function(para, para.id, int, V, bet0, outcome = 'y'){
  
  inv.V <- solve(V)
  
  message('Running Newton-Raphson algorithm...')
  
  i <- 0
  while(i<100){
    s0 <- score.lm(para, para.id, int, inv.V, bet0, outcome)
    #s1 <- grad(obj.lm, para, para.id = para.id, int = int, inv.V = inv.V, bet0 = bet0, outcome = outcome)
    if(all(abs(s0) < 1e-6)){
      break
    }
    
    t0 <- try(inv.h0 <- solve(hess.lm(para, para.id, int, inv.V, bet0, outcome)), silent = TRUE)
    if('try-error' %in% class(t0)){
      return(list(est = rep(NA, length(para)), sc = rep(NA, length(para)), conv =0))
    }
    d0 <- as.vector(inv.h0 %*% s0)
    if(max(abs(d0)) > 1){
      #d0 <- d0/max(abs(d0))
    }
    para <- para - d0
    
    i <- i + 1
    print(s0)
  }
  
  #svd(inv.h0)$d
  
  print(s0)
  
  list(est = para, sc = s0, conv = ifelse(all(abs(s0) < 1e-6), 1, 0))
  
}
