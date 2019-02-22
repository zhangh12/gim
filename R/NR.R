
# Newton-Raphson algorithm
NR <- function(para, map, family, data, ref, V, bet0, outcome){
  
  inv.V <- solve(V)
  np <- length(para)
  para.null <- rep(NA, np)
  
  i <- 0
  while(i<100){
    if(family == 'gaussian'){
      s0 <- score.lm(para, map, data, ref, inv.V, bet0, outcome)
      #s1 <- grad(obj.lm, para, map = map, data = data, ref = ref, inv.V = inv.V, bet0 = bet0, outcome = outcome)
    }
    
    if(family == 'binomial'){
      s0 <- score.lo(para, map, data, ref, inv.V, bet0, outcome)
      #s1 <- grad(obj.lo, para, map = map, data = data, ref = ref, inv.V = inv.V, bet0 = bet0, outcome = outcome)
    }
    
    if(all(abs(s0) < 1e-6)){
      break
    }
    
    if(family == 'gaussian'){
      t0 <- try(inv.h0 <- solve(hess.lm(para, map, data, ref, inv.V, bet0, outcome)), silent = TRUE)
    }
    
    if(family == 'binomial'){
      t0 <- try(inv.h0 <- solve(hess.lo(para, map, data, ref, inv.V, bet0, outcome)), silent = TRUE)
    }
    
    if('try-error' %in% class(t0)){
      return(list(coefficients = para.null, score = para.null, conv =0))
    }
    d0 <- as.vector(inv.h0 %*% s0)
    if(max(abs(d0)) > 1){
      #d0 <- d0/max(abs(d0))
    }
    para <- para - d0
    
    i <- i + 1
    #print(s0)
  }
  
  #svd(inv.h0)$d
  
  #print(s0)
  
  list(coefficients = para, score = s0, conv = ifelse(all(abs(s0) < 1e-6), 1, 0))
  
}
