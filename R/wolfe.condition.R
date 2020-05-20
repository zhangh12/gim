
wolfe.condition <- function(para, map, family, data, ref, inv.V, bet0, sample.info, outcome, type, 
                            p, score, obj){
  
  p.sc <- sum(p * score)
  step <- convex.hull(para, map, family, ref, sample.info, .5^(0:20), p)
  c1 <- 1e-4
  c2 <- .9
  
  cond1 <- NULL
  cond2 <- NULL
  o <- obj
  for(a in step){
    obj1 <- compute.obj(para + a * p, map, family, data, ref, inv.V, bet0, sample.info, outcome, type)
    o <- c(o, obj1)
    cond1 <- c(cond1, ifelse(obj1 - obj >= c1 * a * p.sc, TRUE, FALSE))
    
    score1 <- compute.score(para + a * p, map, family, data, ref, inv.V, bet0, sample.info, outcome, type)
    cond2 <- c(cond2, ifelse(sum(score1 * p) <= c2 * sum(score * p), TRUE, FALSE))
  }
  
  if(all(!cond1)){
    stop('Armijo fails')
  }
  
  if(any(cond1 & cond2)){
    message('Wolfe met')
    a <- max(step[cond1 & cond2])
  }
  
  # all cond2 are FALSE
  message('Armijo met')
  a <- step[which.max(o[-1][cond1])]
  
  #print(cbind(step[1:length(cond1)], cond1, cond2, o[-1]))
  #plot(c(0, step[1:length(cond1)]), o, pch = 20)
  
  return(a)
  
}
