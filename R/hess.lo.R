
hess.lo <- function(para, map, data, ref, inv.V, bet0, outcome){
  
  h <- numDeriv::jacobian(score.lo, para, map = map, data = data, ref = ref, inv.V = inv.V, bet0 = bet0, outcome = outcome)
  colnames(h) <- names(para)
  rownames(h) <- names(para)
  h
  
}
