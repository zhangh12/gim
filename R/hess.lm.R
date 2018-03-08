
hess.lm <- function(para, para.id, data, inv.V, bet0, outcome){
  
  h <- numDeriv::jacobian(score.lm, para, para.id = para.id, data = data, inv.V = inv.V, bet0 = bet0, outcome = outcome)
  colnames(h) <- names(para)
  rownames(h) <- names(para)
  h
  
}
