
hess.lm <- function(para, para.id, int, inv.V, bet0, outcome = 'y'){
  
  h <- numDeriv::jacobian(score.lm, para, para.id = para.id, int = int, inv.V = inv.V, bet0 = bet0, outcome = outcome)
  colnames(h) <- names(para)
  rownames(h) <- names(para)
  h
  
}
