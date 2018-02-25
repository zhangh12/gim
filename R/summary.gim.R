
summary.gim <- function(object, ...){
  
  se <- sqrt(diag(object$vcov))
  tval <- coef(object) / se
  
  TAB <- cbind(Estimate = coef(object), 
               StdErr = se, 
               t.value = tval, 
               p.value = pchisq(tval^2, df = 1, lower.tail = FALSE))
  colnames(TAB) <- c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)')
  
  res <- list(call = object$call, 
              coefficients = TAB)
  
  class(res) <- 'summary.gim'
  res
  
}

