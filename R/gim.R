

gim <- function(formula, family, data, model, 
                nsample = NULL, ncase = NULL, nctrl = NULL,
                ref = NULL, niter = 2){
  
  UseMethod('gim')
  
}
