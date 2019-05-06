

gim <- function(formula, family, data, model, 
                nsample = NULL, ncase = NULL, nctrl = NULL,
                ref = NULL, ...){
  
  UseMethod('gim')
  
}
