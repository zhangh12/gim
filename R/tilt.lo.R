
# compute Delta = exp(X * theta) for internal data
# compute delta_i = exp(X_1i * alp + X_2i * bet) for ith auxiliary model
tilt.lo <- function(para, map, ref){
  
  tilt.cc(para, map, ref)
  
}
