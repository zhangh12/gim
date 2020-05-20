
collapse.model <- function(family, model, nsample, ncase, nctrl){
  
  nmodel <- length(model)
  if(nmodel == 1){
    return(list(model = model, nsample = nsample, ncase = ncase, nctrl = nctrl))
  }
  
  if(family == 'case-control'){
    
    rm.id <- NULL
    for(i in 1:(nmodel - 1)){
      for(j in (i+1):nmodel){
        if(ncase[i, j] == ncase[i, i] && ncase[i, j] == ncase[j, j] 
           && nctrl[i, j] == nctrl[i, i] && nctrl[i, j] == nctrl[j, j] 
           && model[[i]]$form == model[[j]]$form){
          model[[i]]$info <- rbind(model[[i]]$info, model[[j]]$info)
          model[[i]]$info <- model[[i]]$info[!duplicated(model[[i]]$info), ]
          rm.id <- c(rm.id, j)
        }
      }
    }
    
    if(!is.null(rm.id)){
      m <- list()
      k <- 0
      for(i in 1:nmodel){
        if(!(i %in% rm.id)){
          k <- k + 1
          m[[k]] <- model[[i]]
        }
      }
      model <- m
      rm(m)
      ncase <- ncase[-rm.id, -rm.id, drop = FALSE]
      nctrl <- nctrl[-rm.id, -rm.id, drop = FALSE]
    }
    
  }else{
    
    rm.id <- NULL
    for(i in 1:(nmodel - 1)){
      for(j in (i+1):nmodel){
        if(nsample[i, j] == nsample[i, i] && nsample[i, j] == nsample[j, j] 
           && model[[i]]$form == model[[j]]$form){
          model[[i]]$info <- rbind(model[[i]]$info, model[[j]]$info)
          model[[i]]$info <- model[[i]]$info[!duplicated(model[[i]]$info), ]
          rm.id <- c(rm.id, j)
        }
      }
    }
    
    if(!is.null(rm.id)){
      m <- list()
      k <- 0
      for(i in 1:nmodel){
        if(!(i %in% rm.id)){
          k <- k + 1
          m[[k]] <- model[[i]]
        }
      }
      model <- m
      rm(m)
      nsample <- nsample[-rm.id, -rm.id, drop = FALSE]
    }
    
  }
  
  list(model = model, nsample = nsample, ncase = ncase, nctrl = nctrl)
  
}
