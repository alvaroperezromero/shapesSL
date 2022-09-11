
shapesRegistration=function(shapesProc,mean_shape=NULL,shapesRaw,
                            registration_method){
  
  if(FALSE %in%lapply(c('shapes','plyr','roahd'),
                      require, character.only = TRUE)){
    return()
  }
                                #'mean','median','macro'
  if(registration_method %in% c('mean','median')==FALSE){
    return('`registration_method` is not well defined')
  }
  k=dim(shapesProc)[1]  #number of landmarks
  m=dim(shapesProc)[2]  #dimension of each landmark
  n=dim(shapesProc)[3]  #number of Procrustes configuration matrices
  
  # If `shapesRaw` is just a single matrix, i.e, a 2d array, then we 
  # transform it into a 3d array so that the dimensions match with `shapesProc`
  if(length(dim(shapesRaw))==2){
    shapesRaw=simplify2array(list(shapesRaw))
  }
  n_test=dim(shapesRaw)[3]
  
  # Let us check whether the dimensions match
  if(k!=dim(shapesRaw)[1] | m!=dim(shapesRaw)[2] | 
                            length(dim(shapesProc))!=length(dim(shapesRaw))){
    return('Bad definition of `shapesProc` or `shapesRaw`: dimensions do not match')
  }
  
  ## MEAN
  if(registration_method=='mean'){
    
    if(is.null(mean_shape)){
      mean_shape=rowMeans(shapesProc, dims = 2)
    }
    
    registration_mean=list()
    for (i in 1:n_test) {
      registration_mean[[i]]=procOPA(mean_shape,shapesRaw[,,i],scale=T)$Bhat
    }
    registration_mean=simplify2array(registration_mean)
    
    return(registration_mean=registration_mean)
  }
  
  ## MEDIAN
  if(registration_method=='median'){
    
    p_l=alply(shapesProc,2,function(x) t(x))
    Proc_mfData = mfData(1:k,p_l)
    med_mfData = median_mfData(Proc_mfData)
    
    median_shape=matrix(0,k,m)
    for (d in 1:m) {
      median_shape[,d]=med_mfData$fDList[[d]]$values
    }
    
    registration_median=list()
    for (i in 1:n_test) {
      registration_median[[i]]=procOPA(median_shape,shapesRaw[,,i],scale=T)$Bhat
    }
    registration_median=simplify2array(registration_median)
    
    
    return(registration_median=registration_median)
  }
  
 
}




##### Example

data(apes)

shapesProc=procGPA(apes$x[,,c(1:5,30:35)],distances=FALSE,scale=TRUE,pcaoutput=FALSE,eigen2d=F)$rotated
shapesRaw=apes$x[,,36]

###
start.time = Sys.time()##
a1=shapesRegistration(shapesProc,mean_shape=NULL,shapesRaw,
                            registration_method='mean')
end.time = Sys.time()##
end.time-start.time

###
start.time = Sys.time()##
a2=shapesRegistration(shapesProc,mean_shape=NULL,shapesRaw,
                      registration_method='median')
end.time = Sys.time()##
end.time-start.time


