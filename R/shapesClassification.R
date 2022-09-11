
shapesClassification=function(training_shapesRaw,class_training,
                               test_shapesRaw,class_test,
                               remove_outliers=F, max_iter=2,c_IQR=3,
                               registration_method='mean',
                               allometry=F,
                               model_name='knn'){
  
  if(FALSE %in%lapply(c('shapes','geomorph','caret'),
                          require, character.only = TRUE)){
    return()
  }
  
  k=dim(training_shapesRaw)[1] #number of landmarks
  m=dim(training_shapesRaw)[2] #dimension of each landmark
  n_test=dim(test_shapesRaw)[3] #number of configuration matrices in the test set
  
  if(k!=dim(test_shapesRaw)[1] | m!=dim(test_shapesRaw)[2] | 
     length(dim(training_shapesRaw))!=length(dim(test_shapesRaw))){
    return('Bad definition of `training_shapesRaw` or `shapesRaw`: dimensions do not match')
  }
  
  
  # 1) (optional) check for outliers
  if(remove_outliers==TRUE){
    pos_WO=shapesOutliers(shapesRaw=training_shapesRaw,class=class_training,
                           desired_levels=levels(class_training),max_iter,c_IQR,
                                          plot2d_outliers=F,plot2d_procrustes=F)
   # Updating training set
    training_shapesRaw=training_shapesRaw[,,pos_WO]
    class_training=class_training[pos_WO]
  }
  
  
  
  # 2) Obtaining the Procrustes fit in the training set
  proc=procGPA(training_shapesRaw,distances=ifelse(allometry==TRUE, TRUE,FALSE),
                                       scale=TRUE,pcaoutput=FALSE,eigen2d=FALSE) 
  
  
  
  # 3) Registration of `test_shapesRaw` on the mean shape of `proc`
  reg=shapesRegistration(shapesProc=proc$rotated,
                             mean_shape=proc$mshape,
                              shapesRaw=test_shapesRaw,
                                registration_method=registration_method)
  test_shapesProc=reg
  
  
  # 4) Flattening the training and test sets so that they are ready to 
  #    be supplied to the classifiers
  if(allometry==FALSE){
    training_flat=as.data.frame(t(apply(
                                 proc$rotated, 3, function(x) as.vector(t(x)))))
    training_flat$classification=class_training
    test_flat=as.data.frame(t(apply(
                              test_shapesProc, 3, function(x) as.vector(t(x)))))
  
  # 4) (optional) ALLOMETRY 
  }else if(allometry==TRUE){
    allom=procD.lm(rotated~log(size), data=proc)
    training_flat = as.data.frame( allom$residuals) 
    training_flat$classification=class_training
    
    test_flat=matrix(0,n_test,k*m)
    for (i in 1:n_test){
      size=centroid.size(test_shapesRaw[,,i])
      
      coord_est=allom$coefficients[1,]+allom$coefficients[2,]*log(size)
      coord_est=matrix(coord_est,k,m,byrow = TRUE)
      
      residual=test_shapesProc[,,i]-coord_est
      test_flat[i,]=as.vector(t(residual))
    }
  }

  
  
  # 5) Fitting classifiers and testing
  set.seed(1)
  start.time = Sys.time()
  model_fitted=train(classification ~ ., data = training_flat,method=model_name)
  end.time = Sys.time()
  
  predicted=predict(model_fitted, test_flat)
  return(list(actual_class=class_test,predicted=predicted,
            time_training=end.time - start.time))
  
}





shapesClassification_2=function(
                                shapesRaw,class,
                                splitting_method='cv',folds=10,
                                train_prop=0.8,n_partitions=10,
                                remove_outliers=F, max_iter=2,c_IQR=3,
                                registration_method='mean',
                                allometry=F,
                                model_name='knn'){
  
  if(FALSE %in%lapply(c('caret'),require, character.only = TRUE)){
    return()
  }
  
  
  c_actual_class= factor(c(),levels=levels(class))
  c_predicted= factor(c(),levels=levels(class))
  times_training=c()
  if (splitting_method=='cv'){
    partitions=createFolds(class,k=folds,returnTrain = T,list=T)
  }else if (splitting_method=='rp'){
    partitions=caret::createDataPartition(class,p=train_prop,list=TRUE,times=n_partitions)
  }else{
    return('`splitting_method` is not well defined')
  }

  
  for (p in 1:length(partitions)) {
    
    a=shapesClassification(
                            training_shapesRaw=shapesRaw[,,partitions[[p]]],
                            class_training=class[partitions[[p]]],
                            test_shapesRaw=shapesRaw[,,-partitions[[p]]],
                            class_test=class[-partitions[[p]]],
                            remove_outliers=remove_outliers,
                            max_iter=max_iter,
                            c_IQR=c_IQR,
                            registration_method=registration_method,
                            allometry=allometry,
                            model_name=model_name)
    
    c_actual_class=c(c_actual_class,a$actual_class)
    c_predicted=c(c_predicted,a$predicted)
    times_training=c(times_training,a$time_training)
    print(p)
    
  }
  print('Finished estimation')
  return(list(
    conf_matrix=confusionMatrix(data=c_predicted,reference=c_actual_class),
    summary_time_training=summary(times_training)))
}







### Example

data(apes)

training_shapesRaw=apes$x[,,c(1:25,30:55)]
class_training=droplevels(apes$group[c(1:25,30:55)])
test_shapesRaw=apes$x[,,c(26:29,56:59)]
class_test=droplevels(apes$group[c(26:29,56:59)])

a=shapesClassification(training_shapesRaw=training_shapesRaw,
                        class_training=class_training,
                        test_shapesRaw=test_shapesRaw,
                        class_test=class_test,
                               remove_outliers=T, max_iter=2,c_IQR=3,
                               registration_method='median',
                               allometry=F,
                               model_name='rf')
a$actual_class
a$predicted
confusionMatrix(a$actual_class,a$predicted)
