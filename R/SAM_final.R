
#Preprocessing of data
################################################################################

# Let us set the directory
#setwd("")
# The column `Id` is deleted from the original dataset
data=read.csv("(fixed)957_2022_left_arm.csv")[,-1]
# Let us perform reflection over Y=0 of the raw coordinates (odd elements)
#in order to correct the positions of the limbs
data[,seq(2,40,2)]=-data[,seq(2,40,2)]
# We change the levels of $class_global so that:NOR < RIS < MAM < SAM
data$class_global=factor(data$class_global,levels=c("NOR","RIS","MAM","SAM"))
# n_data is the number of children in the dataset
n_data=length(data[,1])

# The age(months) of each child is saved in `agemons`
agemons=data$agemons
# The classification of each child is saved in `class_global`
class_global=data$class_global


#We transform each 40x1 observation into a 20x2 matrix, where 1,3,5,7,9,..,39
#columns correspond with the X values and 2,4,6,8,10,...,40 columns correspond
#with the Y  values. Hence, columns 1,2 are the first landmark, 3,4 the second
#landmark, and so on. The whole set of 20x2 matrices is saved into a list of
# dimension `n_data` (the nb of children) called `limb`
#Additionally, each configuration matrix is centered about the origin...

raw_coordinates=as.matrix(data[,1:40])
limb=list()
for (kid in 1:n_data) {
  limb[[kid]]=cbind(raw_coordinates[kid,seq(1,40,2)], #1,3,5,7,9,..,39
                    raw_coordinates[kid,seq(2,40,2)]) #2,4,6,8,10,...,40

  #Translation about the origin x=0,y=0 of the configuration matrix
  #(see page 63 of Dryden)
  limb[[kid]]=(diag(20)-(1/20)*matrix(1,20,1)%*%matrix(1,1,20))%*%limb[[kid]]
}
rm(raw_coordinates)

# Now, `limb_shapes` is a 3d matrix of dimension 20 x 2 x n_data
#such that limb_shapes[,,j]=limb[[j]] for each j
limb_shapes=simplify2array(limb)

################################################################################




#Final function to test lims' dataset  ***

SAM_final=function(iter=1,
              limb_shapes,class_global,agemons,
                   which_subset='all',
                   landmarks=1:20,
                   byage,month_division,interval_months,
                   list_of_levels,

                   splitting_method='cv',folds=10,
                   train_prop=0.8,n_partitions=50,
                   remove_outliers,max_iter,c_IQR,
                   registration_method,
                   allometry,
                   model_name


                   ){


  #Subset division
  if(which_subset=='all'){
    limb_shapes=limb_shapes[landmarks,,]
  }else if(which_subset=='first'){
    limb_shapes=limb_shapes[landmarks,,1:569]
    class_global=class_global[1:569]
    agemons=agemons[1:569]
  }else if(which_subset=='second'){
    limb_shapes=limb_shapes[landmarks,,-c(1:569)]
    class_global=class_global[-c(1:569)]
    agemons=agemons[-c(1:569)]
  }
  #Age division
  if(byage==T){
    lim_sup_months=Inf
    ind=which(cut(agemons,breaks=sort(c(0,lim_sup_months,month_division)),right = F)
              ==interval_months)
    limb_shapes=limb_shapes[ , ,ind]
    class_global=class_global[ind]
  }

  #Nutritional status division
  new_levels=c()
  for (l_1 in 1:length(list_of_levels)) {
    for(l_1_1 in 1:length(list_of_levels[[l_1]])){
      new_levels=c(new_levels,paste(list_of_levels[[l_1]],collapse="-"))
    }
  }
  ind=which(class_global %in% unlist(list_of_levels))
  limb_shapes=limb_shapes[,,ind]
  class_global=class_global[ind]
  class_global=droplevels(class_global)
  levels(class_global)=new_levels


  a=shapesClassification_2(
                           shapesRaw=limb_shapes,class=class_global,
                           splitting_method=splitting_method,folds=folds,
                           train_prop=train_prop,n_partitions=n_partitions,
                           remove_outliers=remove_outliers, max_iter=max_iter,c_IQR=c_IQR,
                           registration_method=registration_method,
                           allometry=allometry,
                           model_name=model_name)
  return(a$conf_matrix)

} #END OF FUNCTION



#Example

b=SAM_final(limb_shapes=limb_shapes,class_global=class_global,agemons=agemons,
            which_subset = 'All',
            byage=T,month_division=c(24),interval_months='[24,Inf)',
            list_of_levels=list(c('NOR','RIS'),c('MAM','SAM')),
            splitting_method='cv',folds=2,
            train_prop=0.7,n_partitions=50,
            remove_outliers=F,max_iter=2,c_IQR=1.5,
            registration_method='mean',
            allometry=T,
            model_name='lda'
)




###############################
############################### OBTAINING RESULTS WITH PARALLELIZATION #########

### NOR,   RIS ,  MAM,   SAM   (4CLASSES  - TODO EL DATASET)

aa=expand.grid(interval_months=c('[0,24)','[24,Inf)'),
               registration_method=c('mean','median'),
               model_name=c('knn','lda','rpart'))
aa[,1]=as.character(aa[,1])
aa[,2]=as.character(aa[,2])
aa[,3]=as.character(aa[,3])

library(foreach)
library(doParallel)

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

RESULTADOS_all_4class=foreach(i=1:12)%dopar%{


  list(paste(aa[i,3],'-',aa[i,2],'-',aa[i,1]),
       SAM_final(limb_shapes=limb_shapes,class_global=class_global,agemons=agemons,
                 which_subset = 'All',
                 byage=T,month_division=c(24),interval_months=aa[i,1],
                 list_of_levels=list('NOR','RIS','MAM','SAM'),
                 splitting_method='cv',folds=20,
                 remove_outliers=T,max_iter=1,c_IQR=1.5,
                 allometry=T,
                 registration_method=aa[i,2],
                 model_name=aa[i,3])
  )
}
stopCluster(cl)



### NOR,   RIS ,  MAM,   SAM   (4CLASSES  - 1st partition)

aa=expand.grid(interval_months=c('[0,24)','[24,Inf)'),
               registration_method=c('mean','median'),
               model_name=c('knn','lda','rpart'))
aa[,1]=as.character(aa[,1])
aa[,2]=as.character(aa[,2])
aa[,3]=as.character(aa[,3])

library(foreach)
library(doParallel)

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

RESULTADOS_1st_4class=foreach(i=1:12)%dopar%{


  list(paste(aa[i,3],'-',aa[i,2],'-',aa[i,1]),
       SAM_final(limb_shapes=limb_shapes[,,1:569],class_global=class_global[1:569],agemons=agemons[1:569],
                 which_subset = 'All',
                 byage=T,month_division=c(24),interval_months=aa[i,1],
                 list_of_levels=list('NOR','RIS','MAM','SAM'),
                 splitting_method='cv',folds=20,
                 remove_outliers=T,max_iter=1,c_IQR=1.5,
                 allometry=T,
                 registration_method=aa[i,2],
                 model_name=aa[i,3])
  )
}
stopCluster(cl)












### NOR,   RIS ,  MAM,   SAM   (2CLASSES  - TODO EL DATASET)

aa=expand.grid(interval_months=c('[0,24)','[24,Inf)'),
               registration_method=c('mean','median'),
               model_name=c('knn','lda','rpart'))
aa[,1]=as.character(aa[,1])
aa[,2]=as.character(aa[,2])
aa[,3]=as.character(aa[,3])

library(foreach)
library(doParallel)

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

RESULTADOS_all_2class=foreach(i=1:12)%dopar%{


  list(paste(aa[i,3],'-',aa[i,2],'-',aa[i,1]),
       SAM_final(limb_shapes=limb_shapes,class_global=class_global,agemons=agemons,
                 which_subset = 'All',
                 byage=T,month_division=c(24),interval_months=aa[i,1],
                 list_of_levels=list(c('NOR','RIS'),c('MAM','SAM')),
                 splitting_method='cv',folds=20,
                 remove_outliers=T,max_iter=1,c_IQR=1.5,
                 allometry=T,
                 registration_method=aa[i,2],
                 model_name=aa[i,3])
  )
}
stopCluster(cl)



### NOR,   RIS ,  MAM,   SAM   (2CLASSES  - 1st partition)

aa=expand.grid(interval_months=c('[0,24)','[24,Inf)'),
               registration_method=c('mean','median'),
               model_name=c('knn','lda','rpart'))
aa[,1]=as.character(aa[,1])
aa[,2]=as.character(aa[,2])
aa[,3]=as.character(aa[,3])

library(foreach)
library(doParallel)

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

RESULTADOS_1st_2class=foreach(i=1:12)%dopar%{


  list(paste(aa[i,3],'-',aa[i,2],'-',aa[i,1]),
       SAM_final(limb_shapes=limb_shapes[,,1:569],class_global=class_global[1:569],agemons=agemons[1:569],
                 which_subset = 'All',
                 byage=T,month_division=c(24),interval_months=aa[i,1],
                 list_of_levels=list(c('NOR','RIS'),c('MAM','SAM')),
                 splitting_method='cv',folds=20,
                 remove_outliers=T,max_iter=1,c_IQR=1.5,
                 allometry=T,
                 registration_method=aa[i,2],
                 model_name=aa[i,3])
  )
}
stopCluster(cl)
