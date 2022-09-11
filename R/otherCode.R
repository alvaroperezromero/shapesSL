
# CODE IN ORDER TO READ limbs' data




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

################################################################################


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

# (OPTIONAL) Escoger un subset de landmarks
limb_shapes=limb_shapes[1:20,,]


# (OPTIONAL)
limb_shapes=limb_shapes[,,1:569]
class_global=class_global[1:569]
agemons=agemons[1:569]
#limb_shapes=limb_shapes[,,-c(1:569)]
#class_global=class_global[-c(1:569)]
#agemons=agemons[-c(1:569)]


# (OPTIONAL)
interval_months='[0,24)'
month_division=24
lim_sup_months=Inf
ind=which(cut(agemons,breaks=sort(c(0,lim_sup_months,month_division)),right = F)
          ==interval_months)
limb_shapes=limb_shapes[ , ,ind]
class_global=class_global[ind]


# (OPTIONAL)
list_of_levels=list(c('NOR','RIS'),c('MAM','SAM'))
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



















#paramtest::grid_search(SAM_final,
#                       params=list(which_subset=c('all'),
#                                   byage=c(F,T),
#                                   interval_months=c('[0,24)'),
#                                   registration_method=c('mean','median'),
#                                   model_name=c('knn')),

#                       n.iter=1,
#                       limb_shapes=limb_shapes,class_global=class_global,agemons=agemons,
#                       landmarks=1:20,
#                       month_division=c(24),
#                       list_of_levels=list(c('NOR','RIS'),c('MAM','SAM')),

#                       splitting_method='cv',folds=2,
#                       train_prop=0.8,n_partitions=2,
#                       remove_outliers=F,max_iter=1,c_IQR=3,
#                       allometry=T
#)
