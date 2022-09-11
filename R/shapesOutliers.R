
shapesOutliers=function(shapesRaw,class=NULL,
                        desired_levels=NULL,max_iter=2,c_IQR=3,
                        plot2d_outliers=F,joinline=c(1,2,1),
                        plot2d_procrustes=F){

  if(FALSE %in%lapply(c('shapes','pracma','plyr','roahd'),
                      require, character.only = TRUE)){
    return()
  }

  k=dim(shapesRaw)[1]  #number of landmarks
  m=dim(shapesRaw)[2]  #dimension of the landmarks
  n=dim(shapesRaw)[3]  #number of observations/shapes

  pos_WO=c() #saves the positions j in 1:n such that
  #`shapesRaw[,,j]` is an outlier


  #If `desired_levels` is NULL, then we seek outliers in the
  #whole dataset, i.e., without considering any level division. In that case,
  #`class` is set to a constant factor of length `n` equal to `All obs`
  if(is.null(desired_levels)){
    class=factor(rep('All obs',n))
    desired_levels='All obs'
  }else if (is.null(class)){
    return('Error: `desired_levels` exists while `class` does not exist')

  }else if (sum(desired_levels%in%levels(class))!=length(desired_levels)){
    return('Error: desired_levels is not well defined')
  }# (`desired_levels` must be a subset of `levels(class)`)


  #For each level in `desired_levels`, we seek outliers
  for(level in desired_levels){

    pos_WO_level=which(class==level) #same as `pos_WO`, but for each level
    shapes_level=shapesRaw[,,pos_WO_level]
    par(mfrow=c(2,2))

    #At most `max_iter` iterations are used to find outliers in each level
    for (iter in 0:max_iter) {

      n_shapes_level=dim(shapes_level)[3]
      proc=procGPA(shapes_level,distances=FALSE,scale=TRUE,
                   pcaoutput=FALSE,eigen2d=F)

      ##------------------ Plot of Procrustes coord. (level) -------------------
      if(plot2d_outliers==TRUE){

        if(iter==0){
          plot(proc$rotated[,1,],proc$rotated[,2,],axes=FALSE,asp=1,cex=0.3,
               pch=20,xlab="",ylab="",
               main=paste("Initial procrustes coordinates + mean shape (blue)
                            \n level:",level))

          for (j in 1:n_shapes_level){
            lines(proc$rotated[joinline,,j],col="grey")
          }
          lines(proc$mshape [joinline,],col=5,lwd=3)
          points(proc$rotated[,1,],proc$rotated[,2,],cex=0.3, pch=20)
        }
        plot(proc$rotated[,1,],proc$rotated[,2,],axes=FALSE,asp=1,cex=0.3,
             pch=20,xlab="",ylab="",
             main=paste("Updating procrustes coordinates + mean shape (blue)
                         \n level:",level,"| iter:",as.character(iter)))

        for (j in 1:n_shapes_level){
          lines(proc$rotated[joinline,,j],col="grey")
        }
        lines(proc$mshape [joinline,],col=5,lwd=3)
        points(proc$rotated[,1,],proc$rotated[,2,],cex=0.3, pch=20)
      }
      ##--------------------------------------------------------------------

      if(iter==max_iter){ break }


      # We flatten each Procrustes configuration matrix into a k*m dimensional
      # vector. The same applies for the mean shape. Then, we compute the maha-
      # lanobis distances^2  between each Procrustes vector and the mean vector
      sh_flat=t(apply(proc$rotated, 3, function(x) as.vector(t(x))))
      mean_flat=as.vector(t(proc$mshape))
      pseudo_inverse=pinv(cov(sh_flat))
      mahalanobis_dist=c()
      for (j in 1:n_shapes_level) { #squared mahalanobis distances
        mahalanobis_dist[j]=
          mahalanobis_dist[j]=matrix(sh_flat[j,]-mean_flat,1,k*m)%*%
          pseudo_inverse%*%matrix(sh_flat[j,]-mean_flat,k*m,1)
      }
      # Outliers are obtained the same way as in boxplots, except for the fact
      # that the constant multiplying the IQR, `c_IQR`, can be changed
      # to any value >=0 (take into account that the higher values we set in
      # `c_IQR`, the less but more extreme observations we will detect)
      q3_plus_c_IQR=  quantile(mahalanobis_dist,.75) +
        c_IQR*(quantile(mahalanobis_dist,.75)-quantile(mahalanobis_dist,.25))
      outliers_level=which(mahalanobis_dist>q3_plus_c_IQR)

      #If no outliers are found, then we jump to the next level in`desired_levels`
      if(length(outliers_level)==0){ break }

      ##----------------------- Plot of outliers (red) -------------------------
      if(plot2d_outliers==TRUE){

        for (o in outliers_level) {
          lines(proc$rotated[joinline,,o],col=2,lwd=2)
        }
        lines(proc$mshape[joinline,],col=5,lwd=5)
        points(proc$rotated[,1,],proc$rotated[,2,],cex=0.25,pch=20)
      }
      ##--------------------------------------------------------------------

      #If any outliers are found, then they need to be removed from the list of
      #shapes so that we can start finding new outliers in the next iteration
      shapes_level=shapes_level[,,-outliers_level]
      pos_WO_level=pos_WO_level[-outliers_level]

    }  #END of current level
    pos_WO=sort(c(pos_WO,pos_WO_level))
  }  #END of ALL levels


  if(length(pos_WO)==n){
    print('No outliers found')
  }


  #Plot of Procrustes coordinates (by level) before and after removing outliers
  if(plot2d_procrustes==TRUE){

    par(mar = c(1, 1, 3, 1))
    layout(matrix(c(1,2,1,2,1,2,3,3), 4, 2, byrow = TRUE))
    #If `class` has more than 10 classes, then the vector `colorss`
    #must be extended
    colorss=c('limegreen','darkgreen','darkorchid1','blue','deeppink',
    'orange','seagreen1','thistle1','turquoise2','slateblue4')
    #colorss=c('limegreen','red','darkorange','red','deeppink',
    #          'orange','seagreen1','thistle1','turquoise2','slateblue4')

    for(w in c('With','Without')){
      if(w=='With'){
        ind=which(class %in% desired_levels)
        shapesRaw_w=shapesRaw[,,ind]
        class_w=class[ind]
      }else if (w=='Without'){
        ind=which(class[pos_WO] %in% desired_levels)
        shapesRaw_w=(shapesRaw[,,pos_WO])[,,ind]
        class_w=(class[pos_WO])[ind]
      }
      proc=procGPA(shapesRaw_w,distances=FALSE,scale=TRUE,
                   pcaoutput=FALSE,eigen2d=F)

      xmin=min(proc$rotated[,1,]);xmax=max(proc$rotated[,1,]);a=0.1*(xmax-xmin)
      ymin=min(proc$rotated[,2,]);ymax=max(proc$rotated[,2,]);b=0.1*(ymax-ymin)

      plot(proc$mshape[,1],proc$mshape[,2],axes=F,asp=1,cex=1.1,col='white',
           xlim=c(xmin-a,xmax+a),ylim=c(ymin-b,ymax+b),pch=20,xlab="",ylab="",
           main=paste(w,
                      ' outliers| Procrustes coordinates & sample \n mean shape (blue), median shape (black)'))

      for (j in 1:length(class_w)) {
        points(as.matrix(proc$rotated[,,j]),pch=20,cex=1.2,
               col=colorss[as.numeric(class_w[j])])
      }
      #mean shape  (blue)
      lines(proc$mshape[joinline,],lwd=4,col='cyan')

      #median shape  (black)
      p_l=alply(proc$rotated,2,function(x) t(x))
      Proc_mfData = mfData(1:k,p_l)
      med_mfData = median_mfData(Proc_mfData)

      median_shape=matrix(0,k,m)
      for (d in 1:m) {
        median_shape[,d]=med_mfData$fDList[[d]]$values
      }
      lines(median_shape[joinline,],lwd=4,col='black')

      #legend

    }
    par(mar = c(1, 1, 1, 1))
    plot(1:3,1:3,col='white',xlab = "",ylab = "",axes=F)
    legend(c(1,3),c(1,3),
           legend=desired_levels,col=colorss[which(table(class_w)>0)],
           pch=19,box.lty=0,cex=1.5,lty=NA,lwd=1.2,
           ncol=length(desired_levels))
  }
  par(mfrow=c(1,1))
  return(without_outliers=pos_WO)

}# * END of function *





# SOME EXAMPLES

##### Dataset `apes` from `shapes` library
library(shapes)
data(apes)



#1) Finding outliers in all levels and plotting the iterations and the final
#procrustes coordinates. Now, we set the appropriate joinline between landmarks

#**** IN ORDER TO CHECK THE PERFORMANCE, LET US MANUALLY INTRODUCE AN
#*    OUTLIER !!!!!!

outlier_position=168
b=runif(n=8,min=0,max=100)
apesOUTLIER=abind::abind(apes$x, matrix(c(50*cos(b),-3*sin(5*b)),8,2,byrow = F),
                         along = 3)
plot(apesOUTLIER[,,168])
classOUTLIER=apes$group
classOUTLIER[168]='gorm'

#***
#*


shapesOutliers(shapesRaw=apesOUTLIER,class=classOUTLIER,
               desired_levels = levels(classOUTLIER),max_iter=2,c_IQR=3,
               plot2d_outliers=T,joinline=c(2,3,4,5,1,6,7,8,2),
               plot2d_procrustes=T)






#2) Not defining the `joinline` between landmarks nor the `desired_levels`
# (hence, all the shapes are considered to have the same level)
shapesOutliers(shapesRaw=apes$x,class=apes$group,
               plot2d_outliers=T,
               plot2d_procrustes=T)
#equivalent to
shapesOutliers(shapesRaw=apes$x,
               plot2d_outliers=T,
               plot2d_procrustes=T)




#3) Finding outliers in specific levels and  NOT plotting the iterations
# but YES the final proc. coord
shapesOutliers(shapesRaw=apes$x,class=apes$group,
               desired_levels = c('gorm','gorf'),
               plot2d_outliers=F,joinline=c(2,3,4,5,1,6,7,8,2),
               plot2d_procrustes=T)




#4) Considering just 1 level
shapesOutliers(shapesRaw=apes$x,class=apes$group,
               desired_levels = c('gorm'),
               plot2d_outliers=T,joinline=c(2,3,4,5,1,6,7,8,2),
               plot2d_procrustes=T)




#5) Setting `max_iter` and `c_IQR` to different values. The lower `c_IQR`, the
# more observations that we consider as outliers. Ideally, it should be equal to
# 1.5 or 3. By default `max_iter`=2 and `c_IQR`=3
shapesOutliers(shapesRaw=apes$x,class=apes$group,
               max_iter=5,c_IQR=1,desired_levels = levels(apes$group),
               plot2d_outliers=T,joinline=c(2,3,4,5,1,6,7,8,2),
               plot2d_procrustes=T)






#Another examples
data(schizophrenia)
shapesOutliers(shapesRaw=schizophrenia$x,class=schizophrenia$group,
               desired_levels =levels(schizophrenia$group),
               max_iter=2,c_IQR=1.5,
               plot2d_outliers=T,joinline=c(12,11,10,8,7,5,6,1,4,3,2,12),
               plot2d_procrustes=T)

#
data(steroids)
shapesOutliers(shapesRaw=steroids$x,
               max_iter=2,c_IQR=1.5,
               plot2d_outliers=T,joinline=c(1:6,1,6,5,4,7:10,5,4,7,11:14,8,14:17,13),
               plot2d_procrustes=T)

#
data(macaques)
shapesOutliers(shapesRaw=macaques$x,class = macaques$group,
               desired_levels = levels(macaques$group),
               max_iter=2,c_IQR=1.5,
               plot2d_outliers=T,joinline=c(1,2,5,2,3,4,1,6,5,3,7,6,4,7),
               plot2d_procrustes=T)

#
data(mice)
shapesOutliers(shapesRaw=mice$x,class = mice$group,
               desired_levels = levels(mice$group),
               max_iter=2,c_IQR=1.5,
               plot2d_outliers=T,joinline=c(1,6,2:5,1),
               plot2d_procrustes=T)
