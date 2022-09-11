

### -----> *** UNTIDY ** <-------






# Let us set the directory
#setwd("")
# The column `Id` is deleted from the original dataset
data=read.csv("(fixed)957_2022_left_arm")[,-1]
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


limb_shapes=simplify2array(limb)


limb_sizes=procGPA(limb_shapes)$size



# Now, `limb_shapes` is a 3d matrix of dimension 20 x 2 x n_data
#such that limb_shapes[,,j]=limb[[j]] for each j

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






plot((limb[[1]]/limb_sizes[1]),col="white",xlim=c(0.1,1),ylim=c(-0.3,-0.1),asp=1)


kid=949

limb_shapes[,2,kid]=-limb_shapes[,2,kid]

plot(limb_shapes[,,kid],asp=1,col='white',xlab = "",ylab = "",main='937rd child')
lines(limb_shapes[,,kid][c(20,3,19,5,7,9,17,11,13,15,1,4,16,14,12,18,10,8,6,20),],lwd=1,col=2)
text(limb_shapes[,,kid],labels=as.character(1:20),pos=3,cex=1.2)
points(limb_shapes[,,kid],lwd=4,col=1,pch=20)


par(mar=c(5.1, 4.1, 4.1, 2.1))
par(mfrow=c(1,1))




kid=10
lines((limb[[10]]/limb_sizes[10])[c(20,3,19,5,7,9,17,11,13,15,1,4,16,14,12,18,10,8,6,20),],lwd=4,col=2)
points(limb[[10]]/limb_sizes[10],lwd=4,col=2,pch=20)
text(limb[[10]]/limb_sizes[10],labels=as.character(1:20),pos=3,cex=1.2)


kid=1

#plot(0,0,col="white",xlim=c(min(limb[[kid]][,1])-5,max(limb[[kid]][,1])+5),ylim=c(min(limb[[kid]][,2])-5,max(limb[[kid]][,2])+5),xlab="X",ylab="Y",main="")
plot((limb[[kid]]/limb_sizes[kid])[c(20,3,19,5,7,9,17,11,13,15,1,4,16,14,12,18,10,8,6,20),],asp=1,lwd=4,col=1,pch=20,xlab="X",ylab="Y",main="")
points((limb[[kid]]/limb_sizes[kid])[c(2,2),],lwd=4,col=1,pch=20)

points(colMeans(limb[[kid]])[1],colMeans(limb[[kid]])[2],col="red",cex=2,pch=20)
text(colMeans(limb[[kid]])[1],colMeans(limb[[kid]])[2],labels="centroid",pos=4,cex=1.2,col = "red")

text(limb[[kid]]/limb_sizes[kid],labels=as.character(1:20),pos=1,cex=1.2)

limb[[kid]]

lines(rbind(limb[[kid]][18,],colMeans(limb[[kid]])),type="l",lty=2,lwd=2)
lines(rbind(limb[[kid]][2,],colMeans(limb[[kid]])),type="l",lty=2,lwd=2)
lines(rbind(limb[[kid]][10,],colMeans(limb[[kid]])),type="l",lty=2,lwd=2)
lines(rbind(limb[[kid]][9,],colMeans(limb[[kid]])),type="l",lty=2,lwd=2)
points(colMeans(limb[[kid]])[1],colMeans(limb[[kid]])[2],col="red",cex=2.5,pch=20)



par(mfrow = c(1,1))
plot(limb[[1]]%*%matrix(c(0,1,-1,0),2,2,byrow = TRUE))

plot(limb[[573]])


kid=572
kid2=113

plot(((diag(20)-(1/20)*matrix(1,20,1)%*%matrix(1,1,20))%*%limb[[kid2]]),col="white",ylab="",asp=1,xlab="")

lines((((diag(20)-(1/20)*matrix(1,20,1)%*%matrix(1,1,20))%*%limb[[kid2]]))[c(20,3,19,5,7,9,17,11,13,15,1,4,16,14,12,18,10,8,6,20),],lwd=4,col=1)
points((((diag(20)-(1/20)*matrix(1,20,1)%*%matrix(1,1,20))%*%limb[[kid2]]))[c(20,3,19,5,7,9,17,11,13,15,1,4,16,14,12,18,10,8,6,20),],lwd=4,col=1,cex=1.5,pch=20)
points((((diag(20)-(1/20)*matrix(1,20,1)%*%matrix(1,1,20))%*%limb[[kid2]]))[c(2,2),],lwd=4,col=1,cex=1.5,pch=20)

lines((((diag(20)-(1/20)*matrix(1,20,1)%*%matrix(1,1,20))%*%limb[[kid]]))[c(20,3,19,5,7,9,17,11,13,15,1,4,16,14,12,18,10,8,6,20),],lwd=4,col=2)
points((((diag(20)-(1/20)*matrix(1,20,1)%*%matrix(1,1,20))%*%limb[[kid]]))[c(20,3,19,5,7,9,17,11,13,15,1,4,16,14,12,18,10,8,6,20),],lwd=4,col=2,cex=1.5,pch=20)
points((((diag(20)-(1/20)*matrix(1,20,1)%*%matrix(1,1,20))%*%limb[[kid]]))[c(2,2),],lwd=4,col=2,cex=1.5,pch=20)
text((((diag(20)-(1/20)*matrix(1,20,1)%*%matrix(1,1,20))%*%limb[[kid2]])),labels=as.character(1:20),pos=3,cex=1.2)


aa=procOPA(((diag(20)-(1/20)*matrix(1,20,1)%*%matrix(1,1,20))%*%limb[[kid2]]),((diag(20)-(1/20)*matrix(1,20,1)%*%matrix(1,1,20))%*%limb[[kid]]))
lines(aa$Bhat[c(20,3,19,5,7,9,17,11,13,15,1,4,16,14,12,18,10,8,6,20),],lwd=4,col=2)
points(aa$Bhat[c(20,3,19,5,7,9,17,11,13,15,1,4,16,14,12,18,10,8,6,20),],lwd=4,col=2,cex=1.5,pch=20)
points(aa$Bhat[c(2,2),],lwd=4,col=2,cex=1.5,pch=20)
text(((diag(20)-(1/20)*matrix(1,20,1)%*%matrix(1,1,20))%*%limb[[kid2]]),labels=as.character(1:20),pos=3,cex=1.2)




pos=intersect(which(data$class_global=="NOR"),which(data$agemons<24))
plot(limb_shapes[,1,pos],limb_shapes[,2,pos],axes=TRUE,asp=1,cex=0.3,xlab="",ylab="", pch=20)
joinline=c(20,3,19,5,7,9,17,11,13,15,1,4,16,14,12,18,10,8,6,20)
for (i in 1:length(pos)){
  lines(limb_shapes[joinline,,pos[i]],col="grey")
}
points(limb_shapes[,1,pos],limb_shapes[,2,pos],cex=0.8,lwd=0.5, pch=20)
lines(aa$mshape[c(20,3,19,5,7,9,17,11,13,15,1,4,16,14,12,18,10,8,6,20),],lwd=5,col='cyan')
points(aa$mshape,cex=1,lwd=3, pch=20,col='white')

aa=procGPA(limb_shapes[,,pos])

plot(aa$rotated[,1,],aa$rotated[,2,],axes=TRUE,asp=1,cex=0.3,xlab="",ylab="", pch=20)
joinline=c(20,3,19,5,7,9,17,11,13,15,1,4,16,14,12,18,10,8,6,20)
for (i in 1:length(pos)){
  lines(aa$rotated[joinline,,i],col="grey")
}
points(aa$rotated[,1,],aa$rotated[,2,],cex=0.8,lwd=0.5, pch=20)
lines(aa$mshape[c(20,3,19,5,7,9,17,11,13,15,1,4,16,14,12,18,10,8,6,20),],lwd=5,col='cyan')
points(aa$mshape,cex=1,lwd=3, pch=20,col='white')



plot((1.1*limb_shapes[joinline,,1])%*%matrix(c(.7,.7,-.7,.7),2,2,byrow = T),xlim=c(-150,150),ylim = c(-70,115),asp=1,col='white',xlab="",ylab="")

points(aa$rotated[,1,],aa$rotated[,2,],axes=TRUE,asp=1,cex=0.3,xlab="",ylab="", pch=20)
for (i in 1:length(pos)){
  lines(aa$rotated[joinline,,i],col="grey")
}
points(aa$rotated[,1,],aa$rotated[,2,],cex=0.8,lwd=0.5, pch=20)
lines(aa$mshape[c(20,3,19,5,7,9,17,11,13,15,1,4,16,14,12,18,10,8,6,20),],lwd=5,col='cyan')
points(aa$mshape,cex=1,lwd=3, pch=20,col='white')

points((1.1*limb_shapes[,,1])%*%matrix(c(cos(pi/6),.5,-.5,cos(pi/6)),2,2,byrow = T)+matrix(30,20,2),cex=1.5,lwd=0.5, pch=20,col=2)
lines(((1.1*limb_shapes[joinline,,1])%*%matrix(c(cos(pi/6),.5,-.5,cos(pi/6)),2,2,byrow = T)+matrix(30,20,2)),col=2,lwd=2.7)

registration=procOPA(aa$mshape,limb_shapes[,,1])$Bhat
lines(registration[joinline,],col=2,lwd=3)






joinline=c(20,3,19,5,7,9,17,11,13,15,1,4,16,14,12,18,10,8,6,20)
aa=procGPA(limb_shapes[,,1:3])
aa=aa$rotated
aa[,,1]=10*aa[,,1]

aa[,,2]=5*aa[,,2]
aa[,,3]=4*aa[,,3]

plot(aa[,,1],asp=1,col='white',xlab="",ylab="")
lines(aa[joinline,,1],col='grey',lwd=2.7)
points(aa[,,1],cex=0.9,lwd=1.5, pch=20)

lines(aa[joinline,,2],col='grey',lwd=2.7)
points(aa[,,2],cex=0.9,lwd=1.5, pch=20)

lines(aa[joinline,,3],col='grey',lwd=2.7)
points(aa[,,3],cex=0.9,lwd=1.5, pch=20)


lines(aa[joinline,,2],col=1,lwd=3)
points(aa[,,2],cex=1,lwd=0.7, pch=19)


aa=procGPA(aa)
lines(aa$mshape[joinline,],col='cyan',lwd=5)
points(aa$mshape,cex=1.2,lwd=0.7, pch=20)

lines(aa[joinline,,2],col=1,lwd=3)
points(aa[,,2],cex=1.4,lwd=0.7, pch=20)

c=8*aa[,,1]%*%matrix(c(cos(pi/6),.5,-.5,cos(pi/6)),2,2,byrow = T)
lines(c[joinline,],col=2,lwd=4)
points(c,cex=1.4,lwd=0.7, pch=20)


meann=shapesRegistration(aa,shapesRaw=c,
                            registration_method='mean')


lines((meann$registration_mean)[joinline,,1],col=1,lty=2)


points(aa[,,2])
points(aa[,,3])
plot(aa,asp=1)
