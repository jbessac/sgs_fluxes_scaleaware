
rm(list=ls())

library('fields')
library('maps')
library('ggplot2')
library('pals')
library('fda')

source('functions.r')

indnx <- c(3,5,7,9,11,14,21,28,35,42,49)
indnxn <- indnx/max(indnx)
n <- length(indnx)
scalenx <- round(indnx*4/111.11,2)

load(file='lat_list_res.RData')  #lat
load(file='long_list_res.RData') #long
load('index_x_mv_allscales.RData') #index_x_mv
load('index_y_mv_allscales.RData') #index_y_mv
d1 <- length(index_x_mv[[1]])
d2 <- length(index_y_mv[[1]])

#load(file='residuals_noprcp.RData')   #res
#res0 <- res
#load(file='residuals_scale.RData') # res #saved for nx=(5,7,10)
#res3 <- res
#rm(res)
load(file='residuals.RData')    #res
load('prcp.RData')
load('Ures.RData')   # Ur
load('Utrue.RData')  

load(file='long_mv.RData') # long_mv
load(file='lat_mv.RData') # lat_mv

load(file='optim_corr_spacetime_scale/par1_phys.RData') # par1 in parameters space # not from scale-aware regression
load(file='optim_corr_spacetime_scale/par2_phys.RData')
load(file='optim_corr_spacetime_scale/par3_phys.RData')
load(file='optim_corr_spacetime_scale/par4_phys.RData')
load(file='optim_corr_spacetime_scale/par5_phys.RData')


###---------------------------------- plot functions 

map_plot <- function(x,y,z,zlim,title,sc,sc.main,sub,legend=FALSE,xl=1,yl=1,xlim,col=parula(400),breaks=NULL){
  if (missingArg(zlim)){zlim=range(z,na.rm=TRUE)}
  #if (missingArg(xlim)){xlim=range(x,na.rm=TRUE)}
  image(x,y,z,zlim=zlim,main=title,sub=sub,xlab='',ylab='',col=col,breaks=breaks,cex.sub=sc,cex.axis=sc,cex=sc,cex.lab=sc,cex.main=sc.main,mgp=c(4,1.7,0),asp=1,xlim,ylim=range(y,na.rm=TRUE),yaxt='n')
  map("world",add=TRUE,col="grey",border='grey',fill=T)
  map('world',interior=F,add=T)
  if (xl == 1){
    mtext('Longitude',side=1,line=4,cex=sc) }
  if (yl == 1){
    mtext('Latitude',side=2,line=5,cex=sc) }
  if (legend==TRUE){
    image.plot(legend.only=T,col=col,breaks=breaks,zlim=zlim,cex.lab=sc,cex=sc,cex.axis=sc,legend.width=4,legend.cex=2.5,axis.args=list(cex.axis=sc))
  }
}

map_plot0 <- function(x,y,z,zlim,title,sc,sc.main,legend=FALSE,xl=1,yl=1){
  if (missingArg(zlim)){zlim=range(z,na.rm=TRUE)}
  image(x,y,z,zlim=zlim,main=title,xlab='',ylab='',col=parula(400),xaxt='n',yaxt='n',cex.axis=sc,cex=sc,cex.lab=sc,cex.main=sc.main,mgp=c(4,1.7,0),asp=1,xlim=range(x,na.rm=TRUE),ylim=range(y,na.rm=TRUE))
  map("world",add=TRUE,col="grey",border='grey',fill=T)
  map('world',interior=F,add=T)
  if (xl == 1){
    mtext('Longitude',side=1,line=6,cex=sc) }
  if (yl == 1){
    mtext('Latitude',side=2,line=6,cex=sc) }
  if (legend==TRUE){
    image.plot(legend.only=T,col=parula(400),zlim=zlim,cex.lab=sc,cex=sc,cex.axis=sc,legend.width=4,legend.cex=5,axis.args=list(cex.axis=sc))
  }
}

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(1, 0, 1, 0),mar=c(5, 0, 0, 0),new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)  }

break 

###---------------------------------- mse local VS global regression - fig 1
load(file='residuals_movingwindow.RData')  # res1 (local)
load(file='residuals.RData')  #res (global)

i <- 4
border = unlist(index_x_mv[[i]][d1])
lon_axis = round(long[[i]],1)          
lat_axis = round(lat[[i]],1)
mse_res0 <- apply(X=(res1[[i]]^2),FUN=mean,MARGIN=1:2)
mse_res1 <- apply(X=(res[[i]]^2),FUN=mean,MARGIN=1:2)

## mse and resolved precip (fig1 - main manuscript)
zlim0 = range((mse_res1-mse_res0)/mse_res0,na.rm=TRUE)
z0 = (mse_res1-mse_res0)/mse_res0
setEPS()      
postscript('map_mse_localVSglobal_regression.eps',height=4.8,width=27)
par(mar=c(25,20,6,30),oma=c(2,1,.5,8),mai=c(.8,1.2,.6,1.4))
layout(matrix(1:2,1,2))
rbPal <- colorRampPalette(c('blue', 'gray96', 'red')) 
nc <- 100 
max_absolute_value <- max(abs(zlim0)) 
brks <- seq(-max_absolute_value,max_absolute_value,length.out=nc+1)
Col <- rbPal(nc)
n_in_class <- hist(z0, breaks=brks, plot=F)$counts>0
col_to_include <- min(which(n_in_class==T)):max(which(n_in_class==T))
brks_to_include <- min(which(n_in_class==T)):(max(which(n_in_class==T))+1)
image(lon_axis,lat_axis, z0, col=Col[col_to_include], breaks=brks[brks_to_include],main='Relative difference in MSE',cex=2.8,cex.axis=2.8,cex.lab=2.8,cex.main=3.6,xlab='',ylab='',mgp=c(4,1.7,0),asp=1)
image.plot(legend.only=TRUE,col=Col[col_to_include],breaks=brks[brks_to_include],zlim=c(-16,max(zlim0)),legend.width=4,legend.cex=2.8,axis.args=list(cex.axis=2.8))
mtext('Longitude',side=1,line=4,cex=2.8) 
mtext('Latitude',side=2,line=4,cex=2.8) 
map("world",add=TRUE,col="grey",border='grey',fill=T)
map('world',interior=F,add=T)
#
zp <- 86.4*apply(X=P[[i]],FUN=mean,MARGIN=1:2)
map_plot(lon_axis,lat_axis,zp,zlim=range(zp,na.rm=TRUE),title="Coarse-grained precipitation (mm/day)",sub=NULL,sc=2.8,sc.main=3.6,legend=TRUE,yl=0,xl=1)
dev.off()
#


###---------------------------------- R2 versus scale - fig 2
##---- graphs for n=2 (main manuscript)
load(file='residuals_noprcp.RData')   #res
res0 <- res
load(file='residuals.RData')  #res

R2 <- matrix(0,1,11)
R02 <- matrix(0,1,11)
R12 <- matrix(0,1,11)
for (i in 1:length(indnx)){
  err0 <- Ut[[i]] - Ur[[i]]
  err0[err0==0] <- NA 
  err <- log(err0)
  rm(err0)
  R2[i] <- 1-var(c(res[[i]]),na.rm=TRUE)/var(c(err),na.rm=TRUE)
  R02[i] <- 1-var(c(res0[[i]]),na.rm=TRUE)/var(c(err),na.rm=TRUE) }
rm(res0,res)


sc <- 2.4
setEPS()
postscript(file='R2_scale_global.eps',width=10,height=8)
par(mar=c(6,7,4.5,2.5))
plot(x=round((indnx*4/111.11),2),R2,pch=20,xlab='Scale N (degree)',ylab='',cex=2.2,cex.lab=sc,cex.axis=sc,ylim=range(R2,R02,R12))
points(round((indnx*4/111.11),2),R02,pch=20,cex=2.2,col=2)
mtext(expression(R^2),side=2,line=4,cex=sc)
#points(round((indnx*4/111.11),2),R12,pch=20,cex=.3,col=3)
grid(col=1,lwd=.8)
legend(y=.32,x=.3,c('Regression with precipitation','Regression without precipitation'),col=1:2,pch=20,cex=sc,bty='n')
dev.off()



##---- graphs for n=1 (supplemental material)
load(file='residuals_noprcp_n1.RData')   #res
res0 <- res
load(file='residuals_n1.RData')  #res

R2 <- matrix(0,1,11)
R02 <- matrix(0,1,11)
R12 <- matrix(0,1,11)
for (i in 1:length(indnx)){
  err0 <- Ut[[i]] - Ur[[i]]
  err0[err0==0] <- NA 
  err <- log(err0)
  rm(err0)
  R2[i] <- 1-var(c(res[[i]]),na.rm=TRUE)/var(c(err),na.rm=TRUE)
  R02[i] <- 1-var(c(res0[[i]]),na.rm=TRUE)/var(c(err),na.rm=TRUE) }
rm(res0,res)

sc <- 2.4
setEPS()
postscript(file='R2_scale_global_n1.eps',width=10,height=8)
par(mar=c(6,7,4.5,2.5))
plot(x=round((indnx*4/111.11),2),R2,pch=20,xlab='Scale N (degree)',ylab='',cex=2.2,cex.lab=sc,cex.axis=sc,ylim=range(R2,R02,R12))
points(round((indnx*4/111.11),2),R02,pch=20,cex=2.2,col=2)
mtext(expression(R^2),side=2,line=4,cex=sc)
grid(col=1,lwd=.8)
legend(y=.18,x=.3,c('Regression with precipitation','Regression without precipitation'),col=1:2,pch=20,cex=sc,bty='n')
dev.off()

###---------------------------------- regression coefficient - fig 3
##---- graphs for n=2 (main manuscript)
load(file='coeff_reg.RData')
load(file='regression_scaleaware_nx5710.RData') #scale.lm
coeff <- scale.lm$coefficients
thetaN <- matrix(coeff,8,3,byrow=TRUE)
parN <- thetaN %*% t(cbind(rep(1,n),indnx,indnx^2))
parN[1,] <- coeff[1] + coeff[2]*log(indnx) + coeff[3]*indnx*indnx

sc = 2.4
setEPS()
postscript(file='regression_coeff_scaleaware_nx5710.eps',width=15,height=8)
layout(matrix(1:2,1,2))
par(mar=c(6,6,4.5,2.5))
matplot(scalenx,t(coeff_reg[1:4,]),typ='b',pch=20,col=1:4,lty=1,xlab='Scale N (degree)',ylab='',main='Coefficients A',cex.lab=sc,cex.axis=sc,cex.main=sc,lwd=2)
mtext('Regression coefficient',side=2,line=4,cex=sc)
for (i in 1:4){
  lines(scalenx,parN[i,],col=i,lty=3,lwd=2,pch=15,type="b")
  points(scalenx[c(5,7,10)],parN[i,c(5,7,10)],col=i,pch=4,cex=3,lwd=2)}
legend(x='bottomright',lty=c(1,2,NA,1,1,1,1),pch=c(20,15,4,NA,NA,NA,NA),c('Single-scale regression','Scale-aware regression','Training scales',expression(A[0] (B[1])),expression(A[1] (B[2])),expression(A[2] (B[3])),expression(A[3] (B[4]))),col=c(1,1,1,1:4),cex=2,lwd=2,bty='n')
#
matplot(scalenx,t(coeff_reg[c(5,8,7,6),]),typ='b',pch=20,col=1:4,lty=1,xlab='Scale N (degree)',ylab='',main='Coefficients B',cex.lab=sc,cex.axis=sc,cex.main=sc,lwd=2)
for (i in 1:4){
  lines(scalenx,parN[(i+4),],col=i,lty=3,lwd=2,pch=15,type="b") 
  points(scalenx[c(5,7,10)],parN[(i+4),c(5,7,10)],col=i,pch=4,cex=3,lwd=2)}
dev.off()


##---- graphs for n=1 (supplemental material)
load(file='coeff_reg_n1.RData')
load(file='regression_scaleaware_nx5710_n1.RData') #scale.lm
coeff <- scale.lm$coefficients
thetaN <- matrix(coeff,8,3,byrow=TRUE)
parN <- thetaN %*% t(cbind(rep(1,n),indnxn,indnxn^2))
parN[1,] <- coeff[1] + coeff[2]*log(indnxn) + coeff[3]*indnxn*indnxn

sc <- 2.4
setEPS()
postscript(file='regression_coeff_scaleaware_nx5710_n1.eps',width=15,height=8)
layout(matrix(1:2,1,2))
par(mar=c(6,6,4.5,2.5))
matplot(scalenx,t(coeff_reg[1:4,]),typ='b',pch=20,col=1:4,lty=1,xlab='Scale N (degree)',ylab='',main='Coefficients A',cex.lab=sc,cex.axis=sc,cex.main=sc,lwd=2)
mtext('Regression coefficient',side=2,line=4,cex=sc)
for (i in 1:4){
  lines(scalenx,parN[i,],col=i,lty=3,lwd=2,pch=15,type="b")
  points(scalenx[c(5,7,10)],parN[i,c(5,7,10)],col=i,pch=4,cex=3,lwd=2)}
legend(x='bottomright',lty=c(1,2,NA,1,1,1,1),pch=c(20,15,4,NA,NA,NA,NA),c('Single-scale regression','Scale-aware regression','Training scales',expression(A[0] (B[1])),expression(A[1] (B[2])),expression(A[2] (B[3])),expression(A[3] (B[4]))),col=c(1,1,1,1:4),cex=2,lwd=2,bty='n')
#
matplot(scalenx,t(coeff_reg[c(5,8,7,6),]),typ='b',pch=20,col=1:4,lty=1,xlab='Scale N (degree)',ylab='',main='Coefficients B',cex.lab=sc,cex.axis=sc,cex.main=sc,lwd=2)
for (i in 1:4){
  lines(scalenx,parN[(i+4),],col=i,lty=3,lwd=2,pch=15,type="b") 
  points(scalenx[c(5,7,10)],parN[(i+4),c(5,7,10)],col=i,pch=4,cex=3,lwd=2)}
dev.off()



###---------------------------------- mean and SD of residuals - fig 4
load(file='coeff_reg.RData')
res1 <- NULL  # scale-specific regression residuals
for (i in 1:length(indnx)){
  err0 <- Ut[[i]] - Ur[[i]]
  err0[err0==0] <- NA 
  err <- log(err0)
  rm(err0)
  x1 <- log(Ur[[i]])
  x2 <- P[[i]]
  res1[[i]] <- err - coeff_reg[1,i] - coeff_reg[2,i]*x1 - coeff_reg[3,i]*x1^2 - coeff_reg[4,i]*x1^3 - coeff_reg[5,i]*x2 - coeff_reg[6,i]*x2^(1/4) - coeff_reg[7,i]*x2^(1/2) - coeff_reg[8,i]*x2^(3/4)
  rm(err,x1,x2)  }
load(file='residuals_scale.RData') # res #saved for nx=(5,7,10)
res3 <- res
rm(res)

ii0 <- c(3,6,11) 
mu1 <- apply(X=res1[[ii0[1]]],FUN=mean,MARGIN=1:2)
mu3 <- apply(X=res3[[ii0[1]]],FUN=mean,MARGIN=1:2)
mu4 <- apply(X=res1[[ii0[2]]],FUN=mean,MARGIN=1:2)
mu6 <- apply(X=res3[[ii0[2]]],FUN=mean,MARGIN=1:2)
mu7 <- apply(X=res1[[ii0[3]]],FUN=mean,MARGIN=1:2)
mu9 <- apply(X=res3[[ii0[3]]],FUN=mean,MARGIN=1:2)
zlim <- range(mu1,mu3,mu4,mu6,mu7,mu9,na.rm=TRUE)

s1 <- apply(X=res3[[ii0[1]]],FUN=sd,MARGIN=1:2)
s2 <- apply(X=res3[[ii0[2]]],FUN=sd,MARGIN=1:2)
s3 <- apply(X=res3[[ii0[3]]],FUN=sd,MARGIN=1:2)
zlim0 <- range(s1,s2,s3,na.rm=TRUE)

i <- 7
zp <- 86.4*apply(X=P[[i]],FUN=mean,MARGIN=1:2)
lon_prcp <- round(long[[i]],1)          
lat_prcp <- round(lat[[i]],1)
zlevels <- c(1e-3,1e-2,5e-2)

xlim <- c(41,180)
xlim2 <- c(42,180)
fig4 <- 1
if(fig4 == 1){
setEPS()      
postscript('map_mean_sd_res_scaleregression_nx3611.eps',width=37,height=14)
layout(matrix(1:9,3,3,byrow=TRUE))
par(mar=c(35,10,20,11),oma=c(9,.5,5,3),mai=c(.8,1.8,.6,1.2))
p <- 1
border = unlist(index_x_mv[[ii0[p]]][d1])
lon_axis = round(long[[ii0[p]]],1)          
lat_axis = round(lat[[ii0[p]]],1)
title1 = paste('N=',round(scalenx[ii0[p]],2))
title3 = paste('N=',round(scalenx[ii0[p]],2))
title2 = paste('N=',round(scalenx[ii0[p]],2))
map_plot(lon_axis,lat_axis,mu1,zlim,title=title1,sc=3.5,sc.main=5,legend=FALSE,yl=1,xl=0,sub=NULL,xlim=xlim)
contour(lon_prcp,lat_prcp,zp,col=c(rgb(240/255,240/255,240/255),rgb(130/255,130/255,130/255),rgb(15/255,15/255,15/255)),add=TRUE,labcex=1.8,lwd=2,level=zlevels,method='edge',drawlabels=FALSE)
axis(2,at=c(-10,0,10),labels=c(-10,0,10),cex.lab=3.5,cex=3.5,cex.axis=3.5)
map_plot(lon_axis,lat_axis,mu3,zlim,title3,sc=3.5,sc.main=5,legend=TRUE,xl=0,yl=0,sub=NULL,xlim=xlim)
axis(2,at=c(-10,0,10),labels=c(-10,0,10),cex.lab=3.5,cex=3.5,cex.axis=3.5)
map_plot(lon_axis,lat_axis,s1,zlim0,title2,sc=3.5,sc.main=5,legend=TRUE,yl=0,xl=0,sub=NULL,xlim=xlim)
axis(2,at=c(-10,0,10),labels=c(-10,0,10),cex.lab=3.5,cex=3.5,cex.axis=3.5)
#
p <- 2
border = unlist(index_x_mv[[ii0[p]]][d1])
lon_axis = round(long[[ii0[p]]],1)          
lat_axis = round(lat[[ii0[p]]],1)
title1 = paste('N=',round(scalenx[ii0[p]],2))
title3 = paste('N=',round(scalenx[ii0[p]],2))
title2 = paste('N=',round(scalenx[ii0[p]],2))
map_plot(lon_axis,lat_axis,mu4,zlim,title1,sc=3.5,sc.main=5,legend=FALSE,yl=1,xl=0,sub=NULL,xlim=xlim)
axis(2,at=c(-10,0,10),labels=c(-10,0,10),cex.lab=3.5,cex=3.5,cex.axis=3.5)
map_plot(lon_axis,lat_axis,mu6,zlim,title3,sc=3.5,sc.main=5,legend=TRUE,xl=0,yl=0,sub=NULL,xlim=xlim)
axis(2,at=c(-10,0,10),labels=c(-10,0,10),cex.lab=3.5,cex=3.5,cex.axis=3.5)
map_plot(lon_axis,lat_axis,s2,zlim0,title2,sc=3.5,sc.main=5,legend=TRUE,yl=0,xl=0,sub=NULL,xlim=xlim)
axis(2,at=c(-10,0,10),labels=c(-10,0,10),cex.lab=3.5,cex=3.5,cex.axis=3.5)
#
p <- 3
border = unlist(index_x_mv[[ii0[p]]][d1])
lon_axis = round(long[[ii0[p]]],1)          
lat_axis = round(lat[[ii0[p]]],1)
title1 = paste('N=',round(scalenx[ii0[p]],2))
title3 = paste('N=',round(scalenx[ii0[p]],2))
title2 = paste('N=',round(scalenx[ii0[p]],2))
map_plot(lon_axis,lat_axis,mu7,zlim,title1,sc=3.5,sc.main=5,legend=FALSE,yl=1,xl=0,sub=NULL,xlim=xlim2)
axis(2,at=c(-10,0,10),labels=c(-10,0,10),cex.lab=3.5,cex=3.5,cex.axis=3.5)
mtext('Longitude',line=6,side=1,cex=3.5)
map_plot(lon_axis,lat_axis,mu9,zlim,title3,sc=3.5,sc.main=5,legend=TRUE,xl=0,yl=0,sub=NULL,xlim=xlim2)
axis(2,at=c(-10,0,10),labels=c(-10,0,10),cex.lab=3.5,cex=3.5,cex.axis=3.5)
mtext('Longitude',line=6,side=1,cex=3.5)
map_plot(lon_axis,lat_axis,s3,zlim0,title2,sc=3.5,sc.main=5,legend=TRUE,yl=0,xl=0,sub=NULL,xlim=xlim2)
axis(2,at=c(-10,0,10),labels=c(-10,0,10),cex.lab=3.5,cex=3.5,cex.axis=3.5)
mtext('Longitude',line=6,side=1,cex=3.5)
add_legend(x=-1.1,y=-1,c('Precipitation (mm/day)','1e-3','1e-2','5e-2'),col=c(NA,rgb(240/255,240/255,240/255),rgb(130/255,130/255,130/255),rgb(15/255,15/255,15/255)),lty=c(NA,1,1,1),lwd=4,horiz=TRUE,xpd=TRUE, bty="n",cex=5)
dev.off() }

###---------------------------------- boxplot - param across scales - fig 5

## ---- graph for n=2 (main manuscript)
# single-scale estimates
load(file='optim_corr_spacetime_scale/par_ani1_phys_scalereg.RData')
load(file='optim_corr_spacetime_scale/par_ani2_phys_scalereg.RData')
load(file='optim_corr_spacetime_scale/par_temp_phys_scalereg.RData')
load(file='optim_corr_spacetime_scale/par_exp_phys_scalereg.RData')
load(file='optim_corr_spacetime_scale/par_var_phys_scalereg.RData')

# scale-aware estimates
load(file='par_estim_nx5710_phys_scale.RData')   # par1_estim,par2_estim,par3_estim,par4_estim,par5_estim
load('empty_window.RData')
P1 <- Map3Dto2D(par1_estim[-ind_na[,1],-ind_na[,2],])
P2 <- Map3Dto2D(par2_estim[-ind_na[,1],-ind_na[,2],])
P3 <- Map3Dto2D(par3_estim[-ind_na[,1],-ind_na[,2],])
P4 <- Map3Dto2D(par4_estim[-ind_na[,1],-ind_na[,2],])
P5 <- Map3Dto2D(par5_estim[-ind_na[,1],-ind_na[,2],])

sc <- 2.5
pdf('fboxplot_singlemultiscale_allres_theta123.pdf',height=10,width=15)
layout(matrix(1:6,2,3,byrow=TRUE))
par(mar=c(7,6,6,4),las=1)
fbplot(fit=t(P2),prob=0.5,ylim=c(0,5),xlim=c(1,11.25),color=rgb(.1,.6,0,.4),outliercol=0,factor=1.5,barcol=3,cex=sc,cex.axis=sc,cex.lab=sc,cex.main=sc,main=expression('Zonal anisotropy'~theta[Z]~'(degree'^-1~')'),xlab='',xaxt='n',ylab='')
boxplot(par_ani2,factor=1.5,xlab='',xaxs='n',xaxt='n',yaxt='n',pch=20,lwd=1.4,add=TRUE,col=rgb(0.6,0.1,0.1,alpha=0.6),medcol=2, whiskcol=2, staplecol=2, boxcol=2, outcol=2)
title(ylab=expression(theta[Z]),line=3,cex.lab=sc)
axis(1,at=1:n,label=round(scalenx,2),cex.axis=sc,cex.lab=sc,mgp=c(2,2,0))
#
fbplot(fit=t(P1),prob=0.5,ylim=c(0,5),xlim=c(1,11.25),color=rgb(.1,.6,0,.4),outliercol=0,factor=1.5,barcol=3,cex=sc,cex.axis=sc,cex.lab=sc,cex.main=sc,main=expression('Meridional anisotropy'~theta[M]~'(degree'^-1~')'),xlab='',xaxt='n',ylab='')
boxplot(par_ani1,factor=1.5,xlab='',xaxs='n',xaxt='n',yaxt='n',pch=20,lwd=1.4,add=TRUE,col=rgb(0.6,0.1,0.1,alpha=0.6),medcol=2, whiskcol=2, staplecol=2, boxcol=2, outcol=2)
title(ylab=expression(theta[M]),line=3,cex.lab=sc)
axis(1,at=1:n,label=round(scalenx,2),cex.axis=sc,cex.lab=sc,mgp=c(2,2,0))
#
fbplot(fit=t(P3),prob=0.5,ylim=c(0,15),xlim=c(1,11.25),color=rgb(.1,.6,0,.4),outliercol=0,factor=1.5,barcol=3,cex=sc,cex.axis=sc,cex.lab=sc,cex.main=sc,main=expression('Temporal range'~theta[T]~'(hour'^-1~')'),xlab='',xaxt='n',ylab='') 
boxplot(par_temp,factor=1.5,xlab='',xaxs='n',xaxt='n',yaxt='n',pch=20,lwd=1.4,add=TRUE,col=rgb(0.6,0.1,0.1,alpha=0.6),medcol=2, whiskcol=2, staplecol=2, boxcol=2, outcol=2)
title(ylab=expression(theta[T]),line=3,cex.lab=sc)
axis(1,at=1:n,label=round(scalenx,2),cex.axis=sc,cex.lab=sc,mgp=c(2,2,0))
#
fbplot(fit=t(P4),prob=0.5,ylim=c(0,2.5),xlim=c(1,11.25),color=rgb(.1,.6,0,.4),outliercol=0,factor=1.5,barcol=3,cex=sc,cex.axis=sc,cex.lab=sc,cex.main=sc,main=expression('Exponent'~gamma),xlab='',xaxt='n',ylab='')
boxplot(par_exp,factor=1.5,xlab='',xaxs='n',xaxt='n',yaxt='n',pch=20,lwd=1.4,add=TRUE,col=rgb(0.6,0.1,0.1,alpha=0.6),medcol=2, whiskcol=2, staplecol=2, boxcol=2, outcol=2)
title(xlab='Resolution N (degree)',line=5,cex.lab=sc)
title(ylab=expression(gamma),line=3,cex.lab=sc)
axis(1,at=1:n,label=round(scalenx,2),cex.axis=sc,cex.lab=sc,mgp=c(2,2,0))
#
fbplot(fit=t(P5),prob=0.5,ylim=c(0,2.5),xlim=c(1,11.25),color=rgb(.1,.6,0,.4),outliercol=0,factor=1.5,barcol=3,cex=sc,cex.axis=sc,cex.lab=sc,cex.main=sc,main=expression('Variance'~sigma),xlab='',xaxt='n',ylab='')
boxplot(par_var,factor=1.5,xlab='',xaxs='n',xaxt='n',yaxt='n',pch=20,lwd=1.4,add=TRUE,col=rgb(0.6,0.1,0.1,alpha=0.6),medcol=2, whiskcol=2, staplecol=2, boxcol=2, outcol=2)
title(xlab='Resolution N (degree)',line=5,cex.lab=sc)
title(ylab=expression(sigma),line=3,cex.lab=sc)
axis(1,at=1:n,label=round(scalenx,2),cex.axis=sc,cex.lab=sc,mgp=c(2,2,0))
#
plot(1,1,pch=NA,xaxt='n',yaxt='n',xlab='',ylab='',bty="n")
legend(x='center',c('Scale-aware Functional Boxplot','Single-scale Boxplot'),bty='n',col=c(3,2),pch=15,cex=2.2)
dev.off()

## ---- graph for n=1 (supplemental material)
# single-scale estimates
load(file='optim_corr_spacetime_scale/par_ani1_n1.RData')
load(file='optim_corr_spacetime_scale/par_ani2_n1.RData')
load(file='optim_corr_spacetime_scale/par_temp_n1.RData')
load(file='optim_corr_spacetime_scale/par_exp_n1.RData')
load(file='optim_corr_spacetime_scale/par_var_n1.RData')

# scale-aware estimates
load(file='par_estim_nx5710_n1.RData')   # par1_estim,par2_estim,par3_estim,par4_estim,par5_estim

load('empty_window.RData')
P1 <- Map3Dto2D(par1_estim[-ind_na[,1],-ind_na[,2],])
P2 <- Map3Dto2D(par2_estim[-ind_na[,1],-ind_na[,2],])
P3 <- Map3Dto2D(par3_estim[-ind_na[,1],-ind_na[,2],])
P4 <- Map3Dto2D(par4_estim[-ind_na[,1],-ind_na[,2],])
P5 <- Map3Dto2D(par5_estim[-ind_na[,1],-ind_na[,2],])

sc <- 2.5
pdf('fboxplot_singlemultiscale_allres_theta123_n1.pdf',height=10,width=15)
layout(matrix(1:6,2,3,byrow=TRUE))
par(mar=c(7,6,6,4),las=1)
fbplot(fit=t(P2),prob=0.5,ylim=c(0,5),xlim=c(1,11.25),color=rgb(.1,.6,0,.4),outliercol=0,factor=1.5,barcol=3,cex=sc,cex.axis=sc,cex.lab=sc,cex.main=sc,main=expression('Zonal anisotropy'~theta[Z]~'(degree'^-1~')'),xlab='',xaxt='n',ylab='')
boxplot(par_ani2,factor=1.5,xlab='',xaxs='n',xaxt='n',yaxt='n',pch=20,lwd=1.4,add=TRUE,col=rgb(0.6,0.1,0.1,alpha=0.6),medcol=2, whiskcol=2, staplecol=2, boxcol=2, outcol=2)
title(ylab=expression(theta[Z]),line=3,cex.lab=sc)
axis(1,at=1:n,label=round(scalenx,2),cex.axis=sc,cex.lab=sc,mgp=c(2,2,0))
#
fbplot(fit=t(P1),prob=0.5,ylim=c(0,5),xlim=c(1,11.25),color=rgb(.1,.6,0,.4),outliercol=0,factor=1.5,barcol=3,cex=sc,cex.axis=sc,cex.lab=sc,cex.main=sc,main=expression('Meridional anisotropy'~theta[M]~'(degree'^-1~')'),xlab='',xaxt='n',ylab='')
boxplot(par_ani1,factor=1.5,xlab='',xaxs='n',xaxt='n',yaxt='n',pch=20,lwd=1.4,add=TRUE,col=rgb(0.6,0.1,0.1,alpha=0.6),medcol=2, whiskcol=2, staplecol=2, boxcol=2, outcol=2)
title(ylab=expression(theta[M]),line=3,cex.lab=sc)
axis(1,at=1:n,label=round(scalenx,2),cex.axis=sc,cex.lab=sc,mgp=c(2,2,0))
#
fbplot(fit=t(P3),prob=0.5,ylim=c(0,15),xlim=c(1,11.25),color=rgb(.1,.6,0,.4),outliercol=0,factor=1.5,barcol=3,cex=sc,cex.axis=sc,cex.lab=sc,cex.main=sc,main=expression('Temporal range'~theta[T]~'(hour'^-1~')'),xlab='',xaxt='n',ylab='') 
boxplot(par_temp,factor=1.5,xlab='',xaxs='n',xaxt='n',yaxt='n',pch=20,lwd=1.4,add=TRUE,col=rgb(0.6,0.1,0.1,alpha=0.6),medcol=2, whiskcol=2, staplecol=2, boxcol=2, outcol=2)
title(ylab=expression(theta[T]),line=3,cex.lab=sc)
axis(1,at=1:n,label=round(scalenx,2),cex.axis=sc,cex.lab=sc,mgp=c(2,2,0))
#
fbplot(fit=t(P4),prob=0.5,ylim=c(0,2.5),xlim=c(1,11.25),color=rgb(.1,.6,0,.4),outliercol=0,factor=1.5,barcol=3,cex=sc,cex.axis=sc,cex.lab=sc,cex.main=sc,main=expression('Exponent'~gamma),xlab='',xaxt='n',ylab='')
boxplot(par_exp,factor=1.5,xlab='',xaxs='n',xaxt='n',yaxt='n',pch=20,lwd=1.4,add=TRUE,col=rgb(0.6,0.1,0.1,alpha=0.6),medcol=2, whiskcol=2, staplecol=2, boxcol=2, outcol=2)
title(xlab='Resolution N (degree)',line=5,cex.lab=sc)
title(ylab=expression(gamma),line=3,cex.lab=sc)
axis(1,at=1:n,label=round(scalenx,2),cex.axis=sc,cex.lab=sc,mgp=c(2,2,0))
#
fbplot(fit=t(P5),prob=0.5,ylim=c(0,2.5),xlim=c(1,11.25),color=rgb(.1,.6,0,.4),outliercol=0,factor=1.5,barcol=3,cex=sc,cex.axis=sc,cex.lab=sc,cex.main=sc,main=expression('Variance'~sigma),xlab='',xaxt='n',ylab='')
boxplot(par_var,factor=1.5,xlab='',xaxs='n',xaxt='n',yaxt='n',pch=20,lwd=1.4,add=TRUE,col=rgb(0.6,0.1,0.1,alpha=0.6),medcol=2, whiskcol=2, staplecol=2, boxcol=2, outcol=2)
title(xlab='Resolution N (degree)',line=5,cex.lab=sc)
title(ylab=expression(sigma),line=3,cex.lab=sc)
axis(1,at=1:n,label=round(scalenx,2),cex.axis=sc,cex.lab=sc,mgp=c(2,2,0))
#
plot(1,1,pch=NA,xaxt='n',yaxt='n',xlab='',ylab='',bty="n")
legend(x='center',c('Scale-aware Functional Boxplot','Single-scale Boxplot'),bty='n',col=c(3,2),pch=15,cex=2.2)
dev.off()


###---------------------------------- MLE uncertainty - fig 6

filename = 'optim_corr_spacetime_scale/singlescale/corrfit_residuals_nx14_n2_long14lat8_phys_scale_hessian.RData'
load(filename)
h1 = opt_corr$hessian
p1 = opt_corr$par
InvH1 = solve(h1)
S = svd(InvH1)
D = diag(S$d)
U = S$u
V = S$v
covP = U %*% (sqrt(D)) %*% t(V)
s1 = sqrt(S$d)
eps_par1 = array(0,c(length(p1),100))
for (i in 1:100){
  eps_par1[,i] = p1 + (covP) %*% rnorm(n=length(p1),mean=0,sd=1)}


filename = "optim_corr_spacetime_scale/nx5710/corrfit_residuals_multisc_nx5710_n2_long14lat8_phys_scale_hessian.RData"
load(filename)
h0 = optim_multiscale_mle$hessian
p0 = optim_multiscale_mle$par
InvH0 = solve(h0)
S0 = svd(InvH0)
D0 = diag(S0$d)
U0 = S0$u
V0 = S0$v
covP0 = U0 %*% (sqrt(D0)) %*% t(V0)

# multivariate delta method 
N = indnxn[6]
grad1 = c(exp(p0[2]*N),N*p0[1]*exp(p0[2]*N),rep(0,8))
covP1 = t(grad1) %*% covP0 %*% grad1
grad2 = c(0,0,exp(p0[4]*N),N*p0[3]*exp(p0[4]*N),rep(0,6))
covP2 = t(grad2) %*% covP0 %*% grad2
grad3 = c(rep(0,4),exp(p0[5]*N),N*p0[5]*exp(p0[6]*N),rep(0,4))
covP3 = t(grad3) %*% covP0 %*% grad3
grad4 = c(rep(0,6),(1-tanh(p0[7]+N*p0[8])^2),N*(1-tanh(p0[7]+N*p0[8])^2),rep(0,2))
covP4 = t(grad4) %*% covP0 %*% grad4
grad5 = c(rep(0,8),exp(p0[9]*N),N*p0[9]*exp(p0[10]*N))
covP5 = t(grad5) %*% covP0 %*% grad5
s0 = c(covP1,covP2,covP3,covP4,covP5)
theta0 = c(p0[1]*exp(p0[2]*N),p0[3]*exp(p0[4]*N),p0[5]*exp(p0[6]*N),(1-tanh(p0[7]+N*p0[8])),p0[9]*exp(p0[10]*N))
eps_par0 = array(0,c(5,100))
for (i in 1:100){
  eps_par0[,i] = theta0 + (s0) * rnorm(n=length(theta0),mean=0,sd=1)  }

sc = 3
setEPS()
postscript(file='estimation_uncertainty_nx6.eps',width=21,height=5)
par(oma = c(3,3,4,0),mai=c(0.4,0.4,0.4,0.4))
m = matrix(1:4,1,4,byrow=TRUE)
layout(m)
boxplot(t(eps_par0[1:2,]), boxfill = NA, border = NA,ylim=range(eps_par1[1:2,],eps_par0[1:2,],.8),xaxt='n',cex.axis=sc)  #invisible boxes - only axes and plot area
boxplot(t(eps_par1[1:2,]), xaxt = "n", add = TRUE, boxfill="red", boxwex=0.25, at = 1:ncol(t(eps_par1[1:2,])) - 0.15,cex.axis=sc) #shift these left by -0.15
boxplot(t(eps_par0[1:2,]), xaxt = "n", add = TRUE, boxfill="blue", boxwex=0.25, at = 1:ncol(t(eps_par0[1:2,])) + 0.15,cex.axis=sc) #shift to the right by +0.15
axis(1,at=1:2,c(expression(theta[Z]),expression(theta[M])),cex.axis=4,tick=FALSE,line=2)
#
boxplot((eps_par0[3,]), boxfill = NA, border = NA,ylim=range(eps_par1[3,],eps_par0[3,],5),xaxt='n',cex.axis=sc)  #invisible boxes - only axes and plot area
boxplot((eps_par1[3,]), xaxt = "n",add=TRUE,boxfill="red",boxwex=0.25,at=(1-0.15),cex.axis=sc) #shift these left by -0.15
boxplot((eps_par0[3,]), xaxt = "n",add=TRUE,boxfill="blue",boxwex=0.25,at=(1+0.15),cex.axis=sc) #shift to the right by +0.15
axis(1,at=1,c(expression(theta[T])),cex.axis=4,tick=FALSE,line=2)
#
boxplot((eps_par0[4,]), boxfill = NA, border = NA,ylim=range(eps_par1[4,],eps_par0[4,]),xaxt='n',cex.axis=sc)  #invisible boxes - only axes and plot area
boxplot((eps_par1[4,]), xaxt = "n",add=TRUE,boxfill="red",boxwex=0.25,at=(1-0.15),cex.axis=sc) #shift these left by -0.15
boxplot((eps_par0[4,]), xaxt = "n",add=TRUE,boxfill="blue",boxwex=0.25,at=(1+0.15),cex.axis=sc) #shift to the right by +0.15
axis(1,at=1,c(expression(gamma)),cex.axis=4,tick=FALSE,line=2)
#
boxplot((eps_par0[5,]), boxfill = NA, border = NA,ylim=range(eps_par1[5,],eps_par0[5,]),xaxt='n',cex.axis=sc)  #invisible boxes - only axes and plot area
boxplot((eps_par1[5,]), xaxt = "n",add=TRUE,boxfill="red",boxwex=0.25,at=(1-0.15),cex.axis=sc) #shift these left by -0.15
boxplot((eps_par0[5,]), xaxt = "n",add=TRUE,boxfill="blue",boxwex=0.25,at=(1+0.15),cex.axis=sc) #shift to the right by +0.15
axis(1,at=1,c(expression(sigma)),cex.axis=4,tick=FALSE,line=2)
add_legend("topleft", c('Single-scale','Scale-aware'),col=c('red','blue'),pch=15,cex=4,inset=0, xpd=TRUE, horiz=TRUE, bty="n")
add_legend("top", c('Estimation uncertainty'),col=NA,pch=NA,cex=4,inset=0, xpd=TRUE, horiz=TRUE, bty="n")
dev.off()
#




###---------------------------------- rho maps - fig 7
## ---- graphs for n=2 (main manuscript)
par11 = array(NA,c((d1-1),d2)) ; par12 = par11
par21 = array(NA,c((d1-1),d2)) ; par22 = par21
par31 = array(NA,c((d1-1),d2)) ; par32 = par31
par41 = array(NA,c((d1-1),d2)) ; par42 = par41
par51 = array(NA,c((d1-1),d2)) ; par52 = par51
for (i in 1:(d1-1)){
  for (j in 1:d2){
    filename <- paste('optim_corr_spacetime_scale/nx5710/corrfit_residuals_multisc_nx5710_n2_long',i,'lat',j,'.RData',sep='')
    l = tryCatch( load(filename),error=function(e) 'empty') 
    if (l != "empty"){
      if (length(optim_multiscale_mle)>1){
        par_opt <- optim_multiscale_mle$par
        par11[i,j] <- par_opt[1] ; par12[i,j] = par_opt[2]
        par21[i,j] <- par_opt[3] ; par22[i,j] = par_opt[4]
        par31[i,j] <- par_opt[5] ; par32[i,j] = par_opt[6]
        par41[i,j] <- par_opt[7] ; par42[i,j] <- par_opt[8]
        par51[i,j] <- par_opt[9] ; par52[i,j] <- par_opt[10]  } } } }

i <- 7
zp <- 86.4*apply(X=P[[i]],FUN=mean,MARGIN=1:2)
lon_axis <- round(long[[i]],1)          
lat_axis <- round(lat[[i]],1)
zlevels <- c(1e-3,1e-2,5e-2)

setEPS()
postscript(file='map_rho_nx5710_scale.eps',width=27,height=24)
par(mar=c(20,20,45,13),oma=c(2,1,7,5),mai=c(.8,1.2,.6,1))
layout(matrix(1:10,5,2,byrow=TRUE))
map_plot0(long_mv,lat_mv,par21,title=expression(theta[Z][','][1]),sc=3.5,sc.main=5,legend=TRUE,xl=0,yl=1,zlim=range(par21,par11,na.rm=TRUE))
contour(lon_axis,lat_axis,zp,col=c(rgb(240/255,240/255,240/255),rgb(130/255,130/255,130/255),rgb(15/255,15/255,15/255)),add=TRUE,labcex=1.8,lwd=2,level=zlevels,method='edge',drawlabels=FALSE)
map_plot0(long_mv,lat_mv,par22,title=expression(theta[Z][','][2]),sc=3.5,sc.main=5,legend=TRUE,xl=0,yl=0,zlim=range(par22,par12,na.rm=TRUE))
map_plot0(long_mv,lat_mv,par11,title=expression(theta[M][','][1]),sc=3.5,sc.main=5,legend=TRUE,xl=0,yl=1,zlim=range(par21,par11,na.rm=TRUE))
map_plot0(long_mv,lat_mv,par12,title=expression(theta[M][','][2]),sc=3.5,sc.main=5,legend=TRUE,xl=0,yl=0,zlim=range(par22,par12,na.rm=TRUE))
map_plot0(long_mv,lat_mv,par31,title=expression(theta[T][','][1]),sc=3.5,sc.main=5,legend=TRUE,xl=0,yl=1)
map_plot0(long_mv,lat_mv,par32,title=expression(theta[T][','][2]),sc=3.5,sc.main=5,legend=TRUE,xl=0,yl=0)
map_plot0(long_mv,lat_mv,par41,title=expression(gamma[1]),sc=3.5,sc.main=5,legend=TRUE,xl=0,yl=1)
map_plot0(long_mv,lat_mv,par42,title=expression(gamma[2]),sc=3.5,sc.main=5,legend=TRUE,xl=0,yl=0)
map_plot0(long_mv,lat_mv,par51,title=expression(sigma[1]),sc=3.5,sc.main=5,legend=TRUE,xl=1,yl=1)
map_plot0(long_mv,lat_mv,par52,title=expression(sigma[2]),sc=3.5,sc.main=5,legend=TRUE,xl=1,yl=0)
add_legend('topleft',c('Precipitation (mm/day)','1e-3','1e-2','5e-2'),col=c(NA,rgb(240/255,240/255,240/255),rgb(130/255,130/255,130/255),rgb(15/255,15/255,15/255)),lty=c(NA,1,1,1),lwd=4,horiz=TRUE,xpd=TRUE, bty="n",cex=4)
dev.off()
#


## ---- graphs for n=1 (main manuscript)
par11 = array(NA,c((d1-1),d2)) ; par12 = par11
par21 = array(NA,c((d1-1),d2)) ; par22 = par21
par31 = array(NA,c((d1-1),d2)) ; par32 = par31
par41 = array(NA,c((d1-1),d2)) ; par42 = par41
par51 = array(NA,c((d1-1),d2)) ; par52 = par51
for (i in 1:(d1-1)){
  for (j in 1:d2){
    filename <- paste('optim_corr_spacetime_scale/nx5710/corrfit_residuals_multisc_nx5710_n1_long',i,'lat',j,'.RData',sep='')
    l = tryCatch( load(filename),error=function(e) 'empty') 
    if (l != "empty"){
      if (length(optim_multiscale_mle)>1){
        par_opt <- optim_multiscale_mle$par
        par11[i,j] <- par_opt[1] ; par12[i,j] = par_opt[2]
        par21[i,j] <- par_opt[3] ; par22[i,j] = par_opt[4]
        par31[i,j] <- par_opt[5] ; par32[i,j] = par_opt[6]
        par41[i,j] <- par_opt[7] ; par42[i,j] <- par_opt[8]
        par51[i,j] <- par_opt[9] ; par52[i,j] <- par_opt[10]  } } } }

i <- 7
zp <- 86.4*apply(X=P[[i]],FUN=mean,MARGIN=1:2)
lon_axis <- round(long[[i]],1)          
lat_axis <- round(lat[[i]],1)
zlevels <- c(1e-3,1e-2,5e-2)

setEPS()
postscript(file='map_rho_nx5710_scale_n1.eps',width=27,height=24)
par(mar=c(20,20,45,13),oma=c(2,1,7,5),mai=c(.8,1.2,.6,1))
layout(matrix(1:10,5,2,byrow=TRUE))
map_plot0(long_mv,lat_mv,par21,title=expression(theta[Z][','][1]),sc=3.5,sc.main=5,legend=TRUE,xl=0,yl=1,zlim=range(par21,par11,na.rm=TRUE))
contour(lon_axis,lat_axis,zp,col=c(rgb(240/255,240/255,240/255),rgb(130/255,130/255,130/255),rgb(15/255,15/255,15/255)),add=TRUE,labcex=1.8,lwd=2,level=zlevels,method='edge',drawlabels=FALSE)
map_plot0(long_mv,lat_mv,par22,title=expression(theta[Z][','][2]),sc=3.5,sc.main=5,legend=TRUE,xl=0,yl=0,zlim=range(par22,par12,na.rm=TRUE))
map_plot0(long_mv,lat_mv,par11,title=expression(theta[M][','][1]),sc=3.5,sc.main=5,legend=TRUE,xl=0,yl=1,zlim=range(par21,par11,na.rm=TRUE))
map_plot0(long_mv,lat_mv,par12,title=expression(theta[M][','][2]),sc=3.5,sc.main=5,legend=TRUE,xl=0,yl=0,zlim=range(par22,par12,na.rm=TRUE))
map_plot0(long_mv,lat_mv,par31,title=expression(theta[T][','][1]),sc=3.5,sc.main=5,legend=TRUE,xl=0,yl=1)
map_plot0(long_mv,lat_mv,par32,title=expression(theta[T][','][2]),sc=3.5,sc.main=5,legend=TRUE,xl=0,yl=0)
map_plot0(long_mv,lat_mv,par41,title=expression(gamma[1]),sc=3.5,sc.main=5,legend=TRUE,xl=0,yl=1)
map_plot0(long_mv,lat_mv,par42,title=expression(gamma[2]),sc=3.5,sc.main=5,legend=TRUE,xl=0,yl=0)
map_plot0(long_mv,lat_mv,par51,title=expression(sigma[1]),sc=3.5,sc.main=5,legend=TRUE,xl=1,yl=1)
map_plot0(long_mv,lat_mv,par52,title=expression(sigma[2]),sc=3.5,sc.main=5,legend=TRUE,xl=1,yl=0)
add_legend('topleft',c('Precipitation (mm/day)','1e-3','1e-2','5e-2'),col=c(NA,rgb(240/255,240/255,240/255),rgb(130/255,130/255,130/255),rgb(15/255,15/255,15/255)),lty=c(NA,1,1,1),lwd=4,horiz=TRUE,xpd=TRUE, bty="n",cex=4)
dev.off()
#


###---------------------------------- global mse maps - fig 8
ii0 <- c(3,6,11) 

MSE = NULL
CMSE = NULL
BIAS = NULL
# multiscale GP 
load('mse_global_nx3_nx5710_phys_scale.RData')  #mse0, mse1
MSE[[1]] <- mse0
rm(mse0)
load('cmse_global_nx3_nx5710_phys_scale.RData') #cmse0,cmse1
CMSE[[1]] <- cmse0
rm(cmse0)
load('bias_global_nx3_nx5710_phys_scale.RData') #bias0, bias1
BIAS[[1]] <- bias0
rm(bias0)
####
load('mse_global_nx6_nx5710_phys_scale.RData')  #mse0, mse1
MSE[[2]] <- mse0
rm(mse0)
load('cmse_global_nx6_nx5710_phys_scale.RData') #cmse0,cmse1
CMSE[[2]] <- cmse0
rm(cmse0)
load('bias_global_nx6_nx5710_phys_scale.RData') #bias0, bias1
BIAS[[2]] <- bias0
rm(bias0)
####
load('mse_global_nx11_nx5710_phys_scale.RData')  #mse0, mse1
MSE[[3]] <- mse0
rm(mse0)
load('cmse_global_nx11_nx5710_phys_scale.RData') #cmse0,cmse1
CMSE[[3]] <- cmse0
rm(cmse0)
load('bias_global_nx11_nx5710_phys_scale.RData') #bias0, bias1
BIAS[[3]] <- bias0
rm(bias0)


# single-scale GP
####
load('mse_global_nx3_phys_scale.RData')  #mse0, mse1
MSE[[4]] <- mse0
rm(mse0)
load('cmse_global_nx3_phys_scale.RData') #cmse0,cmse1
CMSE[[4]] <- cmse0
rm(cmse0)
load('bias_global_nx3_phys_scale.RData') #bias0, bias1
BIAS[[4]] <- bias0
rm(bias0)
####
load('mse_global_nx6_phys_scale.RData')  #mse0, mse1
MSE[[5]] <- mse0
rm(mse0)
load('cmse_global_nx6_phys_scale.RData') #cmse0,cmse1
CMSE[[5]] <- cmse0
rm(cmse0)
load('bias_global_nx6_phys_scale.RData') #bias0, bias1
BIAS[[5]] <- bias0
rm(bias0)
####
load('mse_global_nx11_phys_scale.RData')  #mse0, mse1
MSE[[6]] <- mse0
rm(mse0)
load('cmse_global_nx11_phys_scale.RData') #cmse0,cmse1
CMSE[[6]] <- cmse0
rm(cmse0)
load('bias_global_nx11_phys_scale.RData') #bias0, bias1
BIAS[[6]] <- bias0
rm(bias0)

round(median(100*(MSE[[1]]-MSE[[4]])/MSE[[4]],na.rm=TRUE),2)
round(median(100*(MSE[[2]]-MSE[[5]])/MSE[[5]],na.rm=TRUE),2)
round(median(100*(MSE[[3]]-MSE[[6]])/MSE[[6]],na.rm=TRUE),2)

i <- 7
zp <- 86.4*apply(X=P[[i]],FUN=mean,MARGIN=1:2)
lon_prcp <- round(long[[i]],1)          
lat_prcp <- round(lat[[i]],1)
zlevels <- c(1e-3,1e-2,5e-2)

zlim <- log(range(unlist(MSE),unlist(CMSE),unlist(BIAS),na.rm=TRUE))
z0 <- log(c(unlist(MSE),unlist(CMSE),unlist(BIAS)))
rbPal <- parula
nc <- 200 
brks <- seq(min(z0,na.rm=TRUE),max(z0,na.rm=TRUE),length.out=(nc+1))
Col <- rbPal(nc)

xlim <- c(42,175)
fig8 <- 1
#par(mar=c(35,10,20,11),oma=c(9,.5,5,3),mai=c(.8,1.8,.6,1.2))
if (fig8 == 1){
setEPS()      
postscript('maps_mse_scaleaware_nx5710.eps',height=18,width=35)
par(mfrow=c(4,3),mar=c(35,15,20,11),oma=c(9,4,2.5,3),mai=c(.7,1,.7,1.4))
p <- 1
border = unlist(index_x_mv[[ii0[p]]][d1])
lon_axis = round(long[[ii0[p]]],1)[-border]          
lat_axis = round(lat[[ii0[p]]],1)
map_plot(lon_axis,lat_axis,log(MSE[[1]][-border,]),zlim,title=paste('Scale-aware - N=',round(scalenx[ii0[p]],2)),sub=NULL,sc=4,sc.main=4.5,legend=FALSE,xl=0,yl=1,xlim,col=Col,breaks=brks)
axis(2,at=c(-10,0,10),labels=c(-10,0,10),cex.lab=4,cex=4,cex.axis=4)
contour(lon_prcp,lat_prcp,zp,col=c(rgb(240/255,240/255,240/255),rgb(130/255,130/255,130/255),rgb(15/255,15/255,15/255)),add=TRUE,labcex=1.8,lwd=2,level=zlevels,method='edge',drawlabels=FALSE)
map_plot(lon_axis,lat_axis,log(CMSE[[1]])[-border,],zlim,title=paste('Scale-aware - N=',round(scalenx[ii0[p]],2)),sub=NULL,sc=4,sc.main=4.5,legend=FALSE,xl=0,yl=0,xlim,col=Col,breaks=brks)
axis(2,at=c(-10,0,10),labels=c(-10,0,10),cex.lab=4,cex=4,cex.axis=4)
map_plot(lon_axis,lat_axis,log(BIAS[[1]][-border,]),zlim,title=paste('Scale-aware - N=',round(scalenx[ii0[p]],2)),sub=NULL,sc=4,sc.main=4.5,legend=TRUE,xl=0,yl=0,xlim,col=Col,breaks=brks)
axis(2,at=c(-10,0,10),labels=c(-10,0,10),cex.lab=4,cex=4,cex.axis=4)
##
p <- 2
border = unlist(index_x_mv[[ii0[p]]][d1])
lon_axis = round(long[[ii0[p]]],1)[-border]          
lat_axis = round(lat[[ii0[p]]],1)
map_plot(lon_axis,lat_axis,log(MSE[[2]][-border,]),zlim,title=paste('Scale-aware - N=',round(scalenx[ii0[p]],2)),sub=NULL,sc=4,sc.main=4.5,legend=FALSE,xl=0,yl=1,xlim,col=Col,breaks=brks)
axis(2,at=c(-10,0,10),labels=c(-10,0,10),cex.lab=4,cex=4,cex.axis=4)
map_plot(lon_axis,lat_axis,log(CMSE[[2]][-border,]),zlim,title=paste('Scale-aware - N=',round(scalenx[ii0[p]],2)),sub=NULL,sc=4,sc.main=4.5,legend=FALSE,xl=0,yl=0,xlim,col=Col,breaks=brks)
axis(2,at=c(-10,0,10),labels=c(-10,0,10),cex.lab=4,cex=4,cex.axis=4)
map_plot(lon_axis,lat_axis,log(BIAS[[2]][-border,]),zlim,title=paste('Scale-aware - N=',round(scalenx[ii0[p]],2)),sub=NULL,sc=4,sc.main=4.5,legend=TRUE,xl=0,yl=0,xlim,col=Col,breaks=brks)
axis(2,at=c(-10,0,10),labels=c(-10,0,10),cex.lab=4,cex=4,cex.axis=4)
##
p <- 3
lon_axis = round(long[[ii0[p]]],1)          
lat_axis = round(lat[[ii0[p]]],1)
map_plot(lon_axis,lat_axis,log(MSE[[3]]),zlim,title=paste('Scale-aware - N=',round(scalenx[ii0[p]],2)),sub=NULL,sc=4,sc.main=4.5,legend=FALSE,xl=0,yl=1,xlim,col=Col,breaks=brks)
axis(2,at=c(-10,0,10),labels=c(-10,0,10),cex.lab=4,cex=4,cex.axis=4)
map_plot(lon_axis,lat_axis,log(CMSE[[3]]),zlim,title=paste('Scale-aware - N=',round(scalenx[ii0[p]],2)),sub=NULL,sc=4,sc.main=4.5,legend=FALSE,xl=0,yl=0,xlim,col=Col,breaks=brks)
axis(2,at=c(-10,0,10),labels=c(-10,0,10),cex.lab=4,cex=4,cex.axis=4)
map_plot(lon_axis,lat_axis,log(BIAS[[3]]),zlim,title=paste('Scale-aware - N=',round(scalenx[ii0[p]],2)),sub=NULL,sc=4,sc.main=4.5,legend=TRUE,xl=0,yl=0,xlim,col=Col,breaks=brks)
axis(2,at=c(-10,0,10),labels=c(-10,0,10),cex.lab=4,cex=4,cex.axis=4)
##
p <- 2
border = unlist(index_x_mv[[ii0[p]]][d1])
lon_axis = round(long[[ii0[p]]],1)[-border]          
lat_axis = round(lat[[ii0[p]]],1)
map_plot(lon_axis,lat_axis,log(MSE[[5]][-border,]),zlim,title=paste('Scale-specific - N=',round(scalenx[ii0[p]],2)),sub=NULL,sc=4,sc.main=4.5,legend=FALSE,xl=0,yl=1,xlim,col=Col,breaks=brks)
axis(2,at=c(-10,0,10),labels=c(-10,0,10),cex.lab=4,cex=4,cex.axis=4)
mtext('Longitude',line=6,side=1,cex=3.5)
map_plot(lon_axis,lat_axis,log(CMSE[[5]][-border,]),zlim,title=paste('Scale-specific - N=',round(scalenx[ii0[p]],2)),sub=NULL,sc=4,sc.main=4.5,legend=FALSE,xl=0,yl=0,xlim,col=Col,breaks=brks)
axis(2,at=c(-10,0,10),labels=c(-10,0,10),cex.lab=4,cex=4,cex.axis=4)
mtext('Longitude',line=6,side=1,cex=3.5)
map_plot(lon_axis,lat_axis,log(BIAS[[5]][-border,]),zlim,title=paste('Scale-specific - N=',round(scalenx[ii0[p]],2)),sub=NULL,sc=4,sc.main=4.5,legend=TRUE,xl=0,yl=0,xlim,col=Col,breaks=brks)
axis(2,at=c(-10,0,10),labels=c(-10,0,10),cex.lab=4,cex=4,cex.axis=4)
mtext('Longitude',line=6,side=1,cex=3.5)
add_legend(x=-1.1,y=-1.03,c('Precipitation (mm/day)','1e-3','1e-2','5e-2'),col=c(NA,rgb(240/255,240/255,240/255),rgb(130/255,130/255,130/255),rgb(15/255,15/255,15/255)),lty=c(NA,1,1,1),lwd=4,horiz=TRUE,xpd=TRUE, bty="n",cex=5)
dev.off() }




## graph - slides
zlim = range(MSE[[5]],CMSE[[5]],BIAS[[5]],MSE[[2]],CMSE[[2]],BIAS[[2]],na.rm=TRUE)
p <- 2
border = unlist(index_x_mv[[ii0[p]]][d1])
lon_axis = round(long[[ii0[p]]],1)[-border]          
lat_axis = round(lat[[ii0[p]]],1)
#
setEPS()      
postscript('maps_mse_scaleaware_nx5710_k2_singlescale.eps',height=2.4,width=15)
par(mfrow=c(1,3),mar=c(5, 4, 4, 14),oma=c(3,3,3,2.5),mai=c(.2,.3,.2,.5))
map_plot(lon_axis,lat_axis,MSE[[5]][-border,],zlim,title=NULL,sc=2,sc.main=3,legend=FALSE,xl=0,yl=0)
map_plot(lon_axis,lat_axis,CMSE[[5]][-border,],zlim,title=NULL,sc=2,sc.main=3,legend=FALSE,xl=0,yl=0)
map_plot(lon_axis,lat_axis,BIAS[[5]][-border,],zlim,title=NULL,sc=2,sc.main=3,legend=TRUE,xl=0,yl=0)
dev.off()
#
setEPS()      
postscript('maps_mse_scaleaware_nx5710_k2_scaleaware.eps',height=2.4,width=15)
par(mfrow=c(1,3),mar=c(5, 4, 4, 14),oma=c(3,3,3,2.5),mai=c(.2,.3,.2,.5))
map_plot(lon_axis,lat_axis,MSE[[2]][-border,],zlim,title=NULL,sc=2,sc.main=3,legend=FALSE,xl=0,yl=0)
map_plot(lon_axis,lat_axis,CMSE[[2]][-border,],zlim,title=NULL,sc=2,sc.main=3,legend=FALSE,xl=0,yl=0)
map_plot(lon_axis,lat_axis,BIAS[[2]][-border,],zlim,title=NULL,sc=2,sc.main=3,legend=TRUE,xl=0,yl=0)
dev.off()

zlim <- range(log(MSE[[5]]),log(MSE[[2]]),na.rm=TRUE)
setEPS()      
postscript('maps_mse_scaleaware_k2_singlescale.eps',height=1.58,width=10)
par(mfrow=c(1,2),mar=c(5, 4, 4, 16),oma=c(.5,.5,.5,7),mai=c(.2,.3,.2,.7))
map_plot0(lon_axis,lat_axis,log(MSE[[5]][-border,]),zlim,title=NULL,sc=1.5,sc.main=3,legend=FALSE,xl=0,yl=0)
map_plot0(lon_axis,lat_axis,log(MSE[[2]][-border,]),zlim,title=NULL,sc=1.5,sc.main=3,legend=TRUE,xl=0,yl=0)
dev.off()
#
setEPS()      
postscript('maps_mse_scaleaware_k2_ratiomse.eps',height=2.1,width=8)
par(mar=c(1,1,1,5),oma=c(.5,.5,.5,4))
map_plot0(lon_axis,lat_axis,100*(MSE[[5]][-border,]- MSE[[2]][-border,])/MSE[[2]][-border,],title=NULL,sc=1.5,sc.main=3,legend=TRUE,xl=0,yl=0)
dev.off()


###---------------------------------- time series - fig 9

ii0 <- c(3,6,11) 
load(file='regression_scaleaware_nx5710.RData')
coeff = scale.lm$coefficients
thetaN = matrix(coeff,8,3,byrow=TRUE)
parN = thetaN %*% t(cbind(rep(1,n),indnx,indnx^2))
parN[1,] = coeff[1] + coeff[2]*log(indnx) + coeff[3]*indnx*indnx

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(1, 0, 0, 0),mar=c(1, 0, 0, 0),new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)  }

timets <- c('Apr. 8','Apr. 9','Apr. 10','Apr. 11','Apr. 12','Apr. 13','Apr. 14','Apr. 15')

ylim <- c(-5.7,4.5)
sc <- 3.5
graph = 1
if (graph ==1){
setEPS()      
postscript('timeseries_scaleaware_nx5710_nx3611.eps',height=20,width=35)
layout(matrix(1:9,3,3,byrow=TRUE))
par(mar=c(24,21,6,10),oma=c(6,5,1,1),mai=c(.7,.7,.7,.7))
p <- 1
load(file='selectedpoints_valid_nx3.RData') # min.mse, max.mse, med.mse
load('Ysim_nx3singlescale_phys_scalereg.RData') #Ysim1
load('Ysim_nx3_multiscale_phys_scalereg.RData') #Ysim0
load('mse_global_nx3_nx5710_phys_scale.RData')
X1 <- log(Ur[[ii0[p]]])
X2 <- P[[ii0[p]]]
N <- indnx[ii0[p]]
N2 <- N*N
Xp0 <- parN[2,ii0[p]]*X1 + parN[3,ii0[p]]*X1^2 + parN[4,ii0[p]]*X1^3
Xp1 <- parN[5,ii0[p]]*X2 + parN[6,ii0[p]]*X2^.75 + parN[7,ii0[p]]*X2^.5 + parN[8,ii0[p]]*X2^.25
Xp <- parN[1,ii0[p]] + Xp0 + Xp1
err0 <- Ut[[ii0[p]]] - Ur[[ii0[p]]]
err0[err0==0] <- NA 
err <- log(err0)
matplot(Y1sim[q25.mse[1],q25.mse[2],,1:10],typ='l',lwd=2,lty=1,col=16,ylab='',xlab='',main=paste('25th-quantile MSE (',round(mse0[q25.mse[1],q25.mse[2]],2),') - N=',round(scalenx[ii0[p]],2),sep=''),xaxs='i',ylim=ylim,cex.main=4.5,cex=sc,xaxt='n',yaxt='n')
matplot(Y0sim[q25.mse[1],q25.mse[2],,1:10],add=TRUE,typ='l',lwd=2,lty=1,col='lightblue')
points(Xp[q25.mse[1],q25.mse[2],],col=1,pch=20,cex=1)
lines(err[q25.mse[1],q25.mse[2],],col=1,lwd=2)
axis(1,at=seq(21,192,by=24),label=timets,srt=45,cex=sc,cex.axis=4,mgp=c(3.5,3,0))
mtext(expression(epsilon[2]~'(N)'),side=2,line=5,cex=4)
axis(2,at=c(-5,-2.5,0,2.5),c(-5,-2.5,0,2.5),cex.axis=4)
#
matplot(Y1sim[med.mse[1],med.mse[2],,1:10],lwd=2,typ='l',lty=1,col=16,ylab='',xlab='',main=paste('Median MSE (',round(mse0[med.mse[1],med.mse[2]],2),') - N=',round(scalenx[ii0[p]],2),sep=''),xaxs='i',ylim=ylim,cex.main=4.5,cex=sc,xaxt='n',yaxt='n')
matplot(Y0sim[med.mse[1],med.mse[2],,1:10],lwd=2,add=TRUE,typ='l',lty=1,col='lightblue')
points(Xp[med.mse[1],med.mse[2],],col=1,pch=20,cex=1)
lines(err[med.mse[1],med.mse[2],],col=1,lwd=2)
axis(1,at=seq(21,192,by=24),label=timets,srt=45,cex=sc,cex.axis=4,mgp=c(3.5,3,0))
axis(2,at=c(-5,-2.5,0,2.5),c(-5,-2.5,0,2.5),cex.axis=4)
#
matplot(Y1sim[q75.mse[1],q75.mse[2],,1:10],lwd=2,typ='l',lty=1,col=16,ylab='',xlab='',main=paste('75th-quantile MSE (',round(mse0[q75.mse[1],q75.mse[2]],2),') - N=',round(scalenx[ii0[p]],2),sep=''),xaxs='i',ylim=ylim,cex.main=4.5,cex=sc,xaxt='n',yaxt='n')
matplot(Y0sim[q75.mse[1],q75.mse[2],,1:10],lwd=2,add=TRUE,typ='l',lty=1,col='lightblue')
points(Xp[q75.mse[1],q75.mse[2],],col=1,pch=20,cex=1)
lines(err[q75.mse[1],q75.mse[2],],col=1,lwd=2)
axis(1,at=seq(21,192,by=24),label=timets,srt=45,cex=sc,cex.axis=4,mgp=c(3.5,3,0))
axis(2,at=c(-5,-2.5,0,2.5),c(-5,-2.5,0,2.5),cex.axis=4)
####
p <- 2
load(file='selectedpoints_valid_nx6.RData') # min.mse, max.mse, med.mse
load('Ysim_nx6singlescale_phys_scalereg.RData') #Ysim1
load('Ysim_nx6_multiscale_phys_scalereg.RData') #Ysim0
load('mse_global_nx6_nx5710_phys_scale.RData')
X1 <- log(Ur[[ii0[p]]])
X2 <- P[[ii0[p]]]
N <- indnx[ii0[p]]
N2 <- N*N
Xp0 <- parN[2,ii0[p]]*X1 + parN[3,ii0[p]]*X1^2 + parN[4,ii0[p]]*X1^3
Xp1 <- parN[5,ii0[p]]*X2 + parN[6,ii0[p]]*X2^.75 + parN[7,ii0[p]]*X2^.5 + parN[8,ii0[p]]*X2^.25
Xp <- parN[1,ii0[p]] + Xp0 + Xp1
err0 <- Ut[[ii0[p]]] - Ur[[ii0[p]]]
err0[err0==0] <- NA 
err <- log(err0)
matplot(Y1sim[q25.mse[1],q25.mse[2],,1:10],typ='l',lwd=2,lty=1,col=16,ylab='',xlab='',main=paste('25th-quantile MSE (',round(mse0[q25.mse[1],q25.mse[2]],2),') - N=',round(scalenx[ii0[p]],2),sep=''),xaxs='i',ylim=ylim,cex.main=4.5,cex=sc,xaxt='n',yaxt='n')
matplot(Y0sim[q25.mse[1],q25.mse[2],,1:10],add=TRUE,typ='l',lwd=2,lty=1,col='lightblue')
points(Xp[q25.mse[1],q25.mse[2],],col=1,pch=20,cex=1)
lines(err[q25.mse[1],q25.mse[2],],col=1,lwd=2)
axis(1,at=seq(21,192,by=24),label=timets,srt=45,cex=sc,cex.axis=4,mgp=c(3.5,3,0))
mtext(expression(epsilon[2]~'(N)'),side=2,line=5,cex=4)
axis(2,at=c(-5,-2.5,0,2.5),c(-5,-2.5,0,2.5),cex.axis=4)
#
matplot(Y1sim[med.mse[1],med.mse[2],,1:10],lwd=2,typ='l',lty=1,col=16,ylab='',xlab='',main=paste('Median MSE (',round(mse0[med.mse[1],med.mse[2]],2),') - N=',round(scalenx[ii0[p]],2),sep=''),xaxs='i',ylim=ylim,cex.main=4.5,cex=sc,xaxt='n',yaxt='n')
matplot(Y0sim[med.mse[1],med.mse[2],,1:10],lwd=2,add=TRUE,typ='l',lty=1,col='lightblue')
points(Xp[med.mse[1],med.mse[2],],col=1,pch=20,cex=1)
lines(err[med.mse[1],med.mse[2],],col=1,lwd=2)
axis(1,at=seq(21,192,by=24),label=timets,srt=45,cex=sc,cex.axis=4,mgp=c(3.5,3,0))
axis(2,at=c(-5,-2.5,0,2.5),c(-5,-2.5,0,2.5),cex.axis=4)
#
matplot(Y1sim[q75.mse[1],q75.mse[2],,1:10],lwd=2,typ='l',lty=1,col=16,ylab='',xlab='',main=paste('75th-quantile MSE (',round(mse0[q75.mse[1],q75.mse[2]],2),') - N=',round(scalenx[ii0[p]],2),sep=''),xaxs='i',ylim=ylim,cex.main=4.5,cex=sc,xaxt='n',yaxt='n')
matplot(Y0sim[q75.mse[1],q75.mse[2],,1:10],lwd=2,add=TRUE,typ='l',lty=1,col='lightblue')
points(Xp[q75.mse[1],q75.mse[2],],col=1,pch=20,cex=1)
lines(err[q75.mse[1],q75.mse[2],],col=1,lwd=2)
axis(1,at=seq(21,192,by=24),label=timets,srt=45,cex=sc,cex.axis=4,mgp=c(3.5,3,0))
axis(2,at=c(-5,-2.5,0,2.5),c(-5,-2.5,0,2.5),cex.axis=4)
####
p<-3
load(file='selectedpoints_valid_nx11.RData') # min.mse, max.mse, med.mse
load('Ysim_nx11singlescale_phys_scalereg.RData') #Ysim1
load('Ysim_nx11_multiscale_phys_scalereg.RData') #Ysim0
load('mse_global_nx11_nx5710_phys_scale.RData')
X1 <- log(Ur[[ii0[p]]])
X2 <- P[[ii0[p]]]
N <- indnx[ii0[p]]
N2 <- N*N
Xp0 <- parN[2,ii0[p]]*X1 + parN[3,ii0[p]]*X1^2 + parN[4,ii0[p]]*X1^3
Xp1 <- parN[5,ii0[p]]*X2 + parN[6,ii0[p]]*X2^.75 + parN[7,ii0[p]]*X2^.5 + parN[8,ii0[p]]*X2^.25
Xp <- parN[1,ii0[p]] + Xp0 + Xp1
err0 <- Ut[[ii0[p]]] - Ur[[ii0[p]]]
err0[err0==0] <- NA 
err <- log(err0)
matplot(Y1sim[q25.mse[1],q25.mse[2],,1:10],typ='l',lwd=2,lty=1,col=16,ylab='',xlab='',main=paste('25th-quantile MSE (',round(mse0[q25.mse[1],q25.mse[2]],2),') - N=',round(scalenx[ii0[p]],2),sep=''),xaxs='i',ylim=ylim,cex.main=4.5,cex=sc,xaxt='n',yaxt='n')
matplot(Y0sim[q25.mse[1],q25.mse[2],,1:10],add=TRUE,typ='l',lwd=2,lty=1,col='lightblue')
points(Xp[q25.mse[1],q25.mse[2],],col=1,pch=20,cex=1)
lines(err[q25.mse[1],q25.mse[2],],col=1,lwd=2)
axis(1,at=seq(21,192,by=24),label=timets,srt=45,cex=sc,cex.axis=4,mgp=c(3.5,3,0))
mtext('Date (dt=hour)',side=1,line=8,cex=sc)
axis(2,at=c(-5,-2.5,0,2.5),c(-5,-2.5,0,2.5),cex.axis=4)
mtext(expression(epsilon[2]~'(N)'),side=2,line=5,cex=4)
legend(x=1,y=-.3,c('True','Scale-aware','Single-scale','Predicted mean'),col=c(1,'lightblue',16,1),lty=c(1,1,1,NA),pch=c(NA,NA,NA,20),lwd=c(3,3,3,NA),cex=4.5,bty='n')
#
matplot(Y1sim[med.mse[1],med.mse[2],,1:10],lwd=2,typ='l',lty=1,col=16,ylab='',xlab='',main=paste('Median MSE (',round(mse0[med.mse[1],med.mse[2]],2),') - N=',round(scalenx[ii0[p]],2),sep=''),xaxs='i',ylim=ylim,cex.main=4.5,cex=sc,xaxt='n',yaxt='n')
matplot(Y0sim[med.mse[1],med.mse[2],,1:10],lwd=2,add=TRUE,typ='l',lty=1,col='lightblue')
points(Xp[med.mse[1],med.mse[2],],col=1,pch=20,cex=1)
lines(err[med.mse[1],med.mse[2],],col=1,lwd=2)
axis(1,at=seq(21,192,by=24),label=timets,srt=45,cex=sc,cex.axis=4,mgp=c(3.5,3,0))
mtext('Date (dt=hour)',side=1,line=8,cex=sc)
axis(2,at=c(-5,-2.5,0,2.5),c(-5,-2.5,0,2.5),cex.axis=4)
#
matplot(Y1sim[q75.mse[1],q75.mse[2],,1:10],lwd=2,typ='l',lty=1,col=16,ylab='',xlab='',main=paste('75th-quantile MSE (',round(mse0[q75.mse[1],q75.mse[2]],2),') - N=',round(scalenx[ii0[p]],2),sep=''),xaxs='i',ylim=ylim,cex.main=4.5,cex=sc,xaxt='n',yaxt='n')
matplot(Y0sim[q75.mse[1],q75.mse[2],,1:10],lwd=2,add=TRUE,typ='l',lty=1,col='lightblue')
points(Xp[q75.mse[1],q75.mse[2],],col=1,pch=20,cex=1)
lines(err[q75.mse[1],q75.mse[2],],col=1,lwd=2)
axis(1,at=seq(21,192,by=24),label=timets,srt=45,cex=sc,cex.axis=4,mgp=c(3.5,3,0))
mtext('Date (dt=hour)',side=1,line=8,cex=sc)
axis(2,at=c(-5,-2.5,0,2.5),c(-5,-2.5,0,2.5),cex.axis=4)
dev.off()
}




###---------------------------------- ACF - fig14

ii0 <- c(3,6,11) 

p <- 3
load(file='selectedpoints_valid_nx11.RData') # min.mse, max.mse, med.mse
load('Ysim_nx11singlescale_phys_scalereg.RData') #Ysim1
load('Ysim_nx11_multiscale_phys_scalereg.RData') #Ysim0
S1 <- Y1sim[med.mse[1],med.mse[2],,]
M1 <- Y0sim[med.mse[1],med.mse[2],,]
rm(Y0sim,Y1sim)
err0 <- Ut[[ii0[p]]] - Ur[[ii0[p]]]
err0[err0==0] <- NA 
err1 <- log(err0)[med.mse[1],med.mse[2],]

p <- 2
load(file='selectedpoints_valid_nx6.RData') # min.mse, max.mse, med.mse
load('Ysim_nx6singlescale_phys_scalereg.RData') #Ysim1
load('Ysim_nx6_multiscale_phys_scalereg.RData') #Ysim0
S2 <- Y1sim[med.mse[1],med.mse[2],,]
M2 <- Y0sim[med.mse[1],med.mse[2],,]
rm(Y0sim,Y1sim)
err0 <- Ut[[ii0[p]]] - Ur[[ii0[p]]]
err0[err0==0] <- NA 
err2 <- log(err0)[med.mse[1],med.mse[2],]

p <- 1
load(file='selectedpoints_valid_nx3.RData') # min.mse, max.mse, med.mse
load('Ysim_nx3singlescale_phys_scalereg.RData') #Ysim1
load('Ysim_nx3_multiscale_phys_scalereg.RData') #Ysim0
S3 <- Y1sim[med.mse[1],med.mse[2],,]
M3 <- Y0sim[med.mse[1],med.mse[2],,]
rm(Y0sim,Y1sim)
err0 <- Ut[[ii0[p]]] - Ur[[ii0[p]]]
err0[err0==0] <- NA 
err3 <- log(err0)[med.mse[1],med.mse[2],]


Nsim <- dim(S1)[2]
acf_sim0 = array(0,c(3,Nsim,23))
acf_sim1 = array(0,c(3,Nsim,23))
for (s in 1:Nsim){
  acf_sim0[1,s,] = acf(M1[,s],plot=FALSE)$acf
  acf_sim0[2,s,] = acf(M2[,s],plot=FALSE)$acf
  acf_sim0[3,s,] = acf(M3[,s],plot=FALSE)$acf
  acf_sim1[1,s,] = acf(S1[,s],plot=FALSE)$acf
  acf_sim1[2,s,] = acf(S2[,s],plot=FALSE)$acf
  acf_sim1[3,s,] = acf(S3[,s],plot=FALSE)$acf }

a1 = acf(err1,plot=FALSE)$acf
a2 = acf(err2,plot=FALSE)$acf
a3 = acf(err3,plot=FALSE)$acf


sc <- 2.6
setEPS()      
postscript('acf_scaleaware_nx5710.eps',height=5,width=15)
par(mfrow=c(1,3),mar=c(14,17,6,15),oma=c(3,4,1,5),mai=c(.5,.5,.5,.5))
matplot(t(acf_sim1[3,,]),typ='l',lty=1,col=16,xlab='',ylab='',main='N=0.25',ylim=range(acf_sim0,acf_sim1),lwd=1.5,xaxs='i',xaxt='n',cex.axis=sc,cex.main=3.5)
matplot(t(acf_sim0[3,,]),add=TRUE,typ='l',lty=1,col='lightblue',lwd=1.5)
lines(a3,lwd=1.5)
axis(1,at=c(5,10,15,20),c(5,10,15,20),cex.axis=sc)
mtext('Time-lag (hr)',side=1,line=4,cex=2.2)
mtext('Temporal ACF',side=2,line=4,cex=2.2)
legend(x='topright',c('Observed','Scale-aware','Single-scale'),col=c(1,'lightblue',16),lty=1,cex=2,bty='n',lwd=2)
abline(h=0)
abline(h=c(1.96/sqrt(214),-1.96/sqrt(214)),lty=2)
#
matplot(t(acf_sim1[2,,]),typ='l',lty=1,col=16,xlab='',ylab='',main='N=0.5',ylim=range(acf_sim0,acf_sim1),lwd=1.5,xaxs='i',xaxt='n',cex.axis=sc,cex.main=3.5)
matplot(t(acf_sim0[2,,]),add=TRUE,typ='l',lty=1,col='lightblue',lwd=1.5)
lines(a2,lwd=1.5)
axis(1,at=c(5,10,15,20),c(5,10,15,20),cex.axis=sc)
mtext('Time-lag (hr)',side=1,line=4,cex=2.2)
legend(x='topright',c('Observed','Scale-aware','Single-scale'),col=c(1,'lightblue',16),lty=1,cex=2,bty='n',lwd=2)
abline(h=0)
abline(h=c(1.96/sqrt(214),-1.96/sqrt(214)),lty=2)
#
matplot(t(acf_sim1[1,,]),typ='l',lty=1,col=16,xlab='',ylab='',main='N=1.76',ylim=range(acf_sim0,acf_sim1),lwd=1.5,xaxs='i',xaxt='n',cex.axis=sc,cex.main=3.5)
matplot(t(acf_sim0[1,,]),add=TRUE,typ='l',lty=1,col='lightblue',lwd=1.5)
lines(a1,lwd=1.5)
axis(1,at=c(5,10,15,20),c(5,10,15,20),cex.axis=sc)
mtext('Time-lag (hr)',side=1,line=4,cex=2.2)
legend(x='topright',c('Observed','Scale-aware','Single-scale'),col=c(1,'lightblue',16),lty=1,cex=2,bty='n',lwd=2)
abline(h=0)
abline(h=c(1.96/sqrt(214),-1.96/sqrt(214)),lty=2)
dev.off()




###---------------------------------- covariance maps - fig 13
## maps of covariance coutours 
load('Ysim_nx6singlescale_phys_scalereg.RData') #Ysim1
load('Ysim_nx6_multiscale_phys_scalereg.RData') #Ysim0

graph = 1
sc = 2.5
if (graph==1){
ii0 <- c(3,6,11) 
p <- 2
setEPS()      
postscript('map_contour_spatialcorr_scaleaware_nx5710_nx6.eps',height=14,width=13)
layout(matrix(1:3,3,1))
par(mar=c(5,5,3,2))
plot(NA,xlim=c(40,174),ylim=range(lat[[ii0[p]]]),xlab='',ylab='Latitude',main='Empirical correlation',cex=sc,cex.axis=sc,cex.main=3.5,cex.lab=sc,asp=1)
for (i0 in 1:34){
  for (j0 in 1:10){
    maxi <- max(unlist(index_x_mv[[ii0[p]]]))
    maxj <- max(unlist(index_y_mv[[ii0[p]]]))
    ix <- unlist(index_x_mv[[ii0[p]]][i0])
    jy <- unlist(index_y_mv[[ii0[p]]][j0])
    err0 <- Ut[[ii0[p]]] - Ur[[ii0[p]]]
    err0[err0==0] <- NA 
    err <- log(err0)[ix,jy,]
    ie <- which(is.na(err[,,1]),arr.ind=TRUE)
    E <- Map3Dto2D(err)
    IE <- which(is.na(E[,1]))
    ix0 <- ix ; jy0 <- jy
    if (i0>1){
      if (j0>1 & j0<10){
        ix0 <- c(max((min(ix)-10),1):(min(ix)-1),ix,(max(ix)+1):min((max(ix)+10),maxi))
        jy0 <- c(max((min(jy)-10),1):(min(jy)-1),jy,(max(jy)+1):min((max(jy)+10),maxj))
        err <- log(err0)[ix0,jy0,]
        ie <- which(is.na(err[,,1]),arr.ind=TRUE)
        E <- Map3Dto2D(err)
      }
    }
    if ( length(IE) == 0 ){ 
      x0 <- floor(length(jy0)/2)*length(ix0) + floor(length(ix0)/2) + 1
      lon0 <- (long[[ii0[p]]][ix0])
      lat0 <- (lat[[ii0[p]]][jy0])
      latgr <- repLat(lat0,times=length(lon0))
      longr <- rep(lon0,times=length(lat0))
      covEmp <- cov(t(E)) 
      corEmp <- covEmp / (sqrt(diag(covEmp)) %*% t(sqrt(diag(covEmp))))
      c0 <- corEmp[x0,]
      mc0 <- matrix(c0,length(ix0),length(jy0))
      contour(lon0,lat0,mc0,level=.65,lwd=2,drawlabels = FALSE,add=TRUE) }
  }
}
map("world",add=TRUE,col="grey",border='grey',fill=T)
map('world',interior=F,add=T)
#####
plot(NA,xlim=c(40,174),ylim=range(lat[[ii0[p]]]),xlab='',ylab='Latitude',main='Simulation correlation - scale-aware model',cex=sc,cex.axis=sc,cex.main=3.5,cex.lab=sc,asp=1)
for (i0 in 1:34){
  for (j0 in 1:10){
    ix <- unlist(index_x_mv[[ii0[p]]][i0])
    jy <- unlist(index_y_mv[[ii0[p]]][j0])
    err0 <- Ut[[ii0[p]]] - Ur[[ii0[p]]]
    err0[err0==0] <- NA 
    err <- log(err0)[ix,jy,]
    E <- Map3Dto2D(err)
    IE <- which(is.na(E[,1]),arr.ind=TRUE)
    ix0 <- ix ; jy0 <- jy
    if (i0>1){
      if (j0>1 & j0<10){
        ix0 <- c(max((min(ix)-10),1):(min(ix)-1),ix,(max(ix)+1):min((max(ix)+10),maxi))
        jy0 <- c(max((min(jy)-10),1):(min(jy)-1),jy,(max(jy)+1):min((max(jy)+10),maxj))
        err <- log(err0)[ix0,jy0,]
        ie <- which(is.na(err[,,1]),arr.ind=TRUE)
        E <- Map3Dto2D(err)
      }
    }
    if (length(IE) == 0){
      x0 <- floor(length(jy0)/2)*length(ix0) + floor(length(ix0)/2) + 1
      lon0 <- (long[[ii0[p]]][ix0])
      lat0 <- (lat[[ii0[p]]][jy0])
      latgr <- repLat(lat0,times=length(lon0))
      longr <- rep(lon0,times=length(lat0))
      covSim0 <- cov(t(Map3Dto2D(Y0sim[ix0,jy0,,1]))) 
      corSim0 <- covSim0 / (sqrt(diag(covSim0)) %*% t(sqrt(diag(covSim0))))
      c0 <- corSim0[x0,]
      mc0 <- matrix(c0,length(ix0),length(jy0))
      contour(lon0,lat0,mc0,level=.65,lwd=2,drawlabels = FALSE,add=TRUE) }
  }
}
map("world",add=TRUE,col="grey",border='grey',fill=T)
map('world',interior=F,add=T)
#####
plot(NA,xlim=c(40,174),ylim=range(lat[[ii0[p]]]),xlab='Longitude',ylab='Latitude',main='Simulation correlation - scale-specific model',cex=sc,cex.axis=sc,cex.main=3.5,cex.lab=sc,asp=1)
for (i0 in 1:34){
  for (j0 in 1:10){
    ix <- unlist(index_x_mv[[ii0[p]]][i0])
    jy <- unlist(index_y_mv[[ii0[p]]][j0])
    err0 <- Ut[[ii0[p]]] - Ur[[ii0[p]]]
    err0[err0==0] <- NA 
    err <- log(err0)[ix,jy,]
    E <- Map3Dto2D(err)
    err <- log(err0)[ix,jy,]
    E <- Map3Dto2D(err)
    IE <- which(is.na(E[,1]),arr.ind=TRUE)
    ix0 <- ix ; jy0 <- jy
    if (i0>1){
      if (j0>1 & j0<10){
        ix0 <- c(max((min(ix)-10),1):(min(ix)-1),ix,(max(ix)+1):min((max(ix)+10),maxi))
        jy0 <- c(max((min(jy)-10),1):(min(jy)-1),jy,(max(jy)+1):min((max(jy)+10),maxj))
        err <- log(err0)[ix0,jy0,]
        ie <- which(is.na(err[,,1]),arr.ind=TRUE)
        E <- Map3Dto2D(err)
      }
    }
    if (length(IE) == 0){
      x0 <- floor(length(jy0)/2)*length(ix0) + floor(length(ix0)/2) + 1
      lon0 <- (long[[ii0[p]]][ix0])
      lat0 <- (lat[[ii0[p]]][jy0])
      latgr <- repLat(lat0,times=length(lon0))
      longr <- rep(lon0,times=length(lat0))
      covSim1 <- cov(t(Map3Dto2D(Y1sim[ix0,jy0,,1]))) 
      corSim1 <- covSim1 / (sqrt(diag(covSim1)) %*% t(sqrt(diag(covSim1))))
      c0 <- corSim1[x0,]
      mc0 <- matrix(c0,length(ix0),length(jy0))
      contour(lon0,lat0,mc0,level=.65,lwd=2,drawlabels = FALSE,add=TRUE) }
  }
}
map("world",add=TRUE,col="grey",border='grey',fill=T)
map('world',interior=F,add=T)
}
dev.off()


