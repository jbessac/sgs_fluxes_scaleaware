
rm(list=ls())

library('fields')
library('maps')
library('ggplot2')
library('pals')

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

load(file='residuals_scale_n1.RData')
load('prcp.RData')
load('Ures_n1.RData')   # Ur
load('Utrue_n1.RData')

load(file='long_mv.RData') # long_mv
load(file='lat_mv.RData') # lat_mv

load('par_estim_nx5710_phys_scale_n1.RData') #par1_estim, par2_estim, par3_estim, par4_estim, par5_estim
load(file='optim_corr_spacetime_scale/par1_phys_scalereg_n1.RData') # par1 in parameters space
load(file='optim_corr_spacetime_scale/par2_phys_scalereg_n1.RData')
load(file='optim_corr_spacetime_scale/par3_phys_scalereg_n1.RData')
load(file='optim_corr_spacetime_scale/par4_phys_scalereg_n1.RData')
load(file='optim_corr_spacetime_scale/par5_phys_scalereg_n1.RData')

load('empty_window.RData') # ind_na
mat_window = matrix(1,(d1-1),d2)
mat_window[ind_na] = 0 

#break

ii0 <- c(3,6,11) # extrapolated scale to validate
H <- 24

k <- 1 # 1,2,3
load(file='regression_scaleaware_nx5710_n1.RData') #scale.lm
coeff <- scale.lm$coefficients
thetaN <- matrix(coeff,8,3,byrow=TRUE)
parN <- thetaN %*% t(cbind(rep(1,n),indnxn,indnxn^2))
parN[1,] <- coeff[1] + coeff[2]*log(indnxn) + coeff[3]*indnxn*indnxn
Nsim <- 30
Nd <- 8

#break 
sim = 1
if (sim==1){
  Y0sim <- array(NA,c(dim(res[[ii0[k]]])[1:2],(Nd*H),Nsim)) 
  Y1sim <- array(NA,c(dim(res[[ii0[k]]])[1:2],(Nd*H),Nsim))
for (i in 1:(d1-1)){
  for (j in 1:d2){ 
    if(mat_window[i,j] == 1){
    theta0 <- c(par1[i,j,ii0[k]],par2[i,j,ii0[k]],par3[i,j,ii0[k]],par4[i,j,ii0[k]],par5[i,j,ii0[k]])    # single-scale model parameters
    #theta0 <- c(par1_estim[i,j,ii0[k]],par2_estim[i,j,ii0[k]],par3_estim[i,j,ii0[k]],par4_estim[i,j,ii0[k]],par5_estim[i,j,ii0[k]])   # scale-aware model parameters
    print(theta0)
    if (!is.na(sum(theta0))){
    ind1 <- unlist(index_x_mv[[ii0[k]]][[i]])
    ind2 <- unlist(index_y_mv[[ii0[k]]][[j]])
    X1 <- log(Ur[[ii0[k]]])[ind1,ind2,]
    X2 <- P[[ii0[k]]][ind1,ind2,]
    N <- indnx[ii0[k]] ; N2 <- N*N
    xp <- parN[1,ii0[k]] + parN[2,ii0[k]]*X1 + parN[3,ii0[k]]*X1^2 + parN[4,ii0[k]]*X1^3 + parN[5,ii0[k]]*X2 + parN[6,ii0[k]]*X2^.75 + parN[7,ii0[k]]*X2^.5 + parN[8,ii0[k]]*X2^.25
    #
    lat_CR_grid <- repLat((lat[[(ii0[k])]])[ind2],times=length(long[[ii0[k]]][ind1]))
    long_CR_grid <- rep((long[[(ii0[k])]])[ind1],times=length(lat[[ii0[k]]][ind2]))
    dLat0 <- outer(rep(lat_CR_grid,each=H),rep(lat_CR_grid,each=H),'-')
    dLong0 <- outer(rep(long_CR_grid,each=H),rep(long_CR_grid,each=H),'-')
    indT <- rep(1:H,times=(length(long_CR_grid))) 
    dt0 <- abs(outer(indT,indT,'-'))
    #
    n0 <- (Nsim*Nd*((H*length(ind1)*length(ind2))))
    ee0 <- rnorm(n=n0,mean=0,sd=1) ; e0 <- array(ee0,c((H*length(ind1)*length(ind2)),Nd,Nsim))
    ee1 <- rnorm(n=n0,mean=0,sd=1) ; e1 <- array(ee1,c((H*length(ind1)*length(ind2)),Nd,Nsim))
    xy <- 1:(length(ind1)*length(ind2))
    xy0 <- xy
    xx <- xy0 %% length(ind1)
    xx[xx==0] <- length(ind1)
    yy <- ceiling(xy0/length(ind1))
    #
    S0 <- array(NA,c(length(ind1),length(ind2),(24*Nd),Nsim))
    cov0 <- corr_Exp_SpaceTimeAlpha_phys(theta0,dLat0,dLong0,dt0)  
    Cfit0 <- t(chol(cov0))
    #
    for (nd in 1:Nd){
      for (s in 1:Nsim){
        s0 <- Cfit0 %*% c(e0[,nd,s])
        for (h in 1:H){
          for (l in 1:length(xx)) {
            S0[xx[l],yy[l],((nd-1)*H+h),s] <- xp[xx[l],yy[l],((nd-1)*H+h)] + (s0[seq(h,length(s0),by=H)])[l]
            }}	    }}
    Y1sim[ind1,ind2,,] <- S0
  } } } }

save(Y1sim,file=paste('Ysim_nx',ii0[k],'singlescale_n1.RData',sep='')) }




### ------ calculation of MSE, squared bias and centered MSE
err0 <- Ut[[ii0[k]]] - Ur[[ii0[k]]]
err0[err0==0] <- NA 
err <- log(err0)
n1 = dim(err)[1]
n2 = dim(err)[2]

## squarred rmse 
s0 = 0 ; s1 = 0
for (s in 1:Nsim) {
  s0 = s0 + (apply(X=(err[,,1:192] - Y1sim[,,,s])^2,FUN=sum,MARGIN=1:2)) }
s0 = s0/(Nsim*192)
mse0 = s0 ; #mse1 = s1
filename = paste('mse_global_nx',ii0[k],'_singlescale_n1.RData',sep='')
save(mse0,file=filename)
## squarred centered rmse 
s0 = 0 ; s1 = 0 
sigf = apply(X=err[,,1:192],FUN=var,MARGIN=1:2)
corfr0 = matrix(NA,n1,n2) ; corfr1 = corfr0
for (s in 1:Nsim) {
  sigr0 = apply(X=Y1sim[,,,s],FUN=var,MARGIN=1:2)
  for (it in 1:dim(Y1sim)[1]){
    for (jt in 1:dim(Y1sim)[2]){
      corfr0[it,jt] = cor(Y1sim[it,jt,,s],err[it,jt,1:192]) }}
  s0 = s0 + (sigf + sigr0 - 2*sqrt(sigf*sigr0)*corfr0) }
s0 = s0/Nsim
cmse0 = s0
filename = paste('cmse_global_nx',ii0[k],'_singlescale_n1.RData',sep='')
save(cmse0,file=filename)
## bias
s0 = 0 ; s1 = 0 
muf = apply(X=err[,,1:192],FUN=mean,MARGIN=1:2)
for (s in 1:Nsim) {
  mur0 = apply(X=Y1sim[,,,s],FUN=mean,MARGIN=1:2)
  s0 = s0 + (muf - mur0)^2  }
s0 = s0/Nsim
bias0 = s0
filename = paste('bias_global_nx',ii0[k],'_singlescale_n1.RData',sep='')
save(bias0,file=filename)


