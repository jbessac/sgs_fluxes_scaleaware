rm(list=ls())

library('fields')

source('functions.r')


### ------ resolution indexes
indnx <- c(3,5,7,9,11,14,21,28,35,42,49) # nx such that nx*nx grid-points are averaged during coarsening
indnxn <- indnx/max(indnx)
n <- length(indnx)
scalenx <- round(indnx*4/111.11,2)
scale_name <- c()
for (i in 1:length(indnx)){
  scale_name <- c(scale_name,paste('N =',scalenx[i])) }


### ------
regression <- 0        # regression <- 1 to perform whole-domain scale-wise regression
regression_scale <- 0  # regression_scale <- 1 to perform whole-domain scale-aware regression
fit_gp <- 0            # fit_gp <- 1 to fit scale-wise Gaussian process on regression residuals
fit_gp_scale <- 0            # fit_gp_scale <- 1 to fit scale-aware Gaussian process on regression residuals



### ------ load data 
load('prcp.RData')
load('Ures_n1.RData')   # Ur
load('Utrue_n1.RData')  # Ut
load(file='lat_list_res.RData')  #lat
load(file='long_list_res.RData') #long
load('index_x_mv_allscales.RData') #index_x_mv
load('index_y_mv_allscales.RData') #index_y_mv
d1 <- length(index_x_mv[[1]]) # number of windows on x-axis
d2 <- length(index_y_mv[[1]]) # number of windows on y-axis



### ------ fit whole-domain scale-wise regression and save parameters and residuals
if (regression == 1){
  res <- NULL
  coeff_reg <- array(0,c(8,length(indnx)))
  for (i in 1:length(indnx)){
    err0 <- Ut[[i]] - Ur[[i]]
    err0[err0==0] <- NA 
    err <- log(err0)
    rm(err0)
    x1 <- log(Ur[[i]])
    x2 <- P[[i]]
    res.lm <- lm(c(err)~1+c(x1)+c(x1^2)+c(x1^3))  # regression without precipitation
    #res.lm <- lm((c(err))~1+(c(x1))+(c(x1^2))+(c(x1^3))+(c(x2))+(c(x2^(1/4)))+(c(x2^(1/2)))+(c(x2^(3/4))))  # regression with precipitation
    c0 <- res.lm$coefficients[1] ; c1 <- res.lm$coefficients[2]
    c2 <- res.lm$coefficients[3] ; c3 <- res.lm$coefficients[4]
    #c4 <- res.lm$coefficients[5] ; c5 <- res.lm$coefficients[6] # regression with precipitation
    #c6 <- res.lm$coefficients[7] ; c7 <- res.lm$coefficients[8] # regression with precipitation
    coeff_reg[,i] <- res.lm$coefficients
    res[[i]] <- err - c0 - c1*x1 - c2*x1^2 - c3*x1^3
    #res[[i]] <- err - c0 - c1*x1 - c2*x1^2 - c3*x1^3 - c4*x2 - c5*x2^(1/4) - c6*x2^(1/2) - c7*x2^(3/4) # regression with precipitation
    rm(err,x1,x2)   
    save(res,file='residuals_noprcp_n1.RData')    }
  save(coeff_reg,file='coeff_reg_noprcp_n1.RData')	 }



### ------ fit global regression with scale-aware parameterization and save parameters and residuals
if (regression_scale == 1) {
    # fit scale-aware regression on selected resolutions
    i1 = 5 ; i2 = 7 ; i3 = 10
    err0 <- c(c(Ut[[i1]]),c(Ut[[i2]]),c(Ut[[i3]])) - c(c(Ur[[i1]]),c(Ur[[i2]]),c(Ur[[i3]]))
    err0[err0==0] <- NA
    err <- log(err0)
    rm(err0)
    X1 <- log(c(c(Ur[[i1]]),c(Ur[[i2]]),c(Ur[[i3]])))
    X2 <- c(c(P[[i1]]),c(P[[i2]]),c(P[[i3]]))
    N <- (c(rep(indnxn[i1],length(c(P[[i1]]))),rep(indnxn[i2],length(c(P[[i2]]))),rep(indnxn[i3],length(c(P[[i3]])))))
    N2 <- (N*N)
    x1 <- (X1) ; x2 <- x1*N ; x3 <- x1*N2
    x4 <- (X1^2) ; x5 <- x4*N ; x6 <- x4*N2
    x7 <- (X1^3) ; x8 <- x7*N ; x9 <- x7*N2
    x10 <- (X2) ; x11 <-x10*N ; x12 <- x10*N2
    x13 <- (X2^.75) ; x14 <- x13*N ; x15 <- x13*N2
    x16 <- (X2^.5) ; x17 <- x16*N ; x18 <- x16*N2
    x19 <- (X2^.25) ; x20 <- x19*N ; x21 <- x19*N2
    scale.lm <- lm(c(err)~1+log(N)+N2+x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20+x21)
    rm(err,X1,X2,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21)
    save(scale.lm,file='regression_scaleaware_nx5710_n1.RData')
  
    # plot fitted parameters
    coeff = scale.lm$coefficients
    thetaN = matrix(coeff,8,3,byrow=TRUE)
    parN = thetaN %*% t(cbind(rep(1,n),indnxn,indnxn^2))
    parN[1,] = coeff[1] + coeff[2]*log(indnxn) + coeff[3]*indnxn*indnxn
    readline('press enter')
    sc = 1
    layout(matrix(1:2,1,2))
    par(mar=c(6,6,4.5,2.5))
    matplot(scalenx,t(coeff_reg[1:4,]),typ='l',col=1:4,lty=1,xlab='Scale N (degre)',ylab='',main='Coefficients A',cex=sc,cex.lab=sc,cex.axis=sc,cex.main=2,lwd=2)
    mtext('Regression coefficient',side=2,line=4,cex=sc)
    for (i in 1:4){
        lines(scalenx,parN[i,],col=i,lty=3,lwd=2)}
    legend(x='bottomright',lty=c(1,2,1,1,1,1),c('Empirical','Scale-aware regression',expression(a[0] (b[1])),expression(a[1] (b[2])),expression(a[2] (b[3])),expression(a[3] (b[4]))),col=c(1,1,1:4),cex=sc,bty='n')
    matplot(scalenx,t(coeff_reg[c(5,8,7,6),]),typ='l',col=1:4,lty=1,xlab='Scale N (degre)',ylab='',main='Coefficients B',cex=sc,cex.lab=sc,cex.axis=sc,cex.main=2,lwd=2)
    for (i in 1:4){
        lines(scalenx,parN[(i+4),],col=i,lty=3,lwd=2) }
  
    # save scale-aware regression residuals
    res3 <- NULL
    for (i in 1:length(indnx)){
        err0 <- Ut[[i]] - Ur[[i]]
        err0[err0==0] <- NA
        err <- log(err0)
        rm(err0)
        X1 <- log(Ur[[i]])
        X2 <- P[[i]]
        N <- indnx[i]
        N2 <- N*N
        res3[[i]] <- err - parN[1,i] - parN[2,i]*X1 - parN[3,i]*X1^2 - parN[4,i]*X1^3 - parN[5,i]*X2 - parN[6,i]*X2^.75 - parN[7,i]*X2^.5 - parN[8,i]*X2^.25
        rm(err,X1,X2)   }
    res <- res3
    save(res,file='residuals_scale_n1.RData') # saved for nx=(5,7,10)
}



### ------ fit scale-wise Gaussian process model
if (fit_gp == 1){
    # ------ fit GP in moving-window fashion
    # this has to be run for all k=1:length(indnx) ; not to run serially, to be parallelized
    # need to create a folder `optim_corr_spacetime_scale/singlescale'
    # 'movingwindow_gp_scalewise.r' is self-contained
    # specify k between 1 and length(indnx) inside
    source('movingwindow_gp_scalewise.r')
    
    # ------ save and plot space-time covariance parameters (single-scale model)
    par_ani1 <- NULL ; par_ani2 <- NULL ; par_temp <- NULL
    par_exp <- NULL ; par_var <- NULL
    for (k in length(indnx):1){
      par_ani1[[k]] <- matrix(NA,d1,d2)
      par_ani2[[k]] <- matrix(NA,d1,d2)
      par_temp[[k]] <- matrix(NA,d1,d2)
      par_exp[[k]] <- matrix(NA,d1,d2)
      par_var[[k]] <- matrix(NA,d1,d2)
      for (i in 1:(d1-1)){
        for (j in 1:d2){
          ind1 <- unlist(index_x_mv[[k]][[i]])
          ind2 <- unlist(index_y_mv[[k]][[j]])
          model = paste('nx',indnx[k],'_n1_long',i,'lat',j,sep='')
          filename = paste('optim_corr_spacetime_scale/singlescale/corrfit_residuals_',model,'_phys_scale.RData',sep='')
          R1 = Map3Dto2D(res[[k]][ind1,ind2,])
          if (sum(is.na(c(R1)))>0){
            indNA <- which(is.na(R1),arr.ind=TRUE)
            if ( length(unique(indNA[,1]))>(length(ind1)*length(ind2)-2) ) {print('empty window') ; l='empty' }
            else {  l = tryCatch( load(filename),error=function(e) 'empty')  }}
          else { l = tryCatch( load(filename),error=function(e) 'empty')}
          if (l != "empty"){
            opt_corr <- as.list(opt_corr)
            if (length(opt_corr)>1){
              par_opt = opt_corr$par
              par_ani1[[k]][i,j] <- par_opt[1] ; par_ani2[[k]][i,j] <- par_opt[2] ; par_temp[[k]][i,j] <- (par_opt[3]) ;
              par_exp[[k]][i,j] <- par_opt[4] ; par_var[[k]][i,j] <- par_opt[5]  }}}}}
    save(par_ani1,file='optim_corr_spacetime_scale/par_ani1_scalereg_n1.RData',version=2)
    save(par_ani2,file='optim_corr_spacetime_scale/par_ani2_scalereg_n1.RData',version=2)
    save(par_temp,file='optim_corr_spacetime_scale/par_temp_scalereg_n1.RData',version=2)
    save(par_exp,file='optim_corr_spacetime_scale/par_exp_scalereg_n1.RData',version=2)
    save(par_var,file='optim_corr_spacetime_scale/par_var_scalereg_n1.RData',version=2)
    #
    # plot boxplot of parameters
    layout(matrix(1:6,2,3,byrow=TRUE))
    boxplot((par_ani1),xlab='scale N',xaxs='n',main='Meridional anisotropy',xaxt='n',ylim=c(0,4.5),cex=2,cex.axis=2.4,cex.lab=2.4,cex.main=2.4,pch=20,lwd=2)
    axis(1,at=1:n,label=scalenx,cex.axis=2.4,cex.axis=2.4)
    #
    boxplot((par_ani2),xlab='scale N',xaxs='n',main='Zonal anisotropy',xaxt='n',ylim=c(0,4.5),cex=2,cex.axis=2.4,cex.lab=2.4,cex.main=2.4)
    axis(1,at=1:n,label=scalenx,cex.axis=2.4,cex.axis=2.4)
    #
    boxplot((par_temp),xlab='scale N',xaxs='n',main='Temporal decay',xaxt='n',ylim=c(0,15),cex=2,cex.axis=2.4,cex.lab=2.4,cex.main=2.4)
    axis(1,at=1:n,label=scalenx,cex.axis=2.4,cex.axis=2.4)
    #
    boxplot((par_exp),xlab='scale N',xaxs='n',main='Exponant',xaxt='n',ylim=range(par_exp,na.rm=TRUE),cex=2,cex.axis=2.4,cex.lab=2.4,cex.main=2.4)
    axis(1,at=1:n,label=scalenx,cex.axis=2.4,cex.axis=2.4)
    #
    boxplot((par_var),xlab='scale N',xaxs='n',main='Variance',xaxt='n',ylim=c(0,4),cex=2,cex.axis=2.4,cex.lab=2.4,cex.main=2.4)
    axis(1,at=1:n,label=scalenx)

    # save parameters in array format for future use in scale-aware GP model
    par1 <- array(NA,c(d1,d2,n)) ; par2 <- par1 ; par3 <- par1
    par4 <- par1 ; par5 <- par1
    for(i in 1:(d1-1)){
      for (j in 1:d2){
        for (k in 4:length(indnx)){
          par1[i,j,k] <- par_ani1[[k]][i,j]
          par2[i,j,k] <- par_ani2[[k]][i,j]
          par3[i,j,k] <- par_temp[[k]][i,j]
          par4[i,j,k] <- par_exp[[k]][i,j]
          par5[i,j,k] <- par_var[[k]][i,j]  }}}
    save(par1,file='optim_corr_spacetime_scale/par1_scalereg_n1.RData')
    save(par2,file='optim_corr_spacetime_scale/par2_scalereg_n1.RData')
    save(par3,file='optim_corr_spacetime_scale/par3_scalereg_n1.RData')
    save(par4,file='optim_corr_spacetime_scale/par4_scalereg_n1.RData')
    save(par5,file='optim_corr_spacetime_scale/par5_scalereg_n1.RData')
}



### ------ fit scale-aware Gaussian process model
if (fit_gp_scale == 1){
    # fit scale-aware GP model in moving-window fashion
    # this runs the fit for the entire domain, I suggest splitting the domain in 4 (specify range of i and j)
    # need to create a folder `optim_corr_spacetime_scale/nx5710'
    # 'movingwindow_gp_scaleaware.r' is self-contained
    source('movingwindow_gp_scaleaware.r')
    #
    # save parameters
    par1_estim = array(NA,c((d1-1),d2,n))
    par2_estim = par1_estim ; par3_estim = par1_estim
    par4_estim = par1_estim ; par5_estim = par1_estim
    for (i in 1:(d1-1)){
        for (j in 1:d2){
            filename <- paste('optim_corr_spacetime_scale/nx5710/corrfit_residuals_multisc_nx5710_n1_long',i,'lat',j,'_phys_scale.RData',sep='')
            load(filename)
            if (length(optim_multiscale_mle)>1){
                par_opt <- optim_multiscale_mle$par
                par1_estim[i,j,] <- (par_opt[1]*exp(par_opt[2]*indnxn))
                par2_estim[i,j,] <- (par_opt[3]*exp(par_opt[4]*indnxn))
                par3_estim[i,j,] <- (par_opt[5]*exp(par_opt[6]*indnxn))
                par4_estim[i,j,] <- (1+tanh((par_opt[7])+((par_opt[8])*indnxn)))
                par5_estim[i,j,] <- (par_opt[9])*exp((par_opt[10])*indnxn) } }}
    save(par1_estim,par2_estim,par3_estim,par4_estim,par5_estim,file='par_estim_nx5710_n1.RData')
    }



### ------ simulate scale-aware GP samples (psi); reconstruct flux erros (epsilon) and compute MSE, squared-bias and c-MSE
source('sample_mse_generation.r')
