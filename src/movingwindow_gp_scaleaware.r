
rm(list=ls())

indnx <- c(3,5,7,9,11,14,21,28,35,42,49)
indnxn <- indnx/max(indnx)
n <- length(indnx)
scalenx <- round(indnx*4/111.11,2)

load(file='residuals_scale_n1.RData')
source('functions.r')

load(file='lat_list_res.RData')  #lat
load(file='long_list_res.RData') #long

load('index_x_mv_allscales.RData') #index_x_mv
load('index_y_mv_allscales.RData') #index_y_mv
d1 <- length(index_x_mv[[1]])
d2 <- length(index_y_mv[[1]])

H <- 8

## Initial conditions from parameters at other scales
load('optim_corr_spacetime_scale/par_ani1_phys.RData') # par_ani1 - in physics space
load('optim_corr_spacetime_scale/par_ani2_phys.RData')
load('optim_corr_spacetime_scale/par_temp_phys.RData')
load('optim_corr_spacetime_scale/par_exp_phys.RData')
load('optim_corr_spacetime_scale/par_var_phys.RData')
med1 <- c() ; med2 <- c() ; med3 <- c() ; med4 <- c() ; med5 <- c()
for (i in 1:10){
  med1[i] = median(par_ani1[[(i+1)]],na.rm=TRUE)
  med2[i] = median(par_ani2[[(i+1)]],na.rm=TRUE)
  med3[i] = median(par_temp[[(i+1)]],na.rm=TRUE)
  med4[i] = median(par_exp[[(i+1)]],na.rm=TRUE)
  med5[i] = median(par_var[[(i+1)]],na.rm=TRUE)    }
ii1 <- indnxn[2:11]
lm1 <- lm(log(c(med1))~1+ii1)
lm2 <- lm(log(c(med2))~1+ii1)
lm3 <- lm(log(c(med3))~1+ii1)
lm4 <- lm((c(med4))~1+ii1)
lm5 <- lm(log(c(med5))~1+ii1)
lm_exp <- lm((abs(atanh(lm4$fitted-1))) ~ ii1)
a0 <- lm_exp$coefficients[1]
a1 <- lm_exp$coefficients[2]
lambda1 <- c(exp(lm1$coefficients[1]),(lm1$coefficients[2]),exp(lm2$coefficients[1]),(lm2$coefficients[2]),exp(lm3$coefficients[1]),(lm3$coefficients[2]),a0,a1,exp(lm5$coefficients[1]),((lm5$coefficients[2])))


load(file='optim_corr_spacetime_scale/par1_phys.RData') # par1 in parameters space
load(file='optim_corr_spacetime_scale/par2_phys.RData')
load(file='optim_corr_spacetime_scale/par3_phys.RData')
load(file='optim_corr_spacetime_scale/par4_phys.RData')
load(file='optim_corr_spacetime_scale/par5_phys.RData')

# choice of 3 resolutions used to fit scale-aware gp
ii <- c(5,7,10)
indnx0 <- indnxn[ii]
i0 <- c()
j0 <- c()
for (i in 1:floor(d1/2)){
    for (j in (1+floor(d2/2)):d2){
        Y <- NULL
        dt <- NULL
        dLat <- NULL
        dLong <- NULL
        for (k in 1:length(indnx0)){
            ind1 <- unlist(index_x_mv[[ii[k]]][[i]])
            ind2 <- unlist(index_y_mv[[ii[k]]][[j]])
            lat_CR_grid <- repLat((lat[[(ii[k])]])[ind2],times=length(long[[ii[k]]][ind1]))
            long_CR_grid <- rep((long[[(ii[k])]])[ind1],times=length(lat[[ii[k]]][ind2]))
            R0 <- Map3Dto2D(res[[ii[k]]][ind1,ind2,])
            if (sum(is.na(c(R0)))>0){
                indNA <- which(is.na(R0),arr.ind=TRUE)
                if (length(unique(indNA[,1]))>(length(ind1)*length(ind2)-2)) {print('empty window') ; par_opt <- rep(NA,5)}
                else { R0 <- R0[-unique(indNA[,1]),]
                    dLat[[k]] <- outer(rep(lat_CR_grid[-unique(indNA[,1])],each=H),rep(lat_CR_grid[-unique(indNA[,1])],each=H),'-')
                    dLong[[k]] <- outer(rep(long_CR_grid[-unique(indNA[,1])],each=H),rep(long_CR_grid[-unique(indNA[,1])],each=H),'-')
                    indT <- rep(1:H,times=(length(long_CR_grid[-unique(indNA[,1])]))) }}
            else { dLat[[k]] <- outer(rep(lat_CR_grid,each=H),rep(lat_CR_grid,each=H),'-')
                dLong[[k]] <- outer(rep(long_CR_grid,each=H),rep(long_CR_grid,each=H),'-')
                indT <- rep(1:H,times=(length(long_CR_grid))) }
            dt[[k]] <- abs(outer(indT,indT,'-'))
            Y[[k]] <- MapSpaceTime2SpaceTimeDay(X=R0,H=H)   }
        
        med1 <- par1[i,j,2:11] ; med2 <- par2[i,j,2:11] ; med3 <- par3[i,j,2:11] ;
        med4 <- par4[i,j,2:11] ; med5 <- par5[i,j,2:11] ;
        ii1 <- indnxn[2:11]
        lm1 <- lm(log(c(med1))~1+ii1)
        lm2 <- lm(log(c(med2))~1+ii1)
        lm3 <- lm(log(c(med3))~1+ii1)
        lm_exp <- lm((abs(atanh(lm4$fitted-1))) ~ ii1)
        a0 <- lm_exp$coefficients[1]
        a1 <- lm_exp$coefficients[2]
        lambda0 <- c(exp(lm1$coefficients[1]),(lm1$coefficients[2]),exp(lm2$coefficients[1]),(lm2$coefficients[2]),exp(lm3$coefficients[1]),(lm3$coefficients[2]),a0,a1,exp(lm5$coefficients[1]),((lm5$coefficients[2])))
        
        xlower = c(0,0,0,0,0,0,-Inf,-Inf,0,-Inf)
        xupper = c(+Inf,+Inf,+Inf,+Inf,+Inf,+Inf,+Inf,+Inf,+Inf,0)
        optim_multiscale_mle <- tryCatch(optim(par=lambda0,fn=MLE_Exp_SpaceTime_multiscale_diag_phys,Y=Y,dLat=dLat,dLong=dLong,dt=dt,indnx=indnx0,method="L-BFGS-B",control=list(trace=1,maxit=1000),lower=xlower,upper=xupper),error=function(e) NA)
        optim_multiscale_mle <- as.list(optim_multiscale_mle)
        if (length(optim_multiscale_mle)==1) {optim_multiscale_mle <- tryCatch(optim(par=lambda1,fn=MLE_Exp_SpaceTime_multiscale_diag_phys,Y=Y,dLat=dLat,dLong=dLong,dt=dt,indnx=indnx0,method="L-BFGS-B",control=list(trace=1,maxit=1000),lower=xlower,upper=xupper),error=function(e) NA)
            optim_multiscale_mle <- as.list(optim_multiscale_mle) }
        if (length(optim_multiscale_mle)==1) { filename <- paste('optim_corr_spacetime_scale/nx79_phys/corrfit_residuals_multisc_nx79_n1_long',i,'lat',j,'_phys.RData',sep='')
                   l = tryCatch( load(filename),error=function(e) 'empty')
                   if (l != "empty"){
                       lambda2 <- optim_multiscale_mle$par
                       optim_multiscale_mle <- tryCatch(optim(par=lambda2,fn=MLE_Exp_SpaceTime_multiscale_diag_phys,Y=Y,dLat=dLat,dLong=dLong,dt=dt,indnx=indnx0,method="L-BFGS-B",control=list(trace=1,maxit=1000),lower=xlower,upper=xupper),error=function(e) NA)
                       optim_multiscale_mle <- as.list(optim_multiscale_mle) }
                   else (optim_multiscale_mle <- NA) }
        if (length(optim_multiscale_mle)==1) { print('optim failed') ; print(c(i,j)) }
        i0 <- c(i0,i)
        j0 <- c(j0,j)
        filename <- paste('optim_corr_spacetime_scale/nx5710/corrfit_residuals_multisc_nx5710_n1_long',i,'lat',j,'_phys_scale.RData',sep='')
        save(optim_multiscale_mle,file=filename)
    }}
save(i0,file='optim_corr_spacetime_scale/i0_p11_nx5710_phys_n1.RData')
save(j0,file='optim_corr_spacetime_scale/j0_p11_nx5710_phys_n1.RData')
