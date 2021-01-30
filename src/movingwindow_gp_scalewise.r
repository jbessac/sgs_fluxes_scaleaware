
rm(list=ls())

indnx <- c(3,5,7,9,11,14,21,28,35,42,49)
n <- length(indnx)
k <- 7

load(file='residuals_scale.RData')
source('functions.r')

load(file='lat_list_res.RData')  #lat
load(file='long_list_res.RData') #long

load('index_x_mv_allscales.RData') #index_x_mv
load('index_y_mv_allscales.RData') #index_y_mv

load('lambda0_phys.RData')
lambda1 = lambda0

failed_optim = 0

H <- 6

i0 <- c()
j0 <- c()
d1 <- length(index_x_mv[[1]])
d2 <- length(index_y_mv[[1]])
for (i in 1:d1){
    for (j in 1:d2){
    ind1 <- unlist(index_x_mv[[k]][[i]])
    ind2 <- unlist(index_y_mv[[k]][[j]])

    model = paste('nx',indnx[k],'_n1_long',i,'lat',j,sep='')
    filename = paste('optim_corr_spacetime_scale/singlescale/corrfit_residuals_',model,'.RData',sep='')
    lat_CR_grid = repLat(lat[[k]][ind2],times=length(long[[k]][ind1]))
    long_CR_grid = rep(long[[k]][ind1],times=length(lat[[k]][ind2]))
    R1 = Map3Dto2D(res[[k]][ind1,ind2,])
    
    par_opt <- rep(0,5)
    if (sum(is.na(c(R1)))>0){
       indNA <- which(is.na(R1),arr.ind=TRUE)
       if (length(unique(indNA[,1]))>(length(ind1)*length(ind2)-2)) { print('empty window') ; par_opt <- rep(NA,5) }
       else {  print('land-point')
       	   R1 <- as.matrix(R1[-unique(indNA[,1]),])
           R1 = MapSpaceTime2SpaceTimeDay(X=R1,H=H)
           dlat_CR_ST = outer(rep(lat_CR_grid[-unique(indNA[,1])],each=H),rep(lat_CR_grid[-unique(indNA[,1])],each=H),'-')
           dlong_CR_ST = outer(rep(long_CR_grid[-unique(indNA[,1])],each=H),rep(long_CR_grid[-unique(indNA[,1])],each=H),'-')
           indT = rep(1:H,times=(length(long_CR_grid[-unique(indNA[,1])])))
           dt = abs(outer(indT,indT,'-'))
           }      }
       
    else { print('sea-points only')
    	 R1 = MapSpaceTime2SpaceTimeDay(X=R1,H=H)
    	 dlat_CR_ST = outer(rep(lat_CR_grid,each=H),rep(lat_CR_grid,each=H),'-')
         dlong_CR_ST = outer(rep(long_CR_grid,each=H),rep(long_CR_grid,each=H),'-')
    	 indT = rep(1:H,times=(length(long_CR_grid)))
         dt = abs(outer(indT,indT,'-'))
         }
    
   opt_corr = tryCatch({optim(par=lambda1,fn=MLE_Exp_SpaceTimeAlpha_phys,Y=R1,dLat=dlat_CR_ST,dLong=dlong_CR_ST,dt=dt,method="L-BFGS-B",control=list(trace=1,maxit=1000),lower=rep(0,5),upper=c(rep(+Inf,3),2,+Inf))}, error = function(e) NA)
	 opt_corr = as.list(opt_corr)
	 if (length(opt_corr)==1){
	     model0 = paste('nx',indnx[(k+1)],'_n2_long',1,'lat',1,sep='')
       filename0 = paste('optim_corr_spacetime_scale/singlescale/corrfit_residuals_',model0,'.RData',sep='')
	     load(filename0)
	     lambda0=opt_corr$par
	     opt_corr = tryCatch({optim(par=lambda0,fn=MLE_Exp_SpaceTimeAlpha_phys,Y=R1,dLat=dlat_CR_ST,dLong=dlong_CR_ST,dt=dt,method="L-BFGS-B",control=list(trace=1,maxit=1000),,lower=rep(0,5),upper=c(rep(+Inf,3),2,+Inf))}, error = function(e) NA)
	     opt_corr = as.list(opt_corr)
	     if (length(opt_corr)==1){
	     	i0 <- c(i0,i)
	    	j0 <- c(j0,j)
	     	model0 = paste('nx',indnx[(k+1)],'_n2_long',max(1,(i-1)),'lat',max(1,(j-1)),sep='')
            filename0 = paste('optim_corr_spacetime_scale/singlescale/corrfit_residuals_',model0,'.RData',sep='')
	     	load(filename0)
	     	lambda0 = opt_corr$par
	     	opt_corr = tryCatch({optim(par=lambda0,fn=MLE_Exp_SpaceTimeAlpha_phys,Y=R1,dLat=dlat_CR_ST,dLong=dlong_CR_ST,dt=dt,method="L-BFGS-B",control=list(trace=1,maxit=1000),,lower=rep(0,5),upper=c(rep(+Inf,3),2,+Inf))}, error = function(e) NA)
            opt_corr = as.list(opt_corr)         }
	 }
     if (length(opt_corr)==1){ failed_optim = failed_optim + 1 }
     save(failed_optim,file=paste('optim_corr_spacetime_scale/singlescale/failed_optim_nx',indnx[k],'_scale.RData',sep=''))
	save(opt_corr,file=filename)	  
    }}
# saving indexes of failed optimizations
save(i0,file=paste('optim_corr_spacetime_scale/singlescale/i0_index_nx',indnx[k],'.RData',sep=''))
save(j0,file=paste('optim_corr_spacetime_scale/singlescale/j0_index_nx',indnx[k],'.RData',sep=''))


