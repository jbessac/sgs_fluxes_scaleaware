
###----------------------- coarsening

coarsen = function(X,nx,ny){
 d = dim(X)
 dx = floor(d[1]/nx)
 dy = floor(d[2]/ny)
 dt = d[3]
 x = array(0,c(dx,dy,dt))
 for (i in 1:dx){
    for (j in 1:dy){
    	x[i,j,] = apply(X=X[(1+(i-1)*nx):(i*nx),(1+(j-1)*ny):(j*ny),],FUN=mean,MARGIN=3,na.rm=TRUE) }}
 x }



###----------------------- coordinates and distances
distance = function(lat1,long1,lat2,long2){
    R = 6371
    dlat = rad(lat1)-rad(lat2)
    dlong = rad(long1)-rad(long2)
    a = sin(dlat/2)*sin(dlat/2)+cos(rad(lat1))*cos(rad(lat2))*sin(dlong/2)*sin(dlong/2)
    c = 2*atan2(sqrt(a),sqrt(1-a))
    d = R * c
    d }

distanceSites = function(lat,long,names=NULL){
    P = cbind(lat,long)
    D = matrix(0,length(lat),length(lat))
    for ( i in 1:(length(lat)-1) ){
        for ( j in (i+1):length(lat) ){
            D[i,j] = distance(P[i,1],P[i,2],P[j,1],P[j,2]) }}
    D = D + t(D)
    if (missing(names)){
        colnames(D) = as.character(1:length(lat))
        rownames(D) = as.character(1:length(lat)) }
    colnames(D) = names
    rownames(D) = names
    D	}

repLat = function(Lat,times){
   l = length(Lat)
   Lat_rep = c()
   for (i in 1:l){
        Lat_rep[(1+(i-1)*times):(i*times)] = rep(Lat[i],times=times) }
   Lat_rep  }

diffCoord = function(coord){
    dcoord = outer(coord,coord,'-')
    dcoord }


###----------------------- data structure 
Map3Dto2D = function(Y){
  d = dim(Y)
  d1 = d[1]
  d2 = d[2]
  d3 = d[3]
  y = array(0,c((d1*d2),d3))
  for (t in 1:d3){
      y[,t] = c(Y[,,t])}
      y	    }

Map2Dto3D = function(Y,dx){
  d = dim(Y)
  d1 = d[1]
  dt = d[2]
  dy = floor(d1/dx)
  y = array(0,c(dx,dy,dt))
  for (t in 1:dt){
      y[,,t] = matrix(Y[1:(dx*dy),t],dx,dy)}
  y	    }

MapSpaceTime2SpaceTimeDay = function(X,H){
 d = dim(X)
 d1 = d[1]
 d2 = d[2]
 Nd = floor(d2/H)
 x = matrix(0,(d1*H),Nd)
 for (i in 1:Nd){
     x[,i] = c(t(X[,(1+(i-1)*H):(i*H)])) }
 x     }


###----------------------- spatio-temoporal covariance

corr_Exp_SpaceTimeAlpha_phys = function(lambda,dLat,dLong,dt){ #parameters coded as in physical space 
  sig0 = lambda[5]
  dd = sqrt((dLat/lambda[1])^2 + (dLong/lambda[2])^2 + (dt/(lambda[3]))^2)^lambda[4]
  c0 = exp(-dd)     
  c0 = c0 * sig0
  #diag(c0) = diag(c0) #+ lambda[5]
  c0 }


###----------------------- spatio-temporal MLE

MLE_Exp_SpaceTimeAlpha_phys <- function(lambda,Y,dLat,dLong,dt){
  K <- dim(Y)[1]
  N <- dim(Y)[2]
  c0 <- corr_Exp_SpaceTimeAlpha_phys(lambda,dLat,dLong,dt)
  c0.c <- t(chol(c0))
  mu0 <- rowMeans(Y)
  out <- forwardsolve(c0.c,(Y-matrix(rep(mu0,N),K,N)))
  log.det3 <- 2*sum(log(diag(c0.c)))
  quadratic3 <- sum(out^2,na.rm=TRUE)
  return( N*log.det3+quadratic3 )  }

###----------------------- spatio-temporal with diagonal multiscale MLE


MLE_Exp_SpaceTime_multiscale_diag_phys <- function(lambda,Y,dLat,dLong,dt,indnx){			  
  n <- length(Y)
  theta1 <- (lambda[1])*exp((lambda[2])*indnx)  
  theta2 <- (lambda[3])*exp((lambda[4])*indnx)
  theta3 <- (lambda[5])*exp((lambda[6])*indnx)
  theta4 <- (1+tanh((lambda[7])+((lambda[8])*indnx))) 
  theta5 <- (lambda[9])*exp((lambda[10])*indnx)  
  theta <- cbind(theta1,theta2,theta3,theta4,theta5)
  log.det <- c()
  quadratic <- c()
  K <- c()
  for ( k in 1:n ){ 
    K[k] <- dim(Y[[k]])[1]     # space-time dimension
    N <- dim(Y[[k]])[2]        # number of samples ('days') - constant across scale k
    c0 <- corr_Exp_SpaceTimeAlpha_phys(theta[k,],dLat[[k]],dLong[[k]],dt[[k]])
    c0.c <- t(chol(c0))
    mu0 <- rowMeans(Y[[k]])
    out <- forwardsolve(c0.c,(Y[[k]]-matrix(rep(mu0,N),K[k],N)))
    quadratic[k] <- sum(out^2,na.rm=TRUE) 
    log.det[k] <- 2*sum(log(diag(c0.c))) } #+ K[k]*log(2*pi) }
  return( N*sum(log.det)+sum(quadratic) )   }  

