rm(list=ls())

### ------ extract u, v and prcp from netcdf files
source('src/extract_u.r') # this will also extract latitude and longitude
source('src/extract_v.r')
source('src/extract_prcp.r')



### ------ resolution indexes
indnx <- c(3,5,7,9,11,14,21,28,35,42,49) # each number nx such that nx*nx grid-points are averaged during coarsening
indnxn <- indnx/max(indnx)
n <- length(indnx)
scalenx <- round(indnx*4/111.11,2)
scale_name <- c()
for (i in 1:length(indnx)){
  scale_name <- c(scale_name,paste('N =',scalenx[i])) }



### ------ coarsening of data
# define flux exponent
n <- 1
# define coarsening level
# this has to be run for all k=1:length(indnx) ; not to run serially, to be parallelized
k <- 7
nx <- indnx[k]  # the coarsening will average nx=21*nx=21 grid-points into a grid-box
source('coarsen_data.r') # need n and nx to be defined outside 



### ------ save data (list of 3D-arrays, each list item corresponds to a given resolution)
Ur <- NULL ; Ut <- NULL
for (i in 1:length(indnx)){
    filename <- paste('prcp_nx',indnx[i],'.RData',sep='')
    load(filename)
    P[[i]] <- prcp[,,]
    rm(prcp)
    filename <- paste('ws_res_nx',indnx[i],'_n',n,'.RData',sep='')
    load(filename)
    Ur[[i]] <- ws_res[,,]
    rm(ws_res)
    filename <- paste('ws_tr_nx',indnx[i],'_n',n,'.RData',sep='')
    load(filename)
    Ut[[i]] <- ws_tr[,,]
    rm(ws_tr)    }
save(P,file='prcp.RData')
save(Ur,file='Ures_n1.RData')
save(Ut,file='Utrue_n1.RData')



### ------ save moving-windows indexes (list of )
ds <- 4  # corresponds to 500km (moving-window size)
load('lat_v.RData')
load('long_u.RData')
lat_range <- max(lat_v) - min(lat_v)
long_range <- max(long_u) - min(long_u)
ind_long <- seq(min(long_u),max(long_u),by=ds)
d1 <- length(ind_long)
ind_lat <- seq(min(lat_v),max(lat_v),by=ds)
d2 <- length(ind_lat)
#
# define indexes on x-axis
for (k in 1:length(indnx)){
index_x_mv_scales <- c()
index_x_mv_scales[[1]] <- list(which((long[[k]]>=(ind_long[1]))&(long[[k]]<=(ind_long[2]+.4))))
for (i in 2:(d1-1)){
    index_x_mv_scales[[i]] <- list(which((long[[k]]>=(ind_long[i]-.4))&(long[[k]]<=(ind_long[(i+1)]+.4))))  }
index_x_mv_scales[[d1]] <- list(which((long[[k]]>=(ind_long[d1]-.4))&(long[[k]]<=max(long_u))))
save(index_x_mv_scales,file=paste('index_x_mv_nx',indnx[k],'.RData',sep=''))    }
#
index_x_mv <- NULL
for(k in 1:length(indnx)){
    load(file=paste('index_x_mv_nx',indnx[k],'.RData',sep=''))
    index_x_mv[[k]] <- index_x_mv_scales  }
save(index_x_mv,file='index_x_mv_allscales.RData')
#
# define indexes on x-axis
for (k in 1:length(indnx)){
index_y_mv_scales <- c()
index_y_mv_scales[[1]] <- list(which((lat[[k]]>=(ind_lat[1]))&(lat[[k]]<=(ind_lat[2]+.4))))
for (i in 2:(d2-1)){
    index_y_mv_scales[[i]] <- list(which((lat[[k]]>=(ind_lat[i]-.4))&(lat[[k]]<=(ind_lat[(i+1)]+.4))))  }
index_y_mv_scales[[d2]] <- list(which((lat[[k]]>=(ind_lat[d2]-.4))&(lat[[k]]<=max(lat_v))))
save(index_y_mv_scales,file=paste('index_y_mv_nx',indnx[k],'.RData',sep=''))    }
#
index_y_mv <- NULL
for(k in 1:length(indnx)){
    load(file=paste('index_y_mv_nx',indnx[k],'.RData',sep=''))
    index_y_mv[[k]] <- index_y_mv_scales  }
save(index_y_mv,file='index_y_mv_allscales.RData')
#
# save lat-lon coordinates of moving-windows
long_mv <- (ind_long[1:(length(ind_long)-1)]+ind_long[2:(length(ind_long))])/2
lat_mv <- (ind_lat[1:(length(ind_lat)-1)]+ind_lat[2:(length(ind_lat))])/2
lat_mv <- c(lat_mv,lat_mv[9]+4)
save(long_mv,file='long_mv.RData')
save(lat_mv,file='lat_mv.RData')

