rm(list=ls())

library('ncdf4')

#------------------------------------ extract v
# average data to account for the staggering of the grid
indV = 14:120
V = array(0,c(3862,1111,(2*length(indV))))
for (i in 1:length(indV)){
    filename = paste('Data/CASCADE_WarmPool-4km_xfhfc_p',indV[i],'_3_sfc.nc',sep='')
    nn = nc_open(filename)
    V[,,(2*i-1):(2*i)] = ncvar_get(nn)
    nc_close(nn)
    rm(nn)    }
V = (V[1:3861,,] + V[2:3862,,])/2
#save(V,file='V_av.RData') # 3861 * 1111 zonal-averaged grid


##--------------------------- mask land
load('landmask.RData')

V = V[1:3861,,]
m = m[1:3861,1:1111]
m[(m!=0)]=1

dimp = dim(V)
Vc = as.matrix(V)
dim(Vc) = c(dimp[1]*dimp[2], dimp[3])
mv = as.logical(as.vector(m))
Vc[mv, ] = NA
dim(Vc) = dimp
save(Vc,file='V_sea.RData')
