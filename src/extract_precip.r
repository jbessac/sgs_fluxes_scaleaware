
rm(list=ls())

library('ncdf4')


##--------------------------- extract from netcdf files
# average data to account for the staggering of the grid
indU = 14:120 # temporal indexes 
prcp = array(0,c(3862,1112,(2*length(indU))))
for (i in 1:length(indU)){
    filename = paste('Data/CASCADE_WarmPool-4km_xfhfc_p',indU[i],'_5216_15min-mean.nc',sep='')
    pp0 = nc_open(filename)
    pp = ncvar_get(pp0)
    prcp[,,(2*i-1)] = apply(X=pp[,,1:4],FUN=mean,MARGIN=1:2)
    prcp[,,(2*i)] = apply(X=pp[,,5:8],FUN=mean,MARGIN=1:2)
    rm(pp,pp0)    }
prcp = (prcp[,1:1111,] + prcp[,2:1112,])/2
#save(prcp,file='prcp_av.RData')


##--------------------------- mask land
load('landmask.RData')

prcp = prcp[1:3861,,]
m = m[1:3861,1:1111]
m[(m!=0)]=1

dimp = dim(prcp)
prcpc= as.matrix(prcp)
dim(prcpc) = c(dimp[1]*dimp[2], dimp[3])
mv = as.logical(as.vector(m))
prcpc[mv, ] = NA
dim(prcpc) = dimp
save(prcpc,file='prcp_sea.RData')



