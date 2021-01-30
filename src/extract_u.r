rm(list=ls())

library('ncdf4')

##------------------------------------ extract u
# average data to account for the staggering of the grid
indU = 14:120 # temporal indexes 
U = array(0,c(3862,1112,(2*length(indU))))
for (i in 1:length(indU)){
    filename = paste('Data/CASCADE_WarmPool-4km_xfhfc_p',indU[i],'_2_sfc.nc',sep='')
    nn = nc_open(filename)
    U[,,(2*i-1):(2*i)] = ncvar_get(nn)
    nc_close(nn)
    rm(nn)    }
U = (U[,1:1111,] + U[,2:1112,])/2
#save(U,file='U_av.RData') # U_av and V_av are over the all domain


load('lat_v.RData') # all domain on rematched-initially staggered grid 
load('long_u.RData') # all domain on rematched-initially staggered grid 


##--------------------------- mask land
load('landmask.RData')

U = U[1:3861,,]
m = m[1:3861,1:1111]
m[(m!=0)]=1

dimp = dim(U)
Uc = as.matrix(U)
dim(Uc) = c(dimp[1]*dimp[2], dimp[3])
mv = as.logical(as.vector(m))
Uc[mv, ] = NA
dim(Uc) = dimp
save(Uc,file='U_sea.RData')


##-------------------------------------- coarsen data

#load('U_sea.RData')
#load('V_sea.RData')
#load('landmask.RData')
#U = U[1:3861,,]
#m = m[1:3861,1:1111]

#m[(m!=0)]=1
#Uc = U
#Vc = V
#for (t in 1:dim(V)[3]){
# uc = Vc[,,t]
# uc[(m!=0)] = NA
# Vc[,,t] = uc
# rm(uc) }
#save(Vc,file='V_sea.RData')
#break

load('U_sea.RData') #Uc
load('V_sea.RData') #Vc

load('lat_v.RData')
load('long_u.RData')
lat=lat_v
long=long_u[1:3861]




#-------------------------------------------------------------------------------- extract data
#------------------------------------ extract latitude and longtidue
#n = nc_open('CASCADE_WarmPool-4km_xfhfc_p1_30.nc')
#lat = ncvar_get(n,'latitude')
#long = ncvar_get(n,'longitude')
#m = ncvar_get(n,'land_binary_mask')
#save(m,file='landmask.RData')
#save(lat,file='latitude.RData')
#save(long,file='longitude.RData')

#image.plot(long,lat,m)
#title('Land-sea Mask')

#nn = nc_open("Data/CASCADE_WarmPool-4km_xfhfc_p13_2_sfc.nc")
#u = ncvar_get(nn)
#lat_u = ncvar_get(nn,'latitude')
#long_u = ncvar_get(nn,'longitude')
#save(lat_u,file='lat_u.RData')

#nn = nc_open("Data/CASCADE_WarmPool-4km_xfhfc_p13_3_sfc.nc")
#v = ncvar_get(nn)
#lat_v = ncvar_get(nn,'latitude')
#long_v = ncvar_get(nn,'longitude')

#plot(lat[1:20],typ='l')
#lines(lat_u[1:20],col=2)
#lines(lat_v[1:20],col=3)

#plot(long[1:20],typ='l')
#lines(long_u[1:20],col=2)
#lines(long_v[1:20],col=3)



