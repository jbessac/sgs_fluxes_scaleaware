source('functions.r')

## ---- precipitation coarsening
load(file='../prcp_Sea.RData') # prcpc has negative values
prcpc[prcpc<0] <- 0
prcp <- coarsen(X=prcpc,nx=nx,ny=nx)
filename <- paste('prcp_nx',nx,'.RData',sep='')
save(prcp,file=filename)
rm(prcpc,prcp,filename)

## ---- U and V coarsening
load('U_sea.RData') #Uc
load('V_sea.RData') #Vc

# compute true flux
ws_fr = sqrt(Uc*Uc + Vc*Vc)
ws_tr = coarsen(X=ws_fr^n,nx=nx,ny=nx)
filename <- paste('ws_tr_nx',nx,'_n',n,'.RData',sep='')
save(ws_tr,file=filename)
rm(ws_fr)

# compute resolved flux and flux error 
Uco <- coarsen(X=Uc,nx=nx,ny=nx)
Vco <-  coarsen(X=Vc,nx=nx,ny=nx)
ws_res = sqrt(Uco*Uco + Vco*Vco)^n
err = log((ws_tr - ws_res))
filename <- paste('ws_res_nx',nx,'_n',n,'.RData',sep='')
save(ws_res,file=filename)
filename <- paste('err_nx',nx,'_n',n,'.RData',sep='')
save(err,file=filename)
rm(ws_tr,ws_res,err,Uco,Vco)







