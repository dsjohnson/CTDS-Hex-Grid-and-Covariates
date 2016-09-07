#cov_parse_interp


# generate large matrix of covariates -- has n_rows = n_adj*n_t 

#### SURFACE WINDS

library(RNetCDF)
library(maps)
library(fields)

# U-COMPONENT 05
fname <- "~/Documents/Seaglider/Station_P/NCEP/uwnd.10m.gauss.2005.nc"
fid <- open.nc(fname)
# print.nc(fid)
dat <- read.nc(fid)
u_05 <- dat$uwnd
lat <- dat$lat
lon <- dat$lon
time_05 <- dat$time
close.nc(fid)

# U-COMPONENT 06
fname <- "~/Documents/Seaglider/Station_P/NCEP/uwnd.10m.gauss.2006.nc"
fid <- open.nc(fname)
# print.nc(fid)
dat <- read.nc(fid)
u_06 <- dat$uwnd
time_06 <- dat$time
close.nc(fid)
time_tot <- c(time_05,time_06)
time_tot_chron <- chron(dates=0,times=time_tot/24,origin=c(month=1,day=1,year=1800))

uwind_full <- abind(u_05,u_06,along=3)

# V-COMPONENT 05

fname <- "~/Documents/Seaglider/Station_P/NCEP/vwnd.10m.gauss.2005.nc"
fid <- open.nc(fname)
# print.nc(fid)
dat <- read.nc(fid)
v_05 <- dat$vwnd
close.nc(fid)

# V-COMPONENT 06
fname <- "~/Documents/Seaglider/Station_P/NCEP/vwnd.10m.gauss.2006.nc"
fid <- open.nc(fname)
# print.nc(fid)
dat <- read.nc(fid)
v_06 <- dat$vwnd
close.nc(fid)

vwind_full <- abind(v_05,v_06,along=3)

surf_winds_05_list <- list(lat=lat,lon=lon,time=time_tot_chron,u_component=uwind_full,v_component=vwind_full)
save(surf_winds_05_list,file="~/Desktop/NMML_Phase_2/R_scripts/surf_winds_05_list.Rda")

#### SURFACE CURRENTS

# AVISO-calculated UV from MADT

fname_1 <- "~/Documents/Seaglider/Station_P/AVISO_ADT/MADT_UV/dataset-duacs-dt-global-allsat-madt-uv_1447441626053.nc"
fid <- open.nc(fname_1)
dat <- read.nc(fid)
sf_1 <- att.get.nc(fid,"u","scale_factor")
u_1 <- dat$u*sf_1
v_1 <- dat$v*sf_1
lat_1 <- dat$lat;
lon_1 <- dat$lon
time_1 <- dat$time
time_1_chron <- chron(dates = time_1,times=0,origin=c(month=1,day=1,year=1950))
close.nc(fid)

fname_2 <- "~/Documents/Seaglider/Station_P/AVISO_ADT/MADT_UV/dataset-duacs-dt-global-allsat-madt-uv_1448308032615.nc"
fid <- open.nc(fname_2)
dat <- read.nc(fid)
sf_2 <- att.get.nc(fid,"u","scale_factor")
u_2 <- dat$u*sf_2
v_2 <- dat$v*sf_2
lat_2 <- dat$lat;
lon_2 <- dat$lon
time_2 <- dat$time
time_2_chron <- chron(dates = time_2,times=0,origin=c(month=1,day=1,year=1950))
close.nc(fid)

fname_3 = "~/Documents/Seaglider/Station_P/AVISO_ADT/MADT_UV/dataset-duacs-dt-global-allsat-madt-uv_1448318305976.nc"
fid <- open.nc(fname_3)
dat <- read.nc(fid)
sf_3 <- att.get.nc(fid,"u","scale_factor")
u_3 <- dat$u*sf_3
v_3 <- dat$v*sf_3
lat_3 <- dat$lat;
lon_3 <- dat$lon
time_3 <- dat$time
time_3_chron <- chron(dates = time_3,times=0,origin=c(month=1,day=1,year=1950))
close.nc(fid)

time_east = c(time_2_chron,time_3_chron)
v_east <- array(data=NA,dim=c(length(lon_3),length(lat_3),length(time_east)))
u_east <- v_east

v_east[,,1:length(time_2)] <- v_2
v_east[,,seq(from=(length(time_2)+1),to=length(time_east),by=1)] <- v_3
u_east[,,1:length(time_2)] <- u_2
u_east[,,seq(from=(length(time_2)+1),to=length(time_east),by=1)] <- u_3
 
start_05 <- match(TRUE,time_1_chron == chron(dates='01/01/2005',times=0,format=c(dates="m/d/y",times="h:m:s")))
end_06 <- match(TRUE,time_1_chron == chron(dates='12/31/2006',times=0,format=c(dates="m/d/y",times="h:m:s")))

n_east <- length(lon_2)
u_select <- abind(u_1[,,start_05:end_06],u_east[2:n_east,,start_05:end_06],along=1)
v_select <- abind(v_1[,,start_05:end_06],v_east[2:n_east,,start_05:end_06],along=1)

time_select <- time_1_chron[start_05:end_06]
lon_tot <- c(lon_1,lon_2[2:length(lon_2)])

surf_currents_05_list <- list(lat=lat_1,lon=lon_tot,time=time_select,u_component=u_select,v_component=v_select)
save(surf_currents_05_list,file="~/Desktop/NMML_Phase_2/R_scripts/surf_currents_05_list.Rda")


#### SST


fname <- "~/Documents/Seaglider/Station_P/NOAA_OI_SST_V2/X161.55.176.220.152.17.49.54.nc"
fid <- open.nc(fname)
print.nc(fid)
dat <- read.nc(fid)
sst <- dat$sst
lat <- dat$lat
lon <- dat$lon
time <- dat$time
close.nc(fid)
time_sst_1 <- chron(dates=time,times=0,origin=c(month=1,day=1,year=1800))
time_vec <- month.day.year(time,c(month=1,day=1,year=1800))
time_R <- julian(time_vec$month,time_vec$day,time_vec$year)


fname2 <- "~/Documents/Seaglider/Station_P/NOAA_OI_SST_V2/X161.55.176.220.159.16.22.4.nc"
fid <- open.nc(fname2)
print.nc(fid)
dat <- read.nc(fid)
sst2 <- dat$sst
lat2 <- dat$lat
lon2 <- dat$lon
time2 <- dat$time
time_sst_2 <- chron(dates=time2,times=0,origin=c(month=1,day=1,year=1800))
close.nc(fid)
time_vec2 <- month.day.year(time2,c(month=1,day=1,year=1800))
time_R2 <- julian(time_vec2$month,time_vec2$day,time_vec2$year)

sst_full <- abind(sst,sst2,along=3)
time_bind <- abind(time,time2,along=1)
time_full <- chron(dates=time_bind,times=0,origin=c(month=1,day=1,year=1800))
nlon <- length(lon2)
nlat <- length(lat2)
lon_mat <- lon2 %*% matrix(data=1,nrow=1,ncol=nlat)
lat_mat <- matrix(data=1,nrow=nlon,ncol=1) %*% t(lat2)
dx_denom_mat <- matrix(data=NA,nrow=nlon,ncol=nlat)
dy_denom_mat <- dx_denom_mat
a_grs80 <- 6378137
f_grs80 <- 1/298.257222100882711243
for (nli in 1:nlat){
  distance_to_next <- distGeo(c(lon_mat[3,nli],lat_mat[3,nli]),c(lon_mat[1,nli],lat_mat[1,nli]),a=a_grs80,f=f_grs80)
  dx_denom_mat[,nli] <- matrix(data=distance_to_next,nrow=nlon,ncol=1)
  if (nli == 1){
    distance_to_next <- distGeo(c(lon_mat[1,1],lat_mat[1,1]),c(lon_mat[1,3],lat_mat[1,3]),a=a_grs80,f=f_grs80)
  } else if (nli == nlat){
    distance_to_next <- distGeo(c(lon_mat[1,(nlat-2)],lat_mat[1,(nlat-2)]),c(lon_mat[1,nlat],lat_mat[1,nlat]),a=a_grs80,f=f_grs80)
  } else
    distance_to_next <- distGeo(c(lon_mat[1,nli-1],lat_mat[1,nli-1]),c(lon_mat[1,nli+1],lat_mat[1,nli+1]),a=a_grs80,f=f_grs80)
  dy_denom_mat[,nli] <- matrix(data=distance_to_next,nrow=nlon,ncol=1)
}

# compute SST gradient
dsst_dx <- sst_full*NA
dsst_dy <- sst_full*NA
for (sti in 1:length(time_full)){
  this_sst <- sst_full[,,sti]
  #### X DERIVATIVE
  # numerator
  d_dx_this_sst <- this_sst[3:nlon,] - this_sst[seq(from=1,to=(nlon-2),by=1),]
  d_dx_this_sst <- rbind(d_dx_this_sst,1*this_sst[(nlon-2),] - 4*this_sst[(nlon-1),] + 3*this_sst[nlon,])
  d_dx_this_sst <- rbind(-3*this_sst[1,] + 4*this_sst[2,] - 1*this_sst[3,],d_dx_this_sst)
  # denominator
  d_dx_this_sst <- d_dx_this_sst/dx_denom_mat
  
  #### Y DERIVATIVE 
  # numerator
  d_dy_this_sst <- this_sst[,3:nlat] - this_sst[,1:(nlat-2)]
  d_dy_this_sst <- cbind(-3*this_sst[,1] + 4*this_sst[,2] - 1*this_sst[,3],d_dy_this_sst)
  d_dy_this_sst <- cbind(d_dy_this_sst,1*this_sst[,(nlat-2)] - 4*this_sst[,(nlat-1)] - 3*this_sst[,nlat])
  # denominator
  d_dy_this_sst <- d_dy_this_sst/dy_denom_mat
  
  dsst_dx[,,sti] <- d_dx_this_sst
  dsst_dy[,,sti] <- d_dy_this_sst
  
}

sst_list_05 <- list(lat=lat,lon=lon,time=time_full,scalar_cov=sst_full,d_dx=dsst_dx,d_dy=dsst_dy)
save(sst_list_05,file="~/Desktop/NMML_Phase_2/R_scripts/sst_list_05.Rda")

#### SLP

fname <- "~/Documents/Seaglider/Station_P/NCEP/slp.2005.nc"
fid <- open.nc(fname)
print.nc(fid)
dat <- read.nc(fid)
slp_05 <- dat$slp/100
lat <- dat$lat
lon <- dat$lon
time_slp_05 <- dat$time
close.nc(fid)

fname <- "~/Documents/Seaglider/Station_P/NCEP/slp.2006.nc"
fid <- open.nc(fname)
print.nc(fid)
dat <- read.nc(fid)
slp_06 <- dat$slp/100
time_slp_06 <- dat$time
close.nc(fid)

slp_full <- abind(slp_05,slp_06,along=3)
time_full <- abind(time_slp_05,time_slp_06,along=1)
time_chron_full <- chron(dates=time_full/24,times=0,origin=c(month=1,day=1,year=1800))
slp_list_05 <- list(lat=lat,lon=lon,scalar_cov=slp_full,time=time_chron_full)
save(slp_list_05,file="~/Desktop/NMML_Phase_2/R_scripts/slp_list_05.Rda")

