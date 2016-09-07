# covar_assign_spsample_plot_check

# run draft_grid_covariate_assign_spsample first

# source("~/Desktop/NMML_Phase_2/R_scripts/basic_quiver.R")
# source("~/Desktop/NMML_Phase_2/R_scripts/rgb_val_gen.R")

time_sel <- as.POSIXlt("2006-01-01 00:00:00",tz="UTC")
id_sel <- which(grid_time == time_sel,arr.ind=TRUE)

these_rows <- which(row_header_mat[,4] == id_sel,arr.ind=TRUE)
row_header_this_time <- row_header_mat[these_rows,]
interp_slp <- master_cov_frame$slp_in_cell[these_rows]
slp_time_ind <- which(as.POSIXlt(slp_list_05$time,tz="UTC") == time_sel,arr.ind=TRUE)

dup_adj_elements <- duplicated(row_header_mat[these_rows,1])
unique_adj_ix <- which(!dup_adj_elements)
unique_interp_slp <- interp_slp[unique_adj_ix]
# unique_interp_
all_hexes <- row_header_mat[these_rows,1]
unique_hexes <- all_hexes[unique_adj_ix]
lon_west <- hex_info$Lon_Center
lon_west[lon_west < 0] <- lon_west[lon_west < 0] + 360

lon_sel <- lon_west[unique_hexes]
lat_sel <- hex_info$Lat_Center[unique_hexes]

#### plot SLP for comparison
# interp_slp_orig <- interp_slp[unique_adj_ix]
# plot(x=poly_center_lon,y=poly_center_lat,bg=rgb_val_gen(interp_slp_orig,c(970,1030))
#   ,xlim=c(190,240),ylim=c(25,65))
# t1 <- par()
# these_lims <- t1$usr
slp_lat <- slp_list_05$lat
slp_lon <- slp_list_05$lon
slp_slice <- slp_list_05$scalar_cov[,,slp_time_ind]
slp_lon_mat <- slp_lon %*% matrix(data=1,nrow=1,ncol=length(slp_lat))
slp_lat_mat <- matrix(data=1,nrow=length(slp_lon),ncol=1) %*% t(slp_lat)
# rev_inds <- seq(from=length(slp_lat),to=1,by=-1)
# slp_lat_rev <- slp_lat[rev_inds]
# 
# image.plot(x=slp_list_05$lon,y=slp_lat_rev,z=slp_list_05$scalar_cov[,rev_inds,slp_time_ind],xlim=c(min(lon_sel)-1,max(lon_sel)+1),ylim=c(min(lat_sel)-1,max(lat_sel)+1))
plot(x=as.vector(slp_lon_mat),y=as.vector(slp_lat_mat),col=rgb_val_gen(as.vector(slp_slice),c(970,1030)),cex=3,pch=15,xlim=c(min(lon_sel)-3,max(lon_sel)+3),ylim=c(min(lat_sel)-3,max(lat_sel)+3))
for (uqi in 1:length(unique_hexes)){
  lines(x=lon_sel[uqi],y=lat_sel[uqi],type="b",pch=24,lwd=1,bg=rgb_val_gen(unique_interp_slp[uqi],c(970,1030)),col="red",cex=2)
}

#### surface winds
surf_wind_lon <- surf_winds_05_list$lon
surf_wind_lat <- surf_winds_05_list$lat
surf_wind_lon_mat <- surf_wind_lon %*% matrix(data=1,nrow=1,ncol=length(surf_wind_lat))
surf_wind_lat_mat <- matrix(data=1,nrow=length(surf_wind_lon),ncol=1) %*% t(surf_wind_lat)
# slp_lon_mat <- slp_lon %*% matrix(data=1,nrow=1,ncol=length(slp_lat))
# slp_lat_mat <- matrix(data=1,nrow=length(slp_lon),ncol=1) %*% t(slp_lat)

wind_time_ind <- which(as.POSIXlt(surf_winds_05_list$time,tz="UTC") == time_sel)

## quiver on same plot as above
basic_quiver(as.vector(surf_wind_lon_mat),as.vector(surf_wind_lat_mat),as.vector(surf_winds_05_list$u_component[,,wind_time_ind])
             ,as.vector(surf_winds_05_list$v_component[,,wind_time_ind]),0.1,"red")

## surface winds, and rotation
loc <- cbind(lon_sel,lat_sel)
obj <- list(x=surf_winds_05_list$lon,y=surf_winds_05_list$lat,z=surf_winds_05_list$u_component[,,wind_time_ind])
interp_uwnd <- interp.surface(obj,loc)
obj <- list(x=surf_winds_05_list$lon,y=surf_winds_05_list$lat,z=surf_winds_05_list$v_component[,,wind_time_ind])
interp_vwnd <- interp.surface(obj,loc)
interp_spd <- sqrt(interp_uwnd^2 + interp_vwnd^2)

basic_quiver(lon_sel,lat_sel,interp_uwnd,interp_vwnd,0.1,"black")

# plot for test polygon

plot(x=as.vector(slp_lon_mat),y=as.vector(slp_lat_mat),col=rgb_val_gen(as.vector(slp_slice),c(970,1030)),cex=3,pch=15,xlim=c(min(lon_sel)-0.1,max(lon_sel)+0.1),ylim=c(min(lat_sel)-0.1,max(lat_sel)+0.1))
for (uqi in 1:length(unique_hexes)){
  lines(x=lon_sel[uqi],y=lat_sel[uqi],type="b",pch=24,lwd=1,bg=rgb_val_gen(unique_interp_slp[uqi],c(970,1030)),col="red",cex=2)
}
basic_quiver(as.vector(surf_wind_lon_mat),as.vector(surf_wind_lat_mat),as.vector(surf_winds_05_list$u_component[,,wind_time_ind])
             ,as.vector(surf_winds_05_list$v_component[,,wind_time_ind]),0.1,"red")
basic_quiver(lon_sel,lat_sel,interp_uwnd,interp_vwnd,0.1,"black")

test_hex <- unique_hexes[1]
sub_rows <- which(row_header_this_time[,1] == test_hex,arr.ind=TRUE)
sub_from <- row_header_this_time[sub_rows,1]
sub_to <- row_header_this_time[sub_rows,2]
sub_wind_to_adj <- master_cov_frame$wind_component_to_next[these_rows[sub_rows]]
sub_wind_left_to_adj <- master_cov_frame$wind_component_90_deg_left_of_next[these_rows[sub_rows]]

for (tji in 1:length(sub_to)){
  this_U_to_next <- sub_wind_to_adj[tji]
  this_U_left_of_next <- sub_wind_left_to_adj[tji]
  index_of_next <- sub_to[tji]
  text(x=lon_west[index_of_next],y=hex_info$Lat_Center[index_of_next],labels=sprintf('%1.2f %s',this_U_to_next,"m s^-1"),col="gray")
  text(x=lon_west[index_of_next],y=hex_info$Lat_Center[index_of_next]-0.07,labels=sprintf('%1.2f %s',this_U_left_of_next,"m s^-1"),col="black")
}

#### surface currents comparison
surf_currents_t_ind <- which(as.POSIXlt(surf_currents_05_list$time,tz="UTC") == time_sel)
surf_currents_lon <- surf_currents_05_list$lon
surf_currents_lat <- surf_currents_05_list$lat
u_geo_sel <- surf_currents_05_list$u_component[,,surf_currents_t_ind]
v_geo_sel <- surf_currents_05_list$v_component[,,surf_currents_t_ind]
SC_lon_mat <- surf_currents_lon %*% matrix(data=1,nrow=1,ncol=length(surf_currents_lat))
SC_lat_mat <- matrix(data=1,nrow=length(surf_currents_lon),ncol=1) %*% t(surf_currents_lat)

obj <- list(x=surf_currents_lon,y=surf_currents_lat,z=u_geo_sel)
interp_u_surf <- interp.surface(obj,loc)
obj <- list(x=surf_currents_lon,y=surf_currents_lat,z=v_geo_sel)
interp_v_surf <- interp.surface(obj,loc)

## plot
plot("n",xlim=c(min(lon_sel)-2,max(lon_sel)+2),ylim=c(min(lat_sel)-1.5,max(lat_sel)+1.5),pin=c(5,5))
basic_quiver(as.vector(SC_lon_mat),as.vector(SC_lat_mat),as.vector(u_geo_sel),as.vector(v_geo_sel)
             ,1,"black")
these_lims <- par("usr")
these_xlim <- these_lims[1:2]
these_ylim <- these_lims[3:4]
basic_quiver(lon_sel,lat_sel,interp_u_surf,interp_v_surf,1,"red")
