rgb_val_gen <- function(interp_test,clim) {
cmap_jet_data <- read.csv(file(description = "/Users/nap2/Desktop/NMML_Phase_2/R_Scripts/cmap_jet.csv"),
                          header=FALSE)

cspace <- seq(clim[1],clim[2],length=64)

rgb_val_col <- matrix(data=NA,nrow=length(interp_test),ncol=3)
cspace_length <- length(cspace)

this_interp <- approx(cspace,cmap_jet_data[,1],interp_test,
                       yleft=cmap_jet_data[1,1],yright=cmap_jet_data[cspace_length,1])
rgb_val_col[,1] <- this_interp$y
this_interp <- approx(cspace,cmap_jet_data[,2],interp_test,
                      yleft=cmap_jet_data[1,2],yright=cmap_jet_data[cspace_length,2])
 rgb_val_col[,2] <- this_interp$y
 this_interp <- approx(cspace,cmap_jet_data[,3],interp_test,
                       yleft=cmap_jet_data[1,3],yright=cmap_jet_data[cspace_length,3])
 rgb_val_col[,3] <- this_interp$y

 na_rows <- is.na(interp_test)
 rgb_val_col[na_rows,] <- 1
 output <- rgb(rgb_val_col[,1],rgb_val_col[,2],rgb_val_col[,3],max=1)
return(output)
}