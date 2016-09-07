scalar_cov_parse <- function(cov_frame,operation,row_header_mat,grid_time,hex_info){
  unique_nci <- unique(row_header_mat[,4])
  n_rows <- length(row_header_mat[,1])
  output_vec <- matrix(data=NA,nrow=n_rows,ncol=1)
  cov_time <- as.POSIXlt(cov_frame$time)
  # loop on all times in the grid
  for (nci in unique_nci){
    this_time <- grid_time[nci]
    # which rows in the output vector correspond to this time?
    these_header_rows <- which(row_header_mat[,4] == nci,arr.ind=TRUE)
    this_from <- row_header_mat[these_header_rows,1]
    from_lat <- hex_info$Lat_Center[this_from]
    from_lon <- hex_info$Lon_Center[this_from]
    from_lon[from_lon < 0] <- from_lon[from_lon < 0] + 360
    loc <- cbind(from_lon,from_lat)
    this_to <- row_header_mat[these_header_rows,2]
    this_bearing_to <- row_header_mat[these_header_rows,3]
    z_to_next <- complex(modulus = matrix(data=1,nrow=length(this_to),ncol=1),argument=this_bearing_to/180*pi)
    #### IMPORTANT -- logic of how the time index for the covariate is chosen
    cov_t_inds <- which(cov_time <= this_time,arr.ind=TRUE)
    # by finding last covariate time that is less than or equal to the grid time...
    this_t_ind <- cov_t_inds[length(cov_t_inds)]
    # this makes this function compatible with covariates with coarser time resolution.
    # if grid time is, say,
    # Nov 1 at 18:00 UTC, the daily covariate time should be Nov 1 00:00
    # this works as long as daily covariates have a day's data assigned to midnight of that day
    # which is true of NOAA OOI SST and aviso surface currents used here
    
    # select operation
    if (identical(operation,"gradient_to_adjacent")){
      # component of the gradient vector into each adjacent hex
      # interpolate gradient components to the hex centers
      this_d_dx <- cov_frame$d_dx[,,this_t_ind]
      this_d_dy <- cov_frame$d_dy[,,this_t_ind]
      obj <- list(x=cov_frame$lon,y=cov_frame$lat,z=this_d_dx)
      interp_d_dx <- interp.surface(obj,loc)
      obj$z <- this_d_dy
      interp_d_dy <- interp.surface(obj,loc)
      # convert to complex form
      interp_z_grad <- complex(real=interp_d_dx,imaginary=interp_d_dy)
      component_to_adjacent <- Im(interp_z_grad*z_to_next)
      
      output_vec[these_header_rows] <- component_to_adjacent
    } else if (identical(operation,"value_in_hex")){
      # the value of the scalar covariate in the hex
      this_cov_pane <- cov_frame$scalar_cov[,,this_t_ind]
      obj <- list(x=cov_frame$lon,y=cov_frame$lat,z=this_cov_pane)
      interp_cov <- interp.surface(obj,loc)
      output_vec[these_header_rows] <- interp_cov
      
    }
  }
  # output is a column vector of length equal to number of rows in row_header_mat
  return(output_vec)
}