vector_cov_parse <- function(cov_frame,operation,row_header_mat,grid_time,hex_info){
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
    z_90_deg_left_of_next <- z_to_next*complex(modulus=matrix(data=1,nrow=length(this_to),ncol=1),argument=(-pi/2))
    #### IMPORTANT -- logic of how the time index for the covariate is chosen
    cov_t_inds <- which(cov_time <= this_time,arr.ind=TRUE)
    # by finding last covariate time that is less than or equal to the grid time...
    this_t_ind <- cov_t_inds[length(cov_t_inds)]
    # this makes this function compatible with covariates with coarser time resolution.
    # if grid time is, say,
    # Nov 1 at 18:00 UTC, the daily covariate time should be Nov 1 00:00
    # this works as long as daily covariates have a day's data assigned to midnight of that day
    # which is true of NOAA OOI SST and aviso surface currents used here
    
    # interpolate components to the hex centers
    this_u <- cov_frame$u_component[,,this_t_ind]
    this_v <- cov_frame$v_component[,,this_t_ind]
    obj <- list(x=cov_frame$lon,y=cov_frame$lat,z=this_u)
    interp_u <- interp.surface(obj,loc)
    obj$z <- this_v
    interp_v <- interp.surface(obj,loc)
    # convert to complex form
    interp_z <- complex(real=interp_u,imaginary=interp_v)
    interp_spd <- sqrt(interp_u^2 + interp_v^2)
    # select operation
    if (identical(operation,"component_to_adjacent")){
      # component of the vector into each adjacent hex
      component_to_adjacent <- Im(interp_z*z_to_next)
      
      output_vec[these_header_rows] <- component_to_adjacent
    } else if (identical(operation,"component_90_deg_left_to_adjacent")){
      component_left <- Im(interp_z*z_90_deg_left_of_next)
      
      output_vec[these_header_rows] <- component_left
    } else if (identical(operation,"magnitude_in_hex")){
      # 2-norm of interpolated vector at each hex
      output_vec[these_header_rows] <- interp_spd
      
    }
  }
  # output is a column vector of length equal to number of rows in row_header_mat
  return(output_vec)
}