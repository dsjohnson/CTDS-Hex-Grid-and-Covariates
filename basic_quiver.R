basic_quiver <- function(x,y,u,v,scale_fac,col_option) {
  # if scale_fac = 1, vector of magnitude 1 has length in inches equal to 1 unit on x axis in plot
  t1 <- par()
  these_lims <- t1$usr
  these_xlim <- these_lims[1:2]
  these_ylim <- these_lims[3:4]
  this_extent <- t1$pin
  x_in_per_unit <- this_extent[1]/diff(these_xlim)
  y_in_per_unit <- this_extent[2]/diff(these_ylim)
  x_scale <- scale_fac
  y_scale <- scale_fac*(x_in_per_unit/y_in_per_unit)
  n_pts <- length(x)
  for (nti in 1:n_pts){
    this_x <- c(0,u[nti])*x_scale + x[nti]
    this_y <- c(0,v[nti])*y_scale + y[nti]
    lines(this_x,this_y,col=col_option)
    lines(x[nti],y[nti],pch=16)
  }
  output <- TRUE
  return(output)
}