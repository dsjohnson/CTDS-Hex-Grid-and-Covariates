hex_grid_spsample <- function(bbox_poly,cellsize_input,n_input,shapefile_fname) {
  
  print('Generating hexagonal grid...')
  hex_cent_test <- spsample(bbox_poly,n=n_input,type="hexagonal",cellsize=cellsize_input)
  hex_poly_test <- HexPoints2SpatialPolygons(hex_cent_test)
  hex_poly_nb <- poly2nb(hex_poly_test)
  
  print('Calculating grid topology...')
  hex_from_vec <- matrix(data=NA,nrow=length(hex_poly_test),ncol=1)
  hex_to_vec <- hex_from_vec
  row_i <- 1
  for (hpi in 1:length(hex_poly_test)) {
    this_nb <- hex_poly_nb[hpi][[1]]
    length_nb <- length(this_nb)
    these_rows <- seq(from=row_i,to=(row_i+length_nb-1),by=1)
    hex_to_vec[these_rows] <- this_nb
    hex_from_vec[these_rows] <- hpi
    row_i <- row_i + length_nb
  }
  
  hex_coord_test <- coordinates(hex_cent_test)
  hex_center_latlon <- project(hex_coord_test,proj=proj_def,inverse=TRUE)
  
  hex_center_latlon_from <- hex_center_latlon[hex_from_vec,]
  hex_center_latlon_to <- hex_center_latlon[hex_to_vec,]
  
  a_grs80 <- 6378137
  f_grs80 <- 1/298.257222100882711243
  bearing_to_next <- bearing(hex_center_latlon_from,hex_center_latlon_to,a=a_grs80,f=f_grs80)
  dist_to_next <- distGeo(hex_center_latlon_from,hex_center_latlon_to,a_grs80,f=f_grs80)
  
  print('Generating shapefile...')
  # write out shapefile
  # create CRS projection definition
  poly_ID_char <- sapply(slot(hex_poly_test, "polygons"), function(x) slot(x, "ID"))
  proj_CRS <- CRS(proj_def)
  hex_poly_test@proj4string <- proj_CRS
  poly_index <- seq(from=1,to=length(hex_poly_test),by=1)
  # create data frame with hex IDs, center latitudes and longitudes
  poly_data_frame <- data.frame(row.names=poly_ID_char,Lat_Center=hex_center_latlon[,2],Lon_Center=hex_center_latlon[,1])
  # convert to SpatialPolygonsDataFrame object
  SP <- SpatialPolygonsDataFrame(hex_poly_test,data=poly_data_frame)
  
  # write out Polygon shapefile
  writePolyShape(SP,fn=shapefile_fname)
  poly_connect_data_frame <- data.frame(hex_from_vec,hex_to_vec,bearing_to_next,dist_to_next)
  print('Writing output...')
  output <- list(Polygon_info_data_frame=poly_data_frame,Polygon_connection_data_fram=poly_connect_data_frame,
                 hex_neighbors=hex_poly_nb,hex_polys_projected=hex_poly_test,hex_cent_xy_projected=hex_coord_test)
}