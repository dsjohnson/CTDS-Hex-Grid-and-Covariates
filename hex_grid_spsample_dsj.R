hex_grid_spsample <- function(spdata, cellsize, buffer, shapefile_fname) {
  message('Generating hexagonal grid...')
  spdata@bbox <- bbox(spdata) + matrix(c(rep(-buffer,2), rep(buffer,2)), 2, 2)
  p4s = proj4string(spdata)
  
  hex_cent <- spsample(spdata, cellsize=cellsize, type="hexagonal", offset=rep(0.5,2))
  hex_poly <- HexPoints2SpatialPolygons(hex_cent)
  hex_nb = dnearneigh(hex_cent, cellsize-0.01*cellsize, cellsize+0.01*cellsize)
  
  message('Calculating grid topology...')
  hex_to <- unlist(hex_nb)
  hex_from = rep(1:length(hex_cent), sapply(hex_nb, length)) 
  to_next =  coordinates(hex_cent[hex_to,]) - coordinates(hex_cent[hex_from,])
  bearing_to_next = (atan2(to_next[,1], to_next[,2])*180/pi) %>% ifelse(.<0, 360+., .) 
  
  # hex_center_latlon <- spTransform(hex_cent, CRS("+init=epsg:4326")) %>% coordinates()
  # hex_center_latlon_from <- hex_center_latlon[hex_from_vec,]
  # hex_center_latlon_to <- hex_center_latlon[hex_to_vec,]
  # bearing_to_next <- bearing(hex_center_latlon_from,hex_center_latlon_to)
  # dist_to_next <- distGeo(hex_center_latlon_from,hex_center_latlon_to)
  
  message('Outputing data...')
  if(!missing(shapefile_fname)){
    # message('Writing shapefile...')
    # poly_ID_char = sapply(slot(hex_poly, "polygons"), function(x) slot(x, "ID"))
    # poly_index = c(1:length(hex_poly))
    # poly_data_frame <- data.frame(row.names=poly_ID_char,Lat_Center=hex_center_latlon[,2],Lon_Center=hex_center_latlon[,1])
    # SP <- SpatialPolygonsDataFrame(hex_poly,data=poly_data_frame)
    # writePolyShape(SP,fn=shapefile_fname)
  }
  # poly_connect_data_frame <- data.frame(hex_from_vec,hex_to_vec,bearing_to_next,dist_to_next)
  
  message('Writing output...')
  output <- list(
    # Polygon_info_data_frame=poly_data_frame,
    # Polygon_connection_data_frame=poly_connect_data_frame,
    # hex_neighbors=hex_poly_nb,
    # hex_polys_projected=hex_poly,
    # hex_cent_xy_projected=hex_coord
    poly=hex_poly,
    neighbor_df = data.frame(hex_from=hex_from, hex_to=hex_to, bearing=bearing_to_next),
    hex_centroids = hex_cent
  )
}