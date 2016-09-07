# sample hexagonal grid generation

# load pup tracks

adult_pup_frame <- read.csv(file='~/Desktop/NMML_Phase_2/Adult_pup.csv')

pup_ind <- adult_pup_frame$age == 'pup'
lat <- adult_pup_frame$SSM_lat_5
lon <- adult_pup_frame$SSM_lon_5

pup_lat <- lat[pup_ind]
pup_lon <- lon[pup_ind]

pup_lon[pup_lon < 0] <- pup_lon[pup_lon < 0] + 360
pup_site <- adult_pup_frame$site[pup_ind]
pup_bering_ind <- pup_site != 'SM'

pup_lat_b <- pup_lat[pup_bering_ind]
pup_lon_b <- pup_lon[pup_bering_ind]

bound_lat <- c(min(pup_lat_b)-2,max(pup_lat_b)+2)
bound_lon <- c(min(pup_lon_b)-2,max(pup_lon_b)+2)

input_lonlat <- list(x=pup_lon_b,y=pup_lat_b)
library(proj4)
library(sp)
library(grid)
library(hexbin)
library(geosphere)
library(maptools)
library(chron)
library(abind)
library(fields)
library(sqldf)
library(spdep)
proj_def <- "+proj=aea +lat_1=40 +lat_2=60 +lat_0=50 +lon_0=192.5 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
proj_xy <- project(input_lonlat,proj=proj_def)


grid_spacing <- 25000 # length of one side of hex [m]
grid_buffer <- 6 # number of hexes to buffer outside of actual data bounds, for imputations that exceed the latter

# create bounding box in projected coordinates
min_x <- min(proj_xy$x)-grid_buffer*grid_spacing
max_x <- max(proj_xy$x)+grid_buffer*grid_spacing
min_y <- min(proj_xy$y)-grid_buffer*grid_spacing
max_y <- max(proj_xy$y)+grid_buffer*grid_spacing
bbox_mat <- rbind(c(min_x,min_y),
                  c(min_x,max_y),
                  c(max_x,max_y),
                  c(max_x,min_y),
                  c(min_x,min_y))
bbox_poly_xy <- cbind(c(min_x,min_x,max_x,max_x,min_x),c(min_y,max_y,max_y,min_y,min_y))
bbox_poly <- Polygon(bbox_poly_xy)
hex_cent_test <- spsample(bbox_poly,n=(130*50),type="hexagonal",cellsize=40000)
hex_poly_test <- HexPoints2SpatialPolygons(hex_cent_test)
hex_poly_nb <- poly2nb(hex_poly_test)

#### STAGGERED HEX GRID IN PROJECTED COOORDINATES
grid_dx <- 3*grid_spacing
grid_dy <- grid_spacing*cos(30/180*pi)

# create center coordinates for first and second rows, which are staggered relative to one another
# these establish the unique x values for the grid
x_first_row <- seq(from = min_x,to = max_x,by = grid_dx)
x_second_row <- x_first_row+grid_dx/3+grid_dx/6

y_first_row <- matrix(data = 1,nrow = 1,ncol=length(x_first_row))*max_y
y_second_row <- matrix(data = 1,nrow =1,ncol=length(x_second_row))*(max_y-grid_dy)

# if the number of y_vals is not even, add an extra row
# this ensures that the 2-row staggered structure is repeated an integer number of times
y_vals <- seq(from = y_first_row[1],to = (min_y-grid_dy),by = -grid_dy)
n_y <- length(y_vals)
if (floor(n_y/2) != n_y/2) {
  y_vals <- append(y_vals,y_vals[n_y]-grid_dy,n_y)
  n_y <- n_y+1;
}
n_x <- length(x_first_row)

# create 2-d matrices of x, y points corresponding to centers of hex grid
x_vals_mat <- matrix(data=NA,nrow=n_y,ncol=n_x)
y_vals_mat <- x_vals_mat
x_vals_vec <- matrix(data=NA,nrow=n_y*n_x,ncol=1)
y_vals_vec <- matrix(data=NA,nrow=n_y*n_x,ncol=1)
poly_index <- x_vals_vec
poly_index_mat <- x_vals_mat
for (nyi in 1:(n_y/2)) {
  this_row_1 <- nyi*2-1
  this_row_2 <- nyi*2
  row_1_start <- (this_row_1-1)*n_x+1
  row_1_end <- (this_row_1)*n_x
  row_2_start <- (this_row_1*n_x)+1
  row_2_end <- (this_row_2*n_x)
  x_vals_vec[row_1_start:row_1_end,1] <- x_first_row
  x_vals_vec[row_2_start:row_2_end,1] <- x_second_row
  y_vals_vec[row_1_start:row_1_end,1] <- y_vals[this_row_1]
  y_vals_vec[row_2_start:row_2_end,1] <- y_vals[this_row_2]
  poly_index[row_1_start:row_1_end,1] <- row_1_start:row_1_end
  poly_index[row_2_start:row_2_end,1] <- row_2_start:row_2_end
  x_vals_mat[this_row_1,] <- x_first_row
  x_vals_mat[this_row_2,] <- x_second_row
  y_vals_mat[this_row_1,] <- y_vals[this_row_1]
  y_vals_mat[this_row_2,] <- y_vals[this_row_2]
  poly_index_mat[this_row_1,] <- row_1_start:row_1_end
  poly_index_mat[this_row_2,] <- row_2_start:row_2_end
}

# generate vertices for the edge of each hex polygon
# poly_edge_x <- matrix(rep(x_vals_vec,each=6),ncol=6,byrow=TRUE) + rep(c(1/2,1,1/2,-1/2,-1,-1/2),each=n_x*n_y)*grid_spacing
# poly_edge_y = matrix(rep(y_vals_vec,each=6),ncol=6,byrow=TRUE) + rep(c(1,0,-1,-1,0,1),each=n_x*n_y)*grid_spacing*cos(30/180*pi);

# convert these to lat/lon
# poly_edge_lat <- matrix(data=NA,nrow=n_x*n_y,ncol=6)
# poly_edge_lon <- poly_edge_lat
# for (col_i in 1:6) {
#   input_xy <- list(x=poly_edge_x[,col_i],y=poly_edge_y[,col_i])
#   this_edge_latlon = project(input_xy,proj=proj_def,inverse=TRUE)
#   poly_edge_lat[,col_i] <- this_edge_latlon$y
#   poly_edge_lon[,col_i] <- this_edge_latlon$x
# }

poly_center_x <- as.vector(t(x_vals_mat)) # in accord with poly_index vector
# using as.vector on the transpose of x_vals_mat means that hex index numbering
# starts at upper-left corner, then proceeds to east by column, repeats for next row, and so on
poly_center_y <- as.vector(t(y_vals_mat))
poly_center_xy <- list(x=poly_center_x,y=poly_center_y)
poly_center_latlon <- project(poly_center_xy,proj=proj_def,inverse=TRUE)
poly_center_lat <- poly_center_latlon$y
poly_center_lon <- poly_center_latlon$x
poly_center_lon[poly_center_lon < 0] = poly_center_lon[poly_center_lon < 0] + 360;


# Index to adjacent hexes (from interior hexes)
# row 3: col 2 to end (go to 2nd-to-last odd row)
# row 4: col 1 to end-1 (go to 2nd-to-last even row)

odd_rows <- seq(from=3,to=(n_y-3),by=2)
odd_cols <- seq(from=2,to=n_x,by=1)
even_rows <- seq(from=4,to=(n_y-2),by=2)
even_cols <- seq(from=1,to=(n_x-1),by=1)

poly_index_to_next_mat <- array(NA,dim=c(n_y,n_x,6))

####
#### ASSIGN INDEX TO NEIGHBOR
####

# ODD rows
for (ori in odd_rows){

#straight up
poly_index_to_next_mat[ori,odd_cols,1] = poly_index_mat[ori-2,odd_cols];
#up and to right
poly_index_to_next_mat[ori,odd_cols,2] = poly_index_mat[ori-1,odd_cols];
#down and to right
poly_index_to_next_mat[ori,odd_cols,3] = poly_index_mat[ori+1,odd_cols];
#straight down
poly_index_to_next_mat[ori,odd_cols,4] = poly_index_mat[ori+2,odd_cols];
#down and to left
poly_index_to_next_mat[ori,odd_cols,5] = poly_index_mat[ori+1,odd_cols-1];
#up and to left
poly_index_to_next_mat[ori,odd_cols,6] = poly_index_mat[ori-1,odd_cols-1];

}

# EVEN rows
for (eri in even_rows) {

#straight up
poly_index_to_next_mat[eri,even_cols,1] = poly_index_mat[eri-2,even_cols]
# up, right
poly_index_to_next_mat[eri,even_cols,2] = poly_index_mat[eri-1,even_cols+1]
# down, right
poly_index_to_next_mat[eri,even_cols,3] = poly_index_mat[eri+1,even_cols+1]
# straight down
poly_index_to_next_mat[eri,even_cols,4] = poly_index_mat[eri+2,even_cols]
# down, left
poly_index_to_next_mat[eri,even_cols,5] = poly_index_mat[eri+1,even_cols]
# up, left
poly_index_to_next_mat[eri,even_cols,6] = poly_index_mat[eri-1,even_cols]

}

# TAKE CARE OF EDGE POLYGONS TOO
#### left edge
left_edge_rows <- seq(from=1,to=n_y-1,by=2)
near_left_rows <- seq(from=2,to=n_y,by=2)
# up
poly_index_to_next_mat[left_edge_rows[2:length(left_edge_rows)],1,1] <- poly_index_mat[left_edge_rows[1:length(left_edge_rows)-1],1]
# up,right
poly_index_to_next_mat[left_edge_rows[2:length(left_edge_rows)],1,2] <- poly_index_mat[near_left_rows[1:length(near_left_rows)-1],1]
# down, right
poly_index_to_next_mat[left_edge_rows,1,3] <- poly_index_mat[near_left_rows,1]
# down
poly_index_to_next_mat[left_edge_rows[1:length(left_edge_rows)-1],1,4] <- poly_index_mat[left_edge_rows[2:length(left_edge_rows)],1]

#### top
## row 1
# down,right
poly_index_to_next_mat[1,1:n_x,3] <- poly_index_mat[2,1:n_x]
# down
poly_index_to_next_mat[1,1:n_x,4] <- poly_index_mat[3,1:n_x]
# down, left
poly_index_to_next_mat[1,2:n_x,5] <- poly_index_mat[2,1:(n_x-1)]

## row 2
# up,right
poly_index_to_next_mat[2,seq(from=1,to=(n_x-1),by=1),2] <- poly_index_mat[1,2:n_x]
# down,right
poly_index_to_next_mat[2,seq(from=1,to=(n_x-1),by=1),3] <- poly_index_mat[3,2:n_x]
# down
poly_index_to_next_mat[2,1:n_x,4] <- poly_index_mat[4,1:n_x]
# down,left
poly_index_to_next_mat[2,1:n_x,5] <- poly_index_mat[3,1:n_x]
# up,left
poly_index_to_next_mat[2,1:n_x,6] <- poly_index_mat[1,1:n_x]

#### right edge
# only have to worry about right-most column (rest taken care of above)
# up
poly_index_to_next_mat[seq(from=4,to=n_y,by=2),n_x,1] <- poly_index_mat[seq(from=2,to=(n_y-2),by=2),n_x]
#down
poly_index_to_next_mat[seq(from=2,to=(n_y-2),by=2),n_x,4] <- poly_index_mat[seq(from=4,to=n_y,by=2),n_x]
# down, left
poly_index_to_next_mat[seq(from=2,to=(n_y-2),by=2),n_x,5] <- poly_index_mat[seq(from=3,to=(n_y-1),by=2),n_x]
# up, left
poly_index_to_next_mat[seq(from=2,to=n_y,by=2),n_x,6] <- poly_index_mat[seq(from=1,to=(n_y-1),by=2),n_x]

#### bottom
# row n_y-1
# up
poly_index_to_next_mat[(n_y-1),1:n_x,1] <- poly_index_mat[(n_y-3),1:n_x]
# up, right
poly_index_to_next_mat[(n_y-1),1:n_x,2] <- poly_index_mat[(n_y-2),1:n_x]
# down, right
poly_index_to_next_mat[(n_y-1),1:n_x,3] <- poly_index_mat[n_y,1:n_x]
# down, left
poly_index_to_next_mat[(n_y-1),2:n_x,5] <- poly_index_mat[n_y,seq(from=1,to=(n_x-1),by=1)]
# up, left
poly_index_to_next_mat[(n_y-1),2:n_x,6] <- poly_index_mat[n_y-2,seq(from=1,to=(n_x-1),by=1)]

# row n_y
# up
poly_index_to_next_mat[n_y,1:n_x,1] <- poly_index_mat[n_y-2,1:n_x]
# up, right
poly_index_to_next_mat[n_y,seq(from=1,to=(n_x-1),by=1),2] <- poly_index_mat[n_y-1,2:n_x]
# up, left
poly_index_to_next_mat[n_y,1:n_x,6] <- poly_index_mat[n_y-1,1:n_x]

# convert to index-to-next vector, one for each of 6 topological directions
poly_index_to_next <- matrix(data=NA,nrow=n_x*n_y,ncol=6)
for (col_i in 1:6){
  # need to transpose here -- as.vector converts to vector by stacking COLUMNS
  # rather than stacking ROWS which is implicitly what's done between poly_index_mat
  # and poly_index
  poly_index_to_next[,col_i] <- as.vector(t(poly_index_to_next_mat[,,col_i]))
}

# create 2-column vector that represents all possible adjacent (from->to) linkages within the 
# hex lattice grid
n_adjacent <- sum(!is.na(poly_index_to_next))
all_adjacent_vec <- matrix(data=NA,nrow=n_adjacent,ncol=2)
adj_row <- 1
for (row_i in 1:length(poly_index)){
  this_true <- poly_index_to_next[row_i,!is.na(poly_index_to_next[row_i,])]
  n_true <- length(this_true)
  from_rows <- matrix(data=1,nrow=n_true,ncol=1)*poly_index[row_i]
  to_rows <- t(this_true)
  all_adjacent_vec[adj_row:(adj_row+n_true-1),1] <- from_rows
  all_adjacent_vec[adj_row:(adj_row+n_true-1),2] <- to_rows
  adj_row <- adj_row + n_true
}

# determine bearing to next (using all_adjacent_vec)

poly_center_lon_0 <- poly_center_lon
east_ind <- poly_center_lon_0 > 180
poly_center_lon_0[east_ind] <- poly_center_lon_0[east_ind] - 360
p1_mat <- cbind(poly_center_lon_0[all_adjacent_vec[,1]],poly_center_lat[all_adjacent_vec[,1]])
p2_mat <- cbind(poly_center_lon_0[all_adjacent_vec[,2]],poly_center_lat[all_adjacent_vec[,2]])


a_grs80 <- 6378137
f_grs80 <- 1/298.257222100882711243
bearing_to_next <- bearing(p1_mat,p2_mat,a=a_grs80,f=f_grs80)

# create edge vertices, convert to shapefile (identical to what came out of matlab -- but in projected coordinates)

poly_edge_x = poly_center_x %*% matrix(data=1,nrow=1,ncol=6) + matrix(data=1,nrow=(n_x*n_y),ncol=1) %*% c(1/2, 1, 1/2, -1/2, -1 ,-1/2)*grid_spacing;
poly_edge_y = poly_center_y %*% matrix(data=1,nrow=1,ncol=6) + matrix(data=1,nrow=(n_x*n_y),ncol=1) %*% c(1, 0, -1, -1, 0, 1)*grid_spacing*cos(30/180*pi);

# create object of class "Polygons" for each hex, store in list
poly_list <- list()
for (pti in 1:length(poly_center_x)){
  this_poly_x <- poly_edge_x[pti,]
  this_poly_y <- poly_edge_y[pti,]
  
  this_poly_coord <- cbind(this_poly_x,this_poly_y)
  poly_coord_loop <- rbind(this_poly_coord,this_poly_coord[1,])
  this_poly <- Polygon(poly_coord_loop,hole = FALSE)
  this_polys <- Polygons(list(this_poly),pti)
  poly_list[[pti]] <- this_polys
}
# create CRS projection definition
proj_CRS <- CRS(proj_def)
# create object of class SpatialPolygons from list of Polygons objects
Sr <- SpatialPolygons(poly_list,proj4string = proj_CRS)
# create data frame with hex IDs, center latitudes and longitudes
poly_data_frame <- data.frame(ID=poly_index,Lat_Center=poly_center_lat,Lon_Center=poly_center_lon)
# convert to SpatialPolygonsDataFrame object
SP <- SpatialPolygonsDataFrame(Sr,data=poly_data_frame)

# write out Polygon shapefile
writePolyShape(SP,fn='~/Desktop/NMML_Phase_2/Movement_I/hex_lattice_NFS_pup')

# create list of grid properties of polygon centers

grid_properties <- list(Cent_x=poly_center_x,Cent_y=poly_center_y,
                  Cent_lat=poly_center_lat,Cent_lon=poly_center_lon,
                  ID=poly_index,ID_to_next=poly_index_to_next,
                  ID_to_all_adjacent=all_adjacent_vec,
                  bearing_to_all_adjacent=bearing_to_next,
                  grid_spacing=grid_spacing,
                  grid_buffer=grid_buffer,
                  projection=proj_def)
# complex arguments positive EAST of NORTH
z_to_next <- complex(modulus = matrix(data=1,nrow=length(bearing_to_next),ncol=1),argument=bearing_to_next/180*pi)
z_left_of_next <- complex(modulus = matrix(data=1,nrow=length(bearing_to_next),ncol=1),argument=(bearing_to_next/180*pi - pi/2))

# grid is defined spatially - now define in time

#### here using previously-described time base (R1 winds time)
hex_time_base <- read.csv('~/Desktop/NMML_Phase_2/Movement_I/hex_time_base_2005.csv',stringsAsFactors = FALSE)
hex_time_chron <- chron(hex_time_base$Date.Base.2005,hex_time_base$Time.Base.2005,format=c(dates="y-m-d",times="h:m:s"),
                       out.format=c(dates="m/d/y",times="h:m:s"))
# lubridate
grid_properties$time_base <- hex_time_chron
save(grid_properties,file="~/Desktop/NMML_Phase_2/Movement_I/hex_grid_properties.Rda")

n_grid <- n_x*n_y

# load csv file for lat/lon data frame?
load("~/Desktop/NMML_Phase_2/Movement_I/slp_list_05.Rda")
load("~/Desktop/NMML_Phase_2/Movement_I/surf_winds_05_list.Rda")
load("~/Desktop/NMML_Phase_2/Movement_I/surf_currents_05_list.Rda")
load("~/Desktop/NMML_Phase_2/Movement_I/sst_list_05.Rda")
slp_time <- slp_list_05$time
sst_time <- sst_list_05$time
surf_winds_time <- surf_winds_05_list$time
surf_currents_time <- surf_currents_05_list$time

loc <- cbind(poly_center_lon,poly_center_lat)
db <- dbConnect(SQLite(),dbname="test_hex_covar.sqlite")
for (hxi in 1:length(hex_time_chron)){
# for (hxi in 1){
  this_time <- hex_time_chron[hxi]
# SLP time in 6-hrly intervals
  slp_t_ind <- which(slp_time == this_time)
  this_slp <- slp_list_05$slp[,,slp_t_ind]
  obj <- list(x=slp_list_05$lon,y=slp_list_05$lat,z=this_slp)
  interp_slp <- interp.surface(obj,loc)
  
  slp_to_next <- interp_slp[all_adjacent_vec[,1]] # slp covariate: simply SLP in given cell
  
# surface winds
  winds_t_ind <- which(surf_winds_time == this_time)
  obj <- list(x=surf_winds_05_list$lon,y=surf_winds_05_list$lat,z=surf_winds_05_list$uwind[,,winds_t_ind])
  interp_uwnd <- interp.surface(obj,loc)
  obj <- list(x=surf_winds_05_list$lon,y=surf_winds_05_list$lat,z=surf_winds_05_list$vwind[,,winds_t_ind])
  interp_vwnd <- interp.surface(obj,loc)
  interp_Uwnd <- sqrt(interp_uwnd^2 + interp_vwnd^2)
  
  # note complex arguments of rotation vectors positive EAST of NORTH
  # component along direction of interest after multiplying by rotation vector
  # is positive imaginary component
  interp_zwnd <- complex(real=interp_uwnd, imaginary=interp_vwnd)
  interp_zwnd_adj <- interp_zwnd[all_adjacent_vec[,1]] # expand to repeat to all adjacent hexes
  wnd_component_to_next <- Im(interp_zwnd_adj*z_to_next)
  wnd_component_left_of_next <- Im(interp_zwnd_adj*z_left_of_next)
  wnd_strength_adj <- interp_Uwnd[all_adjacent_vec[,1]]
  # CHECK U, V INTERPOLATED PROPERLY
  # CHECK COMPONENT TO NEXT IS CORRECT
  # CHECK ROTATED COMPONENT TO NEXT IS CORRECT
  
  
# surface currents
  currents_t_ind <- which(dates(this_time) == dates(surf_currents_time))
  obj <- list(x=surf_currents_05_list$lon,y=surf_currents_05_list$lat,z=surf_currents_05_list$u_surf[,,currents_t_ind])
  interp_u_surf <- interp.surface(obj,loc)
  obj <- list(x=surf_currents_05_list$lon,y=surf_currents_05_list$lat,z=surf_currents_05_list$v_surf[,,currents_t_ind])
  interp_v_surf <- interp.surface(obj,loc)
  interp_z_current <- complex(real=interp_u_surf,imaginary=interp_v_surf)
  interp_z_current_adj <- interp_z_current[all_adjacent_vec[,1]]
  current_component_to_next <- Im(interp_z_current_adj*z_to_next)
  
# sst gradient, plus SST in given cell
  sst_grad_t_ind <- which(dates(this_time) == dates(sst_time))
  obj <- list(x=sst_list_05$lon,y=sst_list_05$lat,z=sst_list_05$dsst_dx[,,sst_grad_t_ind])
  interp_x_grad <- interp.surface(obj,loc)
  obj <- list(x=sst_list_05$lon,y=sst_list_05$lat,z=sst_list_05$dsst_dy[,,sst_grad_t_ind])
  interp_y_grad <- interp.surface(obj,loc)
  interp_z_grad <- complex(real=interp_x_grad,imaginary=interp_y_grad)
  interp_z_grad_adj <- interp_z_grad[all_adjacent_vec[,1]]
  sst_grad_to_next <- Im(interp_z_grad_adj*z_to_next)
  sst_in_cell <- sst_list_05$sst[all_adjacent_vec[,1]]
  
  export_frame <- data.frame(hex_ID_from=all_adjacent_vec[,1],hex_ID_to=all_adjacent_vec[,2]
                             ,time_ID=matrix(data=hxi,nrow=n_adjacent,ncol=1),time_val=matrix(data=this_time,nrow=n_adjacent,ncol=1)
                            ,wind_component_to_next=wnd_component_to_next,wind_component_90_deg_left_of_next=wnd_component_left_of_next
                            ,current_component_to_next=current_component_to_next,sst_grad_to_next=sst_grad_to_next
                            ,sst_in_cell=sst_in_cell,slp_in_cell=slp_to_next,wind_strength_in_cell=wnd_strength_adj)
  
  dbWriteTable(conn=db,name="Covariates_NPAC_Pup_Hexes_2005",value=export_frame,row.names=FALSE,append=TRUE)  
}

dbDisconnect(db)

#### CHECK THAT OUTPUT DB HAS CORRECTLY APPENDED DATA

# check that time index is correct
# CHECK U, V INTERPOLATED PROPERLY
# CHECK COMPONENT TO NEXT IS CORRECT
# CHECK ROTATED COMPONENT TO NEXT IS CORRECT

db <- dbConnect(SQLite(),dbname="test_hex_covar.sqlite")
time_sel <- chron(dates="01-Jan-2006",times=0,format=c(dates="dd-mm-yyyy",times="h:m:s"))

these_rows <- dbGetQuery(db,statement=paste("SELECT * FROM Covariates_NPAC_Pup_Hexes_2005 WHERE time_val == ",as.numeric(time_sel),sep=""))
dbDisconnect(db)
interp_slp <- these_rows$slp_in_cell
slp_time_ind <- which(slp_time == time_sel)

dup_adj_elements <- duplicated(all_adjacent_vec[,1])
unique_adj_ix <- which(!dup_adj_elements)

#### plot SLP for comparison
interp_slp_orig <- interp_slp[unique_adj_ix]
plot(x=poly_center_lon,y=poly_center_lat,col=rgb_val_gen(interp_slp_orig,c(970,1030)),xlim=c(140,240),ylim=c(25,65))
t1 <- par()
these_lims <- t1$usr
slp_lat <- slp_list_05$lat
rev_inds <- seq(from=length(slp_lat),to=1,by=-1)
slp_lat_rev <- slp_lat[rev_inds]

image.plot(x=slp_list_05$lon,y=slp_lat_rev,z=slp_list_05$slp[,rev_inds,slp_time_ind],xlim=c(140,240),ylim=c(25,65))

#### surface winds
surf_wind_lon <- surf_winds_05_list$lon
surf_wind_lat <- surf_winds_05_list$lat
surf_wind_lon_mat <- surf_wind_lon %*% matrix(data=1,nrow=1,ncol=length(surf_wind_lat))
surf_wind_lat_mat <- matrix(data=1,nrow=length(surf_wind_lon),ncol=1) %*% t(surf_wind_lat)
# slp_lon_mat <- slp_lon %*% matrix(data=1,nrow=1,ncol=length(slp_lat))
# slp_lat_mat <- matrix(data=1,nrow=length(slp_lon),ncol=1) %*% t(slp_lat)

wind_time_ind <- which(surf_winds_time == time_sel)

## quiver on same plot as above
basic_quiver(as.vector(surf_wind_lon_mat),as.vector(surf_wind_lat_mat),as.vector(surf_winds_05_list$uwind[,,wind_time_ind])
             ,as.vector(surf_winds_05_list$vwind[,,wind_time_ind]),0.1,"red")

## surface winds, and rotation
obj <- list(x=surf_winds_05_list$lon,y=surf_winds_05_list$lat,z=surf_winds_05_list$uwind[,,wind_time_ind])
interp_uwnd <- interp.surface(obj,loc)
obj <- list(x=surf_winds_05_list$lon,y=surf_winds_05_list$lat,z=surf_winds_05_list$vwind[,,wind_time_ind])
interp_vwnd <- interp.surface(obj,loc)
interp_Uwnd <- sqrt(interp_uwnd^2 + interp_vwnd^2)

basic_quiver(poly_center_lon,poly_center_lat,interp_uwnd,interp_vwnd,0.1,"black")

# plot for test polygon
this_index <- poly_index_mat[2,1] # upper left portion of grid, on left edge
cent_lat <- poly_center_lat[which(poly_index == this_index)]
cent_lon <- poly_center_lon[which(poly_index == this_index)]

xlim <- c(-1,1)+cent_lon
ylim <- c(-1,1)*cos(cent_lat/180*pi)+cent_lat
plot("n",xlim=xlim,ylim=ylim,pin=c(5,5))
poly_plot_ind <- which(poly_center_lon > xlim[1] & poly_center_lon < xlim[2] & poly_center_lat > ylim[1] & poly_center_lat < ylim[2])

for (pti in poly_plot_ind){
  this_proj_x <- poly_edge_x[pti,]
  this_proj_y <- poly_edge_y[pti,]
  this_poly_latlon <- project(cbind(this_proj_x,this_proj_y),proj=proj_def,inverse=TRUE)
  if (this_poly_latlon[1,1] < 0){
    this_poly_latlon[,1] <- this_poly_latlon[,1]+360
  }
  
  lines(x=c(this_poly_latlon[,1],this_poly_latlon[1,1]),y=c(this_poly_latlon[,2],this_poly_latlon[1,2]),col="black")
}

basic_quiver(poly_center_lon[this_index],poly_center_lat[this_index],interp_uwnd[this_index],interp_vwnd[this_index],0.05)
text(poly_center_lon[this_index],poly_center_lat[this_index],labels=sprintf('%1.2f %s',interp_Uwnd[this_index],"m s^-1"),col="red")
these_adj_rows <- which(all_adjacent_vec[,1] == this_index)
for (tji in these_adj_rows){
  this_U_to_next <- these_rows$wind_component_to_next[tji]
  this_U_left_of_next <- these_rows$wind_component_90_deg_left_of_next[tji]
  index_of_next <- all_adjacent_vec[tji,2]
  text(x=poly_center_lon[index_of_next],y=poly_center_lat[index_of_next],labels=sprintf('%1.2f %s',this_U_to_next,"m s^-1"),col="gray")
  text(x=poly_center_lon[index_of_next],y=poly_center_lat[index_of_next]-0.07,labels=sprintf('%1.2f %s',this_U_left_of_next,"m s^-1"),col="black")
}

#### surface currents comparison
surf_currents_t_ind <- which(surf_currents_time == time_sel)
surf_currents_lon <- surf_currents_05_list$lon
surf_currents_lat <- surf_currents_05_list$lat
u_geo_sel <- surf_currents_05_list$u_surf[,,surf_currents_t_ind]
v_geo_sel <- surf_currents_05_list$v_surf[,,surf_currents_t_ind]
SC_lon_mat <- surf_currents_lon %*% matrix(data=1,nrow=1,ncol=length(surf_currents_lat))
SC_lat_mat <- matrix(data=1,nrow=length(surf_currents_lon),ncol=1) %*% t(surf_currents_lat)

obj <- list(x=surf_currents_lon,y=surf_currents_lat,z=u_geo_sel)
interp_u_surf <- interp.surface(obj,loc)
obj <- list(x=surf_currents_lon,y=surf_currents_lat,z=v_geo_sel)
interp_v_surf <- interp.surface(obj,loc)

## plot
plot("n",xlim=c(190,200),ylim=c(49,55),pin=c(5,5))
basic_quiver(as.vector(SC_lon_mat),as.vector(SC_lat_mat),as.vector(u_geo_sel),as.vector(v_geo_sel)
             ,1,"black")
these_lims <- par("usr")
these_xlim <- these_lims[1:2]
these_ylim <- these_lims[3:4]
this_ind <- poly_center_lon > these_xlim[1] & poly_center_lon < these_xlim[2] & poly_center_lat > these_ylim[1] & poly_center_lat < these_ylim[2]
basic_quiver(poly_center_lon[this_ind],poly_center_lat[this_ind],interp_u_surf[this_ind],interp_v_surf[this_ind],1,"red")

