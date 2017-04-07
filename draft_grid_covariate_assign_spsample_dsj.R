# sample hexagonal grid generation, covariate interpolation
# spsample, HexPoints2SpatialPolygon methodology

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

library(tidyverse)
library(nPacMaps)
npac_poly = npac(fortify = FALSE)

#### to be replaced by library(our_package...)
# dir_header <- "~/Desktop/NMML_Phase_2/R_Scripts/NAP_Hex_Grid_and_Covariate_Operations/"
dir_header = "CTDS-Hex-Grid-and-Covariates/"
source(paste(dir_header,"hex_grid_spsample_dsj.R",sep=""))
source(paste(dir_header,"crawl_fit_and_imp_examp.R",sep=""))
source(paste(dir_header,"scalar_cov_parse.R",sep=""))
source(paste(dir_header,"vector_cov_parse.R",sep=""))
####

message('Defining grid boundaries...')
# define the boundaries of the grid:
# load crawl-modeled 2005, 2006 tracks that should be contained within our grid

pup_frame = read_csv("Pups_2005.csv", 
                     col_types=cols(
                       id = col_character(),
                       dbid = col_integer(),
                       instrumenttype = col_character(),
                       ptt = col_integer(),
                       site = col_character(),
                       GMT = col_datetime(format = "%m/%d/%Y %H:%M:%S"),
                       depcode = col_integer(),
                       Habitat = col_character(),
                       DateDeploy = col_date(format = "%m/%d/%y"),
                       loc_class = col_factor(levels = c('3', '2', '1', '0', 'A', 'B','Z')),
                       lat = col_double(),
                       long = col_double(),
                       long2 = col_double(),
                       sex = col_factor(levels = c('F','M'))
                     )
) %>% filter(site!="SMP", loc_class!='Z') %>% select(dbid, site, sex, GMT, long, lat, loc_class) %>% droplevels()

coordinates(pup_frame) = ~ long+lat
proj4string(pup_frame) = CRS("+init=epsg:4326")
pup_frame = spTransform(pup_frame, proj4string(npac_poly))
hex_grid_list = hex_grid_spsample(pup_frame, cellsize=40000, buffer=150000)






# extract properties of this grid to be used in interpolating covariates
hex_conn <- hex_grid_list$Polygon_connection_data_fram # data frame where each row is one hex connection
hex_info <- hex_grid_list$Polygon_info_data_frame # data frame where each row is one hex
hex_nb <- hex_grid_list$hex_neighbors # output from poly2nb on hex polygons
hex_poly_test <- hex_grid_list$hex_polys_projected # polgons in projected coordinates

#### DETERMINE SPACE/TIME SLICE VISITED BY THE TRACK IN QUESTION
# pick imputations
# pick unique hexes visited by all imputations at each time point
# (make sure the imputation output is in the same projection as the grid!)

# this function fits track 355 (input actually not needed)
# and gives back the x/y coords for 30 imputations
print('Fitting crawl to track and generating imputations...')
#>>>>>># it is a PLACEHOLDER -- to be replaced by a function that can operate on a track of choice
samp <- crawl_fit_and_imp_examp(355)
sampled_tracks <- samp$imputations
samp_predTime <- samp$time
#>>>>>>#

# in any case, once we have imputations, rearrange the x/y output from these
n_imp <- length(sampled_tracks)
n_predtime <- length(samp_predTime)
imp_x_mat <- matrix(data=NA,nrow=n_predtime,ncol=n_imp)
imp_y_mat <- imp_x_mat
for (samp_i in 1:length(sampled_tracks)){
  this_xy <- sampled_tracks[[samp_i]]
  imp_x_mat[,samp_i] <- this_xy[,1]
  imp_y_mat[,samp_i] <- this_xy[,2]
}
# subsample this to hourly resolution (necessary? see below. actually, I think it is the large polygons
# that cause over() to run slowly, rather than imputation points)
# subsample_rows <- seq(from=1,to=n_predtime,by=1) # oppor
imp_x_subsample <- imp_x_mat
imp_y_subsample <- imp_y_mat
samp_time_sub <- samp_predTime


## define time grid
# time grid for covariates here chosen to be at 6-hour increments
# consistent with highest time resolution of available covariates (here, NCEP R1 winds)
# covariates with lower (daily) time resolution will have repeated values
load("surf_winds_05_list.Rda")
ncep_time <- surf_winds_05_list$time
ncep_posix <- as.POSIXlt(ncep_time,tz="UTC")

# the ncep 6-hourly time grid (should) extend before/after the chosen track
# find only the needed elements of the time grid
ncep_time_inds_int <- which(ncep_posix >= samp_predTime[1] & ncep_posix <= samp_predTime[length(samp_predTime)],arr.ind=TRUE)
#>>>>>>># should build catch in here in case grid doesn't cover the track!
# i.e., remind the user to load a different time base and try again!
#>>>>>>>#
ncep_time_inds <- c(ncep_time_inds_int[1]-1,ncep_time_inds_int,ncep_time_inds_int[length(ncep_time_inds_int)]+1)
grid_time <- ncep_posix[ncep_time_inds]
# find which hexes are visited at each of these time points, and include the hexes that neighbor these
hexes_visited_list <- vector("list",length(grid_time))
imp_x_vec <- as.vector(imp_x_subsample)
imp_y_vec <- as.vector(imp_y_subsample)
imp_xy_mat <- cbind(imp_x_vec,imp_y_vec)
imp_points <- SpatialPoints(imp_xy_mat,proj4string = CRS(proj_def))
imp_over <- over(imp_points,hex_poly_test)

print('Determining hexes visited in each time interval...')
for (nci in 1:length(ncep_time_inds)){
  this_time_ind <- ncep_time_inds[nci]
  sub_rows <- which(samp_time_sub >= ncep_posix[this_time_ind] & samp_time_sub < ncep_posix[this_time_ind+1],arr.ind=TRUE)
  if(length(sub_rows) == 0){ # no modeled animal positions in this interval
    next
  }
  sub_over <- imp_over[sub_rows]
  # find unique hexes visited by these imputations, at this grid time point
  unique_hexes_sub <- unique(sub_over)
  # find all neighbors of these unique hexes
  these_neighbors <- hex_nb[unique_hexes_sub]
  all_neighbors <- these_neighbors[[1]]
  if(length(unique_hexes_sub) > 1){
    for (tni in 2:length(these_neighbors)){
      all_neighbors <- c(all_neighbors,these_neighbors[[tni]])
    }
  }
  hexes_visited_list[[nci]] <- unique(c(unique_hexes_sub,all_neighbors))
}
# now have list of grid time points, and list of hexes on which the covariates should be defined for each
# of these time points

print('Creating row header matrix...')
# create empty 'big' data matrix, with master list of from-/to- hex connection rows for each time
row_header_mat <- matrix(data=NA,nrow=0,ncol=4)
# columns of row_header_mat: (1) hex FROM, (2) hex TO, (3) bearing of line connecting FROM-TO, (4) time ID,
# where time ID = 1 is the first ncep time greater than or equal to the first track time (defined above)
for (nci in 1:length(ncep_time_inds)){
  this_hex_visit <- hexes_visited_list[[nci]]
  # imputations generated above are exactly those to be used in model
  this_visit_header <- matrix(data=NA,nrow=0,ncol=4)
  for (thi in 1:length(this_hex_visit)){
    these_rows <- which(hex_conn$hex_from_vec == this_hex_visit[thi],arr.ind=TRUE)
    this_hex_from <- hex_conn$hex_from_vec[these_rows]
    this_hex_to <- hex_conn$hex_to_vec[these_rows]
    this_bearing_to <- hex_conn$bearing_to_next[these_rows]
    this_nci_vec <- matrix(data=nci,nrow=length(these_rows),ncol=1)
    this_header <- cbind(this_hex_from,this_hex_to,this_bearing_to,this_nci_vec)
    this_visit_header <- rbind(this_visit_header,this_header)
  }
  row_header_mat <- rbind(row_header_mat,this_visit_header)
}

# pass this row_header_mat, plus hexes_visited_list, and grid_time (=ncep_posix[ncep_time_inds]), 
# to scalar_cov_parse and vector_cov_parse
# these functions also take lists of geophysical data loaded from *.Rda files as inputs, along with desired operation on these data
# (from a limited set of operations)
# each call of *_cov_parse also adds a column to big data frame

print('Interpolating covariates...')
# prepare to interpolate/assign covariates
# SEA SURFACE TEMPERATURE
load("sst_list_05.Rda") # sst_list_05
sst_grad_col <- scalar_cov_parse(sst_list_05,"gradient_to_adjacent",row_header_mat,grid_time,hex_info)
sst_col <- scalar_cov_parse(sst_list_05,"value_in_hex",row_header_mat,grid_time,hex_info)
# SEA LEVEL PRESSURE
load("slp_list_05.Rda") # slp_list_05
slp_col <- scalar_cov_parse(slp_list_05,"value_in_hex",row_header_mat,grid_time,hex_info)
# SURFACE WINDS
load("surf_winds_05_list.Rda") 
wind_to_adjacent_col <- vector_cov_parse(surf_winds_05_list,"component_to_adjacent",row_header_mat,grid_time,hex_info)
wind_left_to_adjacent_col <- vector_cov_parse(surf_winds_05_list,"component_90_deg_left_to_adjacent",row_header_mat,grid_time,hex_info)
wind_spd_col <- vector_cov_parse(surf_winds_05_list,"magnitude_in_hex",row_header_mat,grid_time,hex_info)
# SURFACE GEOSTROPHIC CURRENTS
load("surf_currents_05_list.Rda")
ug_to_adjacent_col <- vector_cov_parse(surf_currents_05_list,"component_to_adjacent",row_header_mat,grid_time,hex_info)

print('Generating kinematic covariates...')
# KINEMATIC COVARIATES
# distance from rookery
# approximated here as dist from first point
mean_x_start <- mean(imp_x_mat[1,])
mean_y_start <- mean(imp_y_mat[1,])
lonlat_start <- project(c(mean_x_start,mean_y_start),proj=proj_def,inverse=TRUE)
a_grs80 <- 6378137
f_grs80 <- 1/298.257222100882711243
hex_lon_visited <- hex_info$Lon_Center[row_header_mat[,1]]
hex_lat_visited <- hex_info$Lat_Center[row_header_mat[,1]]
lon_start_vec <- matrix(data=lonlat_start[1],nrow=length(hex_lon_visited),ncol=1)
lat_start_vec <- matrix(data=lonlat_start[2],nrow=length(hex_lat_visited),ncol=1)
dist_from_rookery_col <- distGeo(cbind(lon_start_vec,lat_start_vec),cbind(hex_lon_visited,hex_lat_visited),a=a_grs80,f=f_grs80)

# output data frame, same format as before
master_cov_frame <- data.frame(hex_ID_from=row_header_mat[,1],hex_ID_to=row_header_mat[,2]
                               ,time_ID=row_header_mat[,4],time_val=grid_time[row_header_mat[,4]],bearing_to_next=row_header_mat[,3]
                               ,wind_component_to_next=wind_to_adjacent_col,wind_component_90_deg_left_of_next=wind_left_to_adjacent_col
                               ,current_component_to_next=ug_to_adjacent_col,sst_grad_to_next=sst_grad_col
                               ,sst_in_cell=sst_col,slp_in_cell=slp_col,wind_strength_in_cell=wind_spd_col
                               ,dist_from_rookery=dist_from_rookery_col)

print('Covariate data frame created, ready for selection and model fitting.')

# model building, fitting goes here

# secondary to above: create vignette for creating *.Rda file for each covariate?
# (secondary vignette should also have option for OpenDap rather than accessing netCDF files on computer)
# or, following DSJ advice, perhaps just avoid? (since my netCDF downloading has sometimes been in chunks over time...
# would need to re-download, plus OpenDap is incredibly slow over AFSC wi-fi)
