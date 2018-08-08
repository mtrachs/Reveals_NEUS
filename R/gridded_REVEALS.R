#######################################################################################################################
# in this script we want to numerically compare REVEALS with the PAciorek et al dtaset
# REVEALS assumes homogenous vegetation around depositional site. 
# we therefore assign all of a 1 by 1 deg grid cell with the same vegetation and later on compare this to the paciorek et al. 
# dataset
#########################################################################################################################
repo_location <- '~/workflow_stepps_calibration/'
setwd(paste(repo_location,'prediction/',sep=''))
help.fun.loc <- 'utils/'
data.loc <- 'output_nb/'
plot.loc <- 'plots_nb/'#'plots_w_full/'
#-----------------------------------------------------------------------------------------------
us.shp <- readShapeLines('data/map_data/us_alb.shp',proj4string=CRS('+init=epsg:3175'))
#-----------------------------------------------------------------------------------------------

#-------------------------------
#first load REVEALS data and its coordinates
#-------------------------------
reveals_reconstruction <- readRDS('~/Reveals_NEUS/output/veg_pred_ppe_literature_fallspeed_measured_composition_no_larch_lake_size_correct_dwm_gpm neutral_50000_regional_vegetation_aggregated.RDS')
reveals_coord <- readRDS('~/Reveals_NEUS/output/veg_pred_ppe_literature_fallspeed_measured_composition_no_larch_lake_size_correct_dwm_gpm neutral_50000_regional_vegetation_aggregated_coord.RDS')

reveals_median <- reveals_reconstruction[,grep('median',colnames(reveals_reconstruction))] 
reveals_mean <- reveals_reconstruction[,grep('mean',colnames(reveals_reconstruction))] 


#Transform reveals coordinates in to US coordinates
sputm <- SpatialPoints(reveals_coord[,c('x','y')], proj4string=CRS("+init=epsg:4326"))
spgeo <- spTransform(sputm, CRS("+init=epsg:3175"))
reveals_coord_us <- spgeo@coords
#------------------------------------------------------------------------------------------------------------------------
#########################################################################################################################
bounding_coords <- matrix(ncol=4, nrow = nrow(reveals_coord))
bounding_coords[,1] <- reveals_coord[,'x']-0.5
bounding_coords[,2] <- reveals_coord[,'x']+0.5
bounding_coords[,3] <- reveals_coord[,'y']-0.5
bounding_coords[,4] <- reveals_coord[,'y']+0.5

colnames(bounding_coords) <- c('lon.min','lon.max','lat.min','lat.max')

########################################################################################################################
#reverse transform US coordinates into degrees
#load vegetation coordinates
load(paste(repo_location,'calibration/data/elicitation_neus_certainty_median_79_sites_only_abies_new_species.RData',sep=''))
taxa <- colnames(y)

spgeo <- SpatialPoints(veg_coords, proj4string=CRS("+init=epsg:3175"))
sputm <- spTransform(spgeo, CRS("+init=epsg:4326"))
veg_coords_utm  <- sputm@coords #transformation is successful
colnames(veg_coords_utm) <- c('lon','lat')

coords_assign <-  matrix(ncol = 1,nrow =nrow(veg_coords))

for(x in 1:nrow(bounding_coords)){
  coords_assign[((veg_coords_utm[,'lon'] > bounding_coords[x,"lon.min"]) & (veg_coords_utm[,'lon'] < bounding_coords[x,"lon.max"]) &
  (veg_coords_utm[,'lat'] > bounding_coords[x,"lat.min"]) & (veg_coords_utm[,'lat'] < bounding_coords[x,"lat.max"]))] <- x
}


#----------------------------------------------------------------------------------------------------------------------
#lets plot reveals

taxa_plot <- c('Beech','Hemlock','Oak','Spruce')
colours <- rev(brewer.pal(10,'RdYlBu'))
east <- sort(unique(veg_coords[,'meters.east']))
north <- sort(unique(veg_coords[,'meters.north']))
breaks <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,1)
breaks1 <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7)
breaks2 <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,NA)


taxa[grep('Other',taxa)] <- c('Other_conifer','Other_hardwood')

sapply(taxa[c(1:3,5:12)],function(y){

  categories <- cut(reveals_mean[,grep(y,colnames(reveals_mean))]/100,breaks,include.lowest = TRUE,labels = FALSE)
  colours.plot <- colours[categories]

  colours.plot.final <- matrix(ncol=1, nrow=length(coords_assign))

  for(i in 1:length(colours.plot)){
    colours.plot.final[coords_assign==i] <- colours.plot[i]
  }  

  colours.plot.final[is.na(colours.plot.final)] <- 'gray'


  plot(us.shp,xlim=range(east),ylim=range(north),main =y,cex.main = 2)  
  points(veg_coords,col=colours.plot.final,pch = 15,cex=0.5)
  #points(pollen_coord_us,col=1,pch = 16,cex = 0.75)
  plot(us.shp,add=TRUE)

})

#########################################################################################################################
# goal is to make a data produc that we can use to estimate dissimilairites
#########################################################################################################################

reveals_gridded <- matrix(ncol=11,nrow=nrow(veg_coords))

for(i in 1:nrow(reveals_mean)){
  for(ii in 1:ncol(reveals_mean)){
    reveals_gridded[coords_assign==i,ii] <-  reveals_mean[i,ii]/100
  }
}

colnames(reveals_gridded) <- colnames(reveals_mean)

saveRDS(reveals_gridded,'~/Reveals_NEUS/data/reveals_mean_gridded.RDS')
