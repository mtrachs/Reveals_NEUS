library(rstan)
library(RColorBrewer)
library(fields)
#-----------------------------------------------------------------------------------------------
setwd('~/workflow_stepps_calibration/prediction/')
help.fun.loc <- 'utils/'
data.loc <- 'data/'
plot.loc <- 'plots/'
#-----------------------------------------------------------------------------------------------
taxa <- sort(c('Maple','Birch','Hickory','Chestnut','Other conifer','Other hardwood',
               'Beech','Ash','Tamarack','Spruce','Pine','Poplar','Oak','Hemlock','Elm'))


#transfrom pollen coordinates to us coordinates
source('R/prepare_data_settlement_era.R')
sputm <- SpatialPoints(pol_table@coords, proj4string=CRS("+init=epsg:4326"))
spgeo <- spTransform(sputm, CRS("+init=epsg:3175"))
pollen_coord_us <- spgeo@coords
load('~/workflow_stepps_calibration/vegetation/data/veg_data_15_taxa_754_cells_260_knots.rdata')


#------------------------------------------------------------------------------------------------------------------
setwd('~/Reveals_NEUS/')
help.fun.loc <- 'utils/'
data.loc <- 'data/'
plot.loc <- 'figures/'






pdf(paste(plot.loc,'settlement_vegetation.pdf',sep=''),width = 5,height = 4)
sapply(taxa, function(x){
  
  
  breaks <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,1)
  east <- sort(unique(coord.agg.final$east))
  north <- sort(unique(coord.agg.final$north))
  
  breaks1 <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7)
  breaks2 <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,NA)
  #----------------------------------------------------------------------------------------------------------------------
  veg_agg_matrix <- as.matrix(veg_agg)
  
  categories <- cut(veg_agg_matrix,breaks)
  cat.sort <- sort(levels(categories))
  categories <- matrix(categories,ncol=length(taxa))
  colours <- rev(brewer.pal(10,'RdYlBu'))
  colours.plot <- matrix(ncol = ncol(categories),nrow=nrow(categories))
  
  for(i in cat.sort) {
    colours.plot[categories==i] <- colours[cat.sort==i]
  }
  colnames(colours.plot) <- taxa
  
  #sapply(taxa, function(x){
  image.plot(east, north,z = matrix(ncol= length(north),rep(-1,length(east)*length(north))),zlim=c(0,0.7),col = colours,
             breaks= breaks1,lab.breaks = breaks2,main = paste(x, 'Paciorek et al. (2016)'),cex.axis = 0.8)
  points(coord.agg.final,col=colours.plot[,x],pch = 15)
  points(pollen_coord_us,col=1,pch = 16,cex = 0.75)
})
dev.off()

