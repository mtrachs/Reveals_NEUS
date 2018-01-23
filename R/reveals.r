library(raster)
library(sp)
library(DISQOVER)
library(reshape2)
library(ggplot2)
library(dplyr)
library(neotoma)
library(maps)
library(RColorBrewer)
library(fields)

devtools::install_github("PalEON-Project/stepps-cal")
library(stepps)

# # load package
# load('DISQOVER/R/sysdata.rda')

setwd('~/Reveals_NEUS/')

#what is this?
#source('r/utils/main.R')

ena <- TRUE

##################################################################################################################################################
## pull pollen data for NA
##################################################################################################################################################

# sum samples within a time bin for each site 
library(plyr)
#pollen_bin = ddply(pollen_trans, c('dataset', 'lat', 'long', 'site.name', 'slice_bin'),  
#                   function(x) colSums(x[tolower(taxa)]))

# make grid for NA (or ENA)
source('R/make_grid.R')
#-------------------------------------------------------------------------------------------------------------------
#find coordinates of NEUS in Degrees
#-------------------------------------------------------------------------------------------------------------------
# load vegetation data prepared for STEPPS
#source('~/workflow_stepps_prediction/prediction/R/pollen_data_settlement_era.R')
stepps_input <- load('~/workflow_stepps_prediction/prediction/data/elicitation_neus_certainty_median.RData')
pollen_coords <- as.data.frame(pollen_coords)
pollen_bin <- data.frame(pollen=y,lon = pollen_coords$lon,lat = pollen_coords$lat)



sapply(c('larch','no_larch'),function(larix){

taxa <- colnames(y)
taxa[taxa%in%c('Tamarack','Other conifer','Other hardwood')] <-c('Fir','Alder','Larch')
taxa <- taxa[(taxa%in%'Chestnut')==FALSE]
if(larix=='no_larch') taxa <- taxa[(taxa%in%'Larch')==FALSE]

grid <- make_grid(pollen_bin, coord_fun = ~ lon + lat, projection = '+init=epsg:4326', resolution = 1)

cell_id <- extract(grid, pollen_bin[,c('lon', 'lat')])

grid_coords <- rasterToPoints(grid)
# 
#-----------------------------------------------------------------------------------------------------------------
#load available lake sizes
plot_index <- read.table('~/workflow_stepps_calibration/expert_elicitation/data/plot_index.txt',header=TRUE)
lake_names <- read.csv('~/workflow_stepps_calibration/expert_elicitation/data/site_names.csv')
lake_names <- lake_names[as.matrix(plot_index),]

lake_size <- read.csv('~/Reveals_NEUS/reveals-na/data/area_lakes_1.1.csv')
lake_index_NEUS <- lake_size$site.name%in%lake_names$site.names
lake_sizes_NEUS <- lake_size[lake_index_NEUS,]

map('state')
points(lake_sizes_NEUS$long,lake_sizes_NEUS$lat,pch = 16,cex = 0.5)
# one site in Minnesota and one in South Carolina 

lake_sizes_NEUS <- lake_sizes_NEUS[((lake_sizes_NEUS$long>(-90))&(lake_sizes_NEUS$lat>38)),]
lake_size <- lake_sizes_NEUS$area


##-----------------------------------------------------------------------------------------------------------------
#function to generate random lake size
random_lake_size <- function(lake_size){
  library(truncnorm)
  lake_radius = rtruncnorm(1,
                           a=5,
                           mean=mean(lake_size, na.rm=TRUE),
                           sd=sd(lake_size, na.rm=TRUE))
  
  return(lake_radius)
}

#---------------------------------------------------------------------------------------------------------------
#produce random lake sizes
lakes_no_size <- length(lake_names$site.names) -length(lake_size)
rand_la_size <- replicate(lakes_no_size,random_lake_size(lake_size = lake_size))


sort_lake_names <-  lake_names[order(lake_names$site.names),]
sort_lake_sizes <- lake_sizes_NEUS[order(lake_sizes_NEUS$site.name),]

#a few sites have two pollen samples
total_lake_sizes <- vector(length = length(lake_names$site.names))

for (i in 1:length(sort_lake_sizes$site.name)) {
  total_lake_sizes[as.character(sort_lake_names$site.names)==as.character(sort_lake_sizes$site.name[i])]   <- sort_lake_sizes$area[i] 
}

total_lake_sizes[total_lake_sizes==0] <- replicate(sum(total_lake_sizes==0),random_lake_size(lake_size = lake_size))


# in HA; convert to radius in m
# pi * r * r
lake_size = sqrt(total_lake_sizes*0.01 / pi)*1000


#------------------------------------------------------------------------------------------------------------------
#have to find the data that we want to use at the end of the day!!!
#data used for calibration plus data used for 
#function that evaluates sites
help.fun.loc <- '~/workflow_stepps_calibration/calibration/calibration_helper_funs/'
source(paste(help.fun.loc,'evaluate_elicitation_certainty.R',sep=''))
lake_radius <- lake_size[site.sample$site.index]
lake_radius[lake_radius<10] <-10


#look at this...
pollen_bin <- data.frame(lake_radius, pollen_bin)


##################################################################################################################################################
## read in and prep ppes and svs
##################################################################################################################################################
ppes = readRDS('data/PPEs_agg.RDS')
ppes = ppes[which(ppes$taxon %in% taxa),] #no PPE for Chestnut

ppes[which(ppes$error == 0),'error'] = 0.01

###############################################################################################################################
## read in SVs
svs = read.csv('reveals-na/data/svs_LC6K.csv', sep=',', header=TRUE, stringsAsFactors=FALSE)
#svs$taxon[svs$taxon=='Fir'] <- 'Abies'
svs[which(is.na(svs$sv)), 'sv'] = 0.01
svs_agg = aggregate(sv ~ taxon, svs, median)

svs_agg = svs_agg[svs_agg$taxon %in% taxa, ]

## construct param input data frame




sapply(c('literature','stepps'),function(ppe_estimated) {
  sapply(c('measured','constant'),function(fall_speed) {
    
  


    if(ppe_estimated == 'literature') {
  
      if(fall_speed=='measured') fallspeed=svs_agg$sv[svs_agg$taxon%in%ppes$taxon]
      if(fall_speed=='constant') fallspeed=rep(mean(svs_agg$sv[svs_agg$taxon%in%ppes$taxon]),length(ppes$taxon))
      reveals_inputs = data.frame(species=ppes$taxon,
                                fallspeed=fallspeed,
                                PPEs=ppes$ppe,
                                PPE.errors=ppes$error)
      reveals_inputs$species <- as.character(reveals_inputs$species)
      if(larix =='larch') reveals_inputs$species[reveals_inputs$species%in%c('Fir','Alder','Larch')] <- c('Other_hardwood','Other_conifer','Tamarack')
      if(larix =='no_larch') reveals_inputs$species[reveals_inputs$species%in%c('Fir','Alder')] <- c('Other_hardwood','Other_conifer')
      reveals_inputs <- reveals_inputs[order(reveals_inputs$species),]
      rownames(reveals_inputs) = NULL
      reveals_inputs <- reveals_inputs[order(reveals_inputs$species),]
    }

    if(ppe_estimated == 'stepps') {
      ppe_stepps <- as.data.frame(readRDS('~/Reveals_NEUS/data/ppe_stepps.RDS'))
      ppe_stepps <- ppe_stepps[rownames(ppe_stepps)!='Chestnut',]
      if(larix=='no_larch') ppe_stepps <- ppe_stepps[rownames(ppe_stepps)!='Tamarack',]
      if(larix=='larch') svs_agg$taxon[svs_agg$taxon%in%c('Alder','Fir','Larch')] <- c('Other hardwood','Other conifer','Tamarack')
      if(larix=='no_larch') svs_agg$taxon[svs_agg$taxon%in%c('Alder','Fir')] <- c('Other hardwood','Other conifer')
      if(fall_speed=='measured') fallspeed=svs_agg$sv[order(svs_agg$taxon)]
      if(fall_speed=='constant') fallspeed=rep(mean(svs_agg$sv[order(svs_agg$taxon)]),length(taxa))
      reveals_inputs = data.frame(species=taxa,
                              fallspeed=fallspeed,
                              PPEs=ppe_stepps$ppe_STEPPS_mean_stand,
                              PPE.errors=ppe_stepps$ppe_STEPPS_sd_stand)
    }
    write.csv(reveals_inputs, "data/reveals_input_params.csv", row.names=FALSE)

##################################################################################################################################################
## run reveals
##################################################################################################################################################
    ids = unique(pollen_bin$dataset)
    pol_dat = pollen_bin


    if(larix=='larch') veg_pred <- matrix(ncol = ncol(pol_dat)-2,nrow = 1) 
    if(larix=='no_larch') veg_pred <- matrix(ncol = ncol(pol_dat)-3,nrow = 1) 
    colnames(veg_pred) <- NULL


    for (i in 1:nrow(pol_dat)){
  
      print(i)
      id = ids[i]
  
      counts_site = pol_dat
      colnames(counts_site)[colnames(counts_site)%in%c("pollen.Other.conifer" ,"pollen.Other.hardwood")] <- 
        c("pollen.Other_conifer" , "pollen.Other_hardwood")
      counts_site <- cbind(matrix(ncol=1, nrow(counts_site)),counts_site)
      colnames(counts_site)[1] <- 'ages'
   
      basin_radius = lake_radius
  
  
      coords_site = pollen_coords
      rownames(coords_site) = NULL
  
  
  
      if(larix=='larch') write.csv((counts_site[i,(colnames(counts_site)%in%c('lake_radius','pollen.Chestnut','lon','lat'))==FALSE]), 'data/reveals_input.csv', row.names=FALSE)
      if(larix=='no_larch') write.csv((counts_site[i,(colnames(counts_site)%in%c('lake_radius','pollen.Chestnut','pollen.Tamarack','lon','lat'))==FALSE]), 'data/reveals_input.csv', row.names=FALSE)
  
  # cycle through and estimate background veg
  # with english csv files
    a <- REVEALSinR(pollenFile = "data/reveals_input.csv",
                  pf         = "data/reveals_input_params.csv",
                  filetype   = "csv",
                  dwm        = "gpm neutral",
                  tBasin     = "lake",
                  dBasin     = 2*round(basin_radius[i]), # diameter!
                  regionCutoff = 100000,
                  repeats      = 1000)
  
    veg = melt(a, id.vars=c('Pollen.file', 'Parameter.file', 'Distance.weighting', 'Basin.type', 'ages'))
  
    veg$taxon = unlist(lapply(veg$variable, function(x) unlist(strsplit(as.character(x), ".", fixed=TRUE))[1]))
    veg$type  = unlist(lapply(veg$variable, function(x) unlist(strsplit(as.character(x), ".", fixed=TRUE))[2]))
  
    #veg_cast = dcast(veg, taxon ~ ages + type)
    veg_mean <- veg$value[grep('mean',veg$variable)]
  
    veg_inter <- data.frame(coords_site[i,], t(veg_mean))
    colnames(veg_inter) <- NULL
  
    veg_pred = rbind(veg_pred,as.matrix(veg_inter))
    }

    colnames(veg_pred) <- c('lon','lat',taxa[taxa!='Chestnut'])
    veg_pred <- veg_pred[-1,] 

    saveRDS(veg_pred, 'data/veg_pred.RDS')
  
    veg_pred <- cbind(veg_pred,cell_id)
    veg_pred_agg <- aggregate(veg_pred[,-c(1:2,ncol(veg_pred))],by = list(veg_pred[,ncol(veg_pred)]),FUN = mean)

    coord_veg_pred <- grid_coords[grid_coords[,'layer']%in%veg_pred_agg$Group.1,]

    veg_pred_agg <- matrix(ncol = ncol(veg_pred_agg),unlist(veg_pred_agg))
    veg_pred_agg <- veg_pred_agg[,-1]
    #----------------------------------------------------------------------------------------------------------------------
    #look at figures
    breaks <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,1)
    categories <- cut(veg_pred_agg/100,breaks,include.lowest=TRUE)
    cat.sort <- sort(levels(categories))
    categories <- matrix(categories,ncol=ncol(veg_pred_agg))
    colours <- rev(brewer.pal(10,'RdYlBu'))
    colours.plot <- matrix(ncol = ncol(categories),nrow=nrow(categories))

    for(i in cat.sort) {
      colours.plot[categories==i] <- colours[cat.sort==i]
    }
    colnames(colours.plot) <- colnames(veg_pred)[-c(1:2,ncol(veg_pred))]

    breaks1 <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7)
    breaks2 <- c(0,0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,NA)

    east <- sort(unique(grid_coords[,'x']))
    north <- sort(unique(grid_coords[,'y']))

    pdf(paste('~/Reveals_NEUS/figures/ppe_comparison/ppe_',ppe_estimated,'_fallspeed_',fall_speed,'_composition_',larix,'.pdf',sep=''),height = 3.9,width =5)
    for (i in 1:ncol(colours.plot)){
    image.plot(east, north,z = matrix(ncol= length(north),rep(-1,length(east)*length(north))),zlim=c(0,0.7),col = colours,
             breaks= breaks1,lab.breaks = breaks2,main = colnames(colours.plot)[i],cex.axis = 0.8)
    map('state',add=TRUE)
    points(coord_veg_pred,col=colours.plot[,i],pch = 15,cex = 3)
    }
    dev.off()
  })
})
})







##################################################################################################################################################
## process output
##################################################################################################################################################
library(maptools)
library(fields)

# how to get the coords for grid cells
coords = xyFromCell(grid, veg_pred$cell_id)
veg_pred = cbind(coords, veg_pred)

na_shp <- readShapeLines('data/map_data/na/NA_States_Provinces_Albers.shp',
                         proj4string=CRS('+init=epsg:3175'))
na_shp@data$id <- rownames(na_shp@data)
# na_shp_ll <- spTransform(na_shp, CRS("+proj=longlat +datum=WGS84"))
na_shp_ll <- spTransform(na_shp, CRS("+init=epsg:4326"))
na_fort <- fortify(na_shp_ll, region='id') 

veg_grid = aggregate(mediansim ~ taxon + ages + cell_id + x+ y, veg_pred, sum)

veg_cast = dcast(veg_grid, cell_id + x + y + ages ~ taxon, value.var='mediansim')
veg_cast[,5:ncol(veg_cast)] = t(apply(veg_cast[,5:ncol(veg_cast)], 1, function(x) x/sum(x)))

veg_grid = melt(veg_cast, id.vars=c('cell_id', 'x', 'y', 'ages'))

saveRDS(veg_grid, 'data/cache/veg_grid.RDS')

latlimits  <- c(40, 60) 
longlimits <- c(-100, -60) 

# c(-100, 40, -60, 60)

# # make time bins
# # start with: 0.1-0.35k BP, 0.35-0.7k BP, 5.5-6.2k BP
# breaks = c(0.1, 0.35, 0.7, 2.7, 3.2, 5.7, 6.2)
# slice_bins = seq(1, 6)
# slice_labels = c(50, 200, 500, 1500, 3000, 6000)

nslices = length(slice_labels)

pdf('figures/reveals_taxon_slices.pdf')
for (i in 1:nslices){
  
  age = slice_labels[i]
  veg_sub = veg_grid[which(veg_grid$ages==age),]
  
  p <- ggplot(data=veg_sub)
  p <- p + geom_tile(aes(x=x, y=y, fill=value), data=veg_grid)
  p <- p + scale_fill_gradientn(colours=tim.colors(10))
  p <- p + geom_path(data=na_fort, aes(x=long, y=lat, group=group),  colour='grey55')
  # p <- add_map_albers(p, us.shp.ll, limits)+ coord_fixed()
  p <- p + theme_bw()
  p <- p + theme(axis.text = element_blank(),
                 axis.title = element_blank(),
                 axis.ticks = element_blank())
  p <- p + facet_wrap(~variable, ncol=3)
  p <- p + coord_cartesian(xlim = longlimits, ylim = latlimits)
  p <- p + ggtitle(paste0(age, ' YBP slice'))
  print(p)
  # ggsave('figures/reveals_taxon_50.pdf')
}
dev.off()


veg_sub = veg_grid[which((veg_grid$ages==age)&(veg_grid$variable == 'Larch')),]

p <- ggplot(data=veg_sub)
p <- p + geom_tile(aes(x=x, y=y, fill=value), data=veg_grid)
p <- p + scale_fill_gradientn(colours=tim.colors(10))
p <- p + geom_path(data=na_fort, aes(x=long, y=lat, group=group),  colour='grey55')
# p <- add_map_albers(p, us.shp.ll, limits)+ coord_fixed()
p <- p + theme_bw()
p <- p + theme(axis.text = element_blank(),
               axis.title = element_blank(),
               axis.ticks = element_blank())
p <- p + coord_cartesian(xlim = longlimits, ylim = latlimits)
print(p)
ggsave('figures/larch.pdf')

###################################################################################################################################


library(scatterpie)
library(fields)
library(reshape2)
library(tidyr)
library(ggforce)

cols = unique(veg_pred$taxon)

# geom_scatterpie <- function(mapping=NULL, data, cols, ...) {
#   if (is.null(mapping))
#     mapping <- aes_(x=~x, y=~y)
#   mapping <- modifyList(mapping, aes_(r0=0, fill=~type,
#                                       amount=~value))
#   
#   if (!'r' %in% names(mapping)) {
#     xvar <- as.character(mapping)["x"]
#     size <- diff(range(data[, xvar]))/50
#     mapping <- modifyList(mapping, aes_(r=size))
#   }
#   
#   names(mapping)[match(c("x", "y"), names(mapping))] <- c("x0", "y0")
#   
#   df <- gather_(data, "type", "value", cols)
#   ## df$type <- factor(df$type, levels=cols)
#   geom_arc_bar(mapping, data=df, stat='pie', inherit.aes=FALSE, ...)
# }


veg_cast = dcast(veg_pred, id + x + y + ages ~ taxon, value.var='meansim')

ages_list = ages*100
pdf(file='figures/reveals_UMW_v2.pdf')
for (age in ages_list){
  veg_sub = subset(veg_cast, ages == age)
  veg_sub$x = veg_sub$x*1e6
  veg_sub$y = veg_sub$y*1e6
  
  p <- ggplot()
  p <- p + geom_scatterpie(aes(x=x, y=y),
                           data=veg_sub, cols=cols, alpha=.8)
  p <- add_map_albers(p, us.shp, limits)+ coord_fixed()
  p <- p + theme_bw()
  p <- p + theme(axis.text = element_blank(),
                 axis.title = element_blank(),
                 axis.ticks = element_blank())
  p <- p + ggtitle(paste0(age, " YBP"))
  print(p)
  # ggsave(file=paste0('figures/reveals_UMW_', age, 'ybp.pdf'))
}
dev.off()
####################################################################################################################################
cols = unique(veg_pred$taxon)

veg_cast = dcast(veg_pred, id + x + y + ages ~ taxon, value.var='meansim')
source('r/utils/pie_funs.R')

centers   = centers_veg*1e6
colnames(centers) = c('x', 'y')

xlo = min(centers[,1])
xhi = max(centers[,1])
ylo = min(centers[,2])
yhi = max(centers[,2])
y_shift=40000
x_shift=100000


ages_list = ages*100
pdf(file='figures/reveals_UMW_ppies.pdf', width=12, height=10)
for (age in ages_list){
  veg_sub = subset(veg_cast, ages == age)
  veg_sub$x = veg_sub$x*1e6
  veg_sub$y = veg_sub$y*1e6
  
  
  veg_props  = t(apply(veg_sub[,5:ncol(veg_sub)], 1, function(x) if (sum(x) != 0){x/sum(x)} else {x}))
  coords = veg_sub[,c('x', 'y')]
  
  par(mfrow=c(1,1))
  pieMap(proportions = veg_props, 
         centers  = coords,
         restrict = FALSE,
         inputRestricted = FALSE,
         xlim   = c(xlo+x_shift, xhi-x_shift),
         ylim   = c(ylo+y_shift, yhi-y_shift),
         radius = 18000,
         scale  = 1,
         xlab   = 'x',
         ylab   = 'y',
         add_legend = TRUE, 
         main_title='')
}
dev.off()

ages_list = ages*100
pdf(file='figures/reveals_UMW_v2.pdf')
for (age in ages_list){
  veg_sub = subset(veg_cast, ages == age)
  veg_sub$x = veg_sub$x*1e6
  veg_sub$y = veg_sub$y*1e6
  
  p <- ggplot()
  p <- p + geom_scatterpie(aes(x=x, y=y),
                           data=veg_sub, cols=cols, alpha=.8)
  p <- add_map_albers(p, us.shp, limits)+ coord_fixed()
  p <- p + theme_bw()
  p <- p + theme(axis.text = element_blank(),
                 axis.title = element_blank(),
                 axis.ticks = element_blank())
  p <- p + ggtitle(paste0(age, " YBP"))
  print(p)
  # ggsave(file=paste0('figures/reveals_UMW_', age, 'ybp.pdf'))
}
dev.off()

################################################################################################################################################
## make the grid
################################################################################################################################################

make_grid(meta_pol)

# BROKEN AFTER HERE

# by grid cell
veg_cast = aggregate(meansim ~ idx_veg + ages + taxon, veg_pred, function(x) mean(x))

# veg_cast = dcast(veg_pred, idx_veg + ages ~ taxon, value.var='meansim', fun.aggregate = function(x) sum(x))
# veg_cast[, 3:ncol(veg_cast)] = t(apply(veg_cast[, 3:ncol(veg_cast)], 1, function(x) x/sum(x)))
# 
veg_grid = data.frame(x=centers_veg$x[veg_cast$idx_veg], y=centers_veg$y[veg_cast$idx_veg], veg_cast)
# 
# veg_gridm = melt(veg_grid)


ages_list = ages*100
age=ages_list[1]
veg_sub = subset(veg_grid, ages == age)
veg_sub$x = veg_sub$x*1e6
veg_sub$y = veg_sub$y*1e6

p <- ggplot(data=veg_sub)
p <- p + geom_tile(aes(x=x, y=y, fill=meansim), data=veg_sub) 
p <- p + scale_fill_gradientn(colours=tim.colors())
p <- add_map_albers(p, us.shp, limits)+ coord_fixed()
p <- p + theme_bw()
p <- p + theme(axis.text = element_blank(),
               axis.title = element_blank(),
               axis.ticks = element_blank())
p <- p + facet_wrap(~taxon, ncol=3)
print(p)


veg_sub = subset(veg_grid, (taxon %in% c('HEMLOCK', 'PINE', 'OAK', 'ELM')) & (ages %in% c(200, 500, 1000, 1500, 2000)))
veg_sub = subset(veg_grid, (ages %in% c(200, 500, 1000, 1500, 2000)))

veg_sub$x = veg_sub$x*1e6
veg_sub$y = veg_sub$y*1e6
p <- ggplot(data=veg_sub)
p <- p + geom_tile(aes(x=x, y=y, fill=meansim), data=veg_sub) 
p <- p + geom_point(aes(x=x, y=y, colour=meansim), data=veg_sub, size=1.2) 
p <- p + scale_fill_gradientn(colours=tim.colors())
p <- p + scale_colour_gradientn(colours=tim.colors())
p <- add_map_albers(p, us.shp, limits)+ coord_fixed()
p <- p + theme_bw()
p <- p + theme(axis.text = element_blank(),
               axis.title = element_blank(),
               axis.ticks = element_blank())
p <- p + facet_grid(ages~taxon)
print(p)
ggsave('figures/reveals_UMW_grid.pdf')

# # plot pie maps
# # postscript('r/data/figs/pie_plot_pls_UMW_v0.2.eps', width=8, height=6)
# pdf(paste0('figures/pie_plot_pls_ALL_v', v, '.pdf'), width=12, height=10)
# par(mfrow=c(1,1))
# pieMap(proportions = pls_props, 
#        centers  = knots_in,
#        restrict = FALSE,
#        inputRestricted = FALSE,
#        xlim   = c(xlo+shift, xhi-shift),
#        ylim   = c(ylo+shift, yhi-shift),
#        radius = 14000,
#        scale  = 1,
#        xlab   = 'x',
#        ylab   = 'y',
#        add_legend = FALSE, 
#        main_title='')
# dev.off()


##################################################################################################################################################




dat = read.csv('reveals/REVEALS.Lake.BigJohnPond.csv', row.names=NULL, header=FALSE)
dat = dat[-c(2,3),]
dat_flip = t(dat[, -1])
colnames(dat_flip) = as.vector(dat[,1])
# write.table(dat_flip, 'revealsR/reveals_BigJohnPond.csv', col.names=FALSE, sep=',')
write.csv(dat_flip, 'revealsR/reveals_BigJohnPond.csv', row.names=FALSE)


# with english csv files
a <- REVEALSinR(pollenFile = "revealsR/reveals_BigJohnPond.csv", 
                pf="revealsR/reveals_inputs.csv", 
                dwm="LSM unstable", 
                tBasin="lake", 
                dBasin=600, 
                regionCutoff=100000, 
                repeats=1000)

foo = melt(a, id.vars=c('Pollen.file', 'Parameter.file', 'Distance.weighting', 'Basin.type', 'ages'))
foo$variable = gsub('OTHER.', 'OTHER_', foo$variable)

foo$taxon = unlist(lapply(foo$variable, function(x) unlist(strsplit(as.character(x), ".", fixed=TRUE))[1]))
foo$type = unlist(lapply(foo$variable, function(x) unlist(strsplit(as.character(x), ".", fixed=TRUE))[2]))

ggplot(data=subset(foo, type %in% 'meansim')) + geom_line(aes(x=ages, y=value, colour=taxon)) + theme_bw()


# with english csv files
gpm <- REVEALSinR(pollenFile = "revealsR/reveals_BigJohnPond.csv", 
                pf="revealsR/reveals_inputs.csv", 
                dwm="GPM neutral", 
                tBasin="lake", 
                dBasin=600, 
                regionCutoff=100000, 
                repeats=1000)

foo = melt(gpm, id.vars=c('Pollen.file', 'Parameter.file', 'Distance.weighting', 'Basin.type', 'ages'))
foo$variable = gsub('OTHER.', 'OTHER_', foo$variable)

foo$taxon = unlist(lapply(foo$variable, function(x) unlist(strsplit(as.character(x), ".", fixed=TRUE))[1]))
foo$type = unlist(lapply(foo$variable, function(x) unlist(strsplit(as.character(x), ".", fixed=TRUE))[2]))

ggplot(data=subset(foo, type %in% 'meansim')) + geom_line(aes(x=ages, y=value, colour=taxon)) + theme_bw()


# after REVEALS on all sites do pie plot as in STEPPS calibration paper







# # with english csv files
# a <- REVEALSinR(pollenFile = "REVEALSinR Lake Viitina_en.csv", 
#                 pf="REVEALSinR PPEs_en.csv", 
#                 dwm="LSM unstable", 
#                 tBasin="lake", 
#                 dBasin=600, 
#                 regionCutoff=100000, 
#                 repeats=1000)
