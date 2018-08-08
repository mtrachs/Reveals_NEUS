#---------------------------------------------------------------------------------------------------------------------
#find lake sizes 
#---------------------------------------------------------------------------------------------------------------------

lake_size <- read.csv('~/Reveals_NEUS/data/area_lakes_1.2.csv')
lake_size <- read.csv('~/Reveals_NEUS/reveals-na/data/area_lakes_1.1.csv')
lake_names <- read.csv('~/workflow_stepps_calibration/expert_elicitation/data/site_names.csv')


lake_available_index <- lake_names$site.names%in%lake_size$site.name
lakes_available <- lake_names[lake_available_index,]
lakes_not_available <- lake_names[lake_available_index==FALSE,]

write.csv(lakes_not_available ,'~/Reveals_NEUS/data/lake_without_lake_size_1_10_2018.csv',row.names=FALSE)  
write.csv(lakes_available ,'~/Reveals_NEUS/data/lake_available.csv',row.names=FALSE)  


lake_with_size <- lake_size[lake_size$stid%in%lakes_available$site.id,]
lake_with_size <- lake_with_size[order(lake_with_size$stid),]
write.csv(lakes_available ,'~/Reveals_NEUS/data/lakes_with_size.csv',row.names=FALSE)  
