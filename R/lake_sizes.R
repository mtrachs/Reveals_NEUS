#---------------------------------------------------------------------------------------------------------------------
#find lake sizes 
#---------------------------------------------------------------------------------------------------------------------

lake_size <- read.csv('~/Reveals_NEUS/reveals-na/data/area_lakes_1.1.csv')
lake_names <- read.csv('~/workflow_stepps_calibration/expert_elicitation/data/site_names.csv')


lake_available_index <- lake_names$site.names%in%lake_size$site.name
lakes_available <- lake_names[lake_available_index,]
lakes_not_available <- lake_names[lake_available_index==FALSE,]

write.csv(lakes_not_available ,'~/Reveals_NEUS/data/lake_without_lake_size_1_10_2018.csv',row.names=FALSE)  
