###########################################################################################################################
# in this file, we compare pedictions of REVEALS for different distances of regional vegetation
###########################################################################################################################
library(rioja)

#load the data not aggregated at regional level so that we have 79 sites to compare
reveals_pred <-
lapply(c(50000,1e+05,4e+05),function(x){
  reveals_reconstruction <- readRDS(paste('~/Reveals_NEUS/output/veg_pred_ppe_literature_fallspeed_measured_composition_no_larch_lake_size_correct_dwm_gpm neutral_',x,'_regional_vegetation.RDS',sep=''))           
  reveals_coord <- readRDS('~/Reveals_NEUS/output/veg_pred_ppe_literature_fallspeed_measured_composition_no_larch_lake_size_correct_dwm_gpm neutral_50000_regional_vegetation_aggregated_coord.RDS')

  reveals_median <- reveals_reconstruction[,grep('median',colnames(reveals_reconstruction))] 
  reveals_mean <- reveals_reconstruction[,grep('mean',colnames(reveals_reconstruction))] 
  reveals_mean <- reveals_mean/100
})



#----------------------------------------------------------------------------------------------------------------------
#chord distances

reveals_distance <- matrix(ncol=3,nrow=nrow(reveals_pred[[1]]))

for(x in 1:nrow(reveals_pred[[1]])){
    reveals_distance[x,1] <- as.matrix(paldist((rbind(reveals_pred[[1]][x,],reveals_pred[[2]][x,]))))[1,2]
    reveals_distance[x,2] <- as.matrix(paldist((rbind(reveals_pred[[1]][x,],reveals_pred[[3]][x,]))))[1,2]
    reveals_distance[x,3] <- as.matrix(paldist((rbind(reveals_pred[[3]][x,],reveals_pred[[2]][x,]))))[1,2]
}

reveals_distance <- round(reveals_distance,3)

colMeans(reveals_distance)


