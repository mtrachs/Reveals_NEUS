library(reshape2)
library(ggplot2)
setwd('~/Reveals_NEUS/')

PPEs = read.csv('reveals-na/data/PPEs.csv')
PPEs = PPEs[which(!is.na(PPEs$ppe)),]
PPEs = PPEs[which(PPEs$use == 'y'),] # expert elicitation of ppes to use

N_datasets = length(unique(PPEs$dataset))
pubs = as.vector(unique(PPEs$pub))

# no taxon that is common to all data sets
datasets = unique(PPEs$dataset)

#table(PPEs$taxon)
PPEs_cast = dcast(PPEs, dataset ~ taxon, value.var="ppe")


ppe_error_mean = mean(PPEs$error, na.rm=TRUE)

# only datasets 3 and 4 do not have birch
#if there is no error estimate, we assign the error mean, is that a good idea?
PPEs[which(is.na(PPEs$error)),'error'] = ppe_error_mean
PPEs[which(PPEs$error==0),'error'] = ppe_error_mean #sometimes errors are huge compared to estimate, is unrealistic

#set reference taxon
ref = 'Oak'
unique(PPEs[PPEs$taxon == ref,'dataset'])

# rescale using Oak as the reference taxon
#for (i in 1:N_datasets){ # is this correct? I would replace it by unique(PPEs$dataset) 
#with this code we miss three datasets 19,28,29
for (i in unique(PPEs$dataset)){
    if (any((PPEs$dataset == i) & (PPEs$taxon == ref))){
    PPEs[PPEs$dataset == i, 'ppe_scaled'] = PPEs[PPEs$dataset == i, 'ppe'] / 
      PPEs[which((PPEs$dataset == i) & (PPEs$taxon == ref)), 'ppe']
    PPEs[PPEs$dataset == i, 'error_scaled'] = PPEs[PPEs$dataset == i, 'error'] / 
      PPEs[((PPEs$dataset == i) & (PPEs$taxon == ref)), 'error']
  } else {
    PPEs[which(PPEs$dataset == i), 'ppe_scaled'] = NA
    PPEs[which(PPEs$dataset == i), 'error_scaled'] = NA
  }
}


# PPEs[,c('dataset', 'taxon', 'ppe', 'error', 'ppe_scaled', 'error_scaled')]

# mean second ref
ref2 = 'Birch'
# mean_ref2 = mean(PPEs[which(PPEs$taxon == 'Willow'), 'ppe_scaled'], na.rm=TRUE)
mean_ref2 = mean(PPEs[which(PPEs$taxon == ref2), 'ppe_scaled'], na.rm=TRUE)

# rescale using Birch as the reference taxon (in fact this is only filling in NAs...)if there is no oak,and everything is rescaled to oak)
for (i in 1:N_datasets){
  dataset = datasets[i] # is missing above!!!!!
  if (any((PPEs$dataset == dataset) & (is.na(PPEs$ppe_scaled)))){
    ppe_rescale = PPEs[((PPEs$dataset == dataset)&(PPEs$taxon == ref2)), 'ppe']/mean_ref2 #remove which
    error_rescale = PPEs[((PPEs$dataset == dataset)&(PPEs$taxon == ref2)), 'error']/mean_ref2 # I think this is error not ppe 
  
    PPEs[PPEs$dataset == dataset, 'ppe_scaled'] = PPEs[PPEs$dataset == dataset, 'ppe']/ppe_rescale
    PPEs[PPEs$dataset == dataset, 'error_scaled'] = PPEs[PPEs$dataset == dataset, 'error']/error_rescale
  }
}


PPE_cast = dcast(PPEs, taxon ~ dataset, value.var = 'ppe_scaled', fun.aggregate = mean, fill = 0)
PPE_cast[,2:ncol(PPE_cast)] = round(PPE_cast[,2:ncol(PPE_cast)], 2)
PPE_cast[(PPE_cast == 0.00)] = NA
write.csv(PPE_cast, 'data/PPE_compare_table.csv')


ggplot(data=PPEs[PPEs$taxon!='Cedar',]) + geom_point(aes(x=ppe_scaled, y=reorder(taxon,ppe_scaled), colour=factor(dataset)))
ggsave('figures/ppes_all.pdf')
# PPEs_sub = PPEs[,c('taxon', 'ppe_scaled', 'error_scaled')]
# PPEs_sub[which(PPEs_sub$taxon == 'Larch'), 'error_scaled'] = 0.01
# PPEs_sub = subset(PPEs_sub, is.finite(error_scaled))
# PPEs_agg = aggregate(ppe_scaled ~ taxon, PPEs_sub, median)
# PPEs_agg$error = aggregate(error_scaled ~ taxon, PPEs_sub, median, na.rm=TRUE)$error_scaled
# 
# saveRDS(PPEs_agg,  'data/PPEs_agg.RDS')



# Apply Furong's PPE selection protocol
PPEs_sub = PPEs[,c('taxon', 'ppe_scaled', 'error_scaled')]
PPEs_sub[which(PPEs_sub$taxon == 'Larch'), 'error_scaled'] = 0.01
# PPEs_sub = subset(PPEs_sub, is.finite(error_scaled))

taxon_list = sort(unique(PPEs_sub$taxon))
PPE_st3 = data.frame(taxon=taxon_list, ppe=rep(NA), error=rep(NA))


#if there are more than 2 values for taxon exclude the value with largest absolute distance to the mean
for (k in 1:length(taxon_list)){
  idx_taxon = which(PPEs_sub$taxon == taxon_list[k])
  
  if(length(idx_taxon)==1){
    PPE_st3[which(PPE_st3$taxon == taxon_list[k]), 'ppe']   = PPEs_sub[idx_taxon, 'ppe_scaled'] 
    PPE_st3[which(PPE_st3$taxon == taxon_list[k]), 'error'] = PPEs_sub[idx_taxon, 'error_scaled'] 
  } else if ((length(idx_taxon) %in% c(2))){
    PPE_st3[which(PPE_st3$taxon == taxon_list[k]), 'ppe']   = mean(PPEs_sub[idx_taxon, 'ppe_scaled'], na.rm=TRUE)
    PPE_st3[which(PPE_st3$taxon == taxon_list[k]), 'error'] = mean(PPEs_sub[idx_taxon, 'error_scaled'], na.rm=TRUE) 
  } else if ((length(idx_taxon) > 2)){
    idx_remove = which.max(abs(PPEs_sub[idx_taxon,'ppe_scaled'] - mean(PPEs_sub[idx_taxon,'ppe_scaled'], na.rm=TRUE)))
    idx_taxon = idx_taxon[-idx_remove]
    PPE_st3[which(PPE_st3$taxon == taxon_list[k]), 'ppe']   = mean(PPEs_sub[idx_taxon, 'ppe_scaled'], na.rm=TRUE)
    PPE_st3[which(PPE_st3$taxon == taxon_list[k]), 'error'] = mean(PPEs_sub[idx_taxon, 'error_scaled'], na.rm=TRUE) 
  }
}

saveRDS(PPE_st3,  'data/PPEs_agg.RDS')

ggplot(data=PPE_st3[PPE_st3$taxon!='Cedar',]) + geom_point(aes(x=ppe, y=reorder(taxon,ppe)))
ggsave('figures/ppes_averaged.pdf')

# 
# PPEs_agg = aggregate(ppe_scaled ~ taxon, PPEs_sub, median)
# PPEs_agg$error = aggregate(error_scaled ~ taxon, PPEs_sub, median, na.rm=TRUE)$error_scaled
# 
# saveRDS(PPEs_agg,  'data/PPEs_agg.RDS')

##################################################################################################################################


# # PPEs = PPEs[which(!(PPEs$taxon %in% c('Chenopods', 'Cedar'))),]
# ggplot(data=PPEs) + geom_point(aes(x=ppe_oak, y=taxon, colour=factor(dataset)))
# ggsave('figures/ppes_remove_large.pdf')





