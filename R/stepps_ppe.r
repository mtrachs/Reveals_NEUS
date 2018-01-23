#----------------------------------------------------------------------------------------------------------------------
#prepare STEPPS PPEs for use in REVEALS
#----------------------------------------------------------------------------------------------------------------------
library(rstan)
library(abind)

setwd('~/workflow_stepps_calibration/results/')
path_data <- 'data/stepps_median/'

file.names <- list.files(data.location)
files.use <- file.names[grep('Ka_Kgamma_EPs',file.names)]
file.name <- strsplit(files.use[1],'0.')[[1]][1]

taxa <- sort(c('Chestnut','Ash','Elm','Beech','Birch','Oak','Pine','Spruce','Hemlock','Hickory',
               'Other hardwood','Other conifer','Tamarack','Poplar','Maple'))  
# load PPEs of Ka_Kgamma
for (i in 0:11) {
  fname = paste(path_data,file.name,i,'.csv',sep='') #sprintf('%s/%s.csv', path_out, suff_fit)
  fit <- read_stan_csv(fname)
  #fit1 <- read_stan_csv(fname1)
  if (i==0)  {post = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)}
  else{post_new <-  rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
  post <- abind(post,post_new,along=1)
  }
}

par_names <- colnames(post[,1,])
phi.index <- grep('phi',par_names)

ppe_STEPPS_mean <- colMeans(post[,1,phi.index])
ppe_STEPPS_sd <- apply(post[,1,phi.index],2,sd)

names(ppe_STEPPS_mean) <- names(ppe_STEPPS_sd) <- taxa

ppe_STEPPS_mean_stand <- ppe_STEPPS_mean/ppe_STEPPS_mean['Oak']
ppe_STEPPS_sd_stand <- ppe_STEPPS_sd/ppe_STEPPS_sd['Oak']

ppe<- as.data.frame(cbind(ppe_STEPPS_mean_stand,ppe_STEPPS_sd_stand))

saveRDS(ppe,'~/Reveals_NEUS/data/ppe_stepps.RDS')

#----------------------------------------------------------------------------------------------------------------
#compare_ppe literature and STEPPS
#----------------------------------------------------------------------------------------------------------------
taxa1 <- taxa
taxa1[taxa1%in%c('Other conifer','Other hardwood','Tamarack')] <- c('Fir','Alder','Larch')
taxa2 <- taxa[taxa!='Chestnut']

ppe_lit <- readRDS('data/PPEs_agg.RDS')
ppe_lit <- ppe_lit[which(ppe_lit$taxon %in% taxa1),] #no PPE for Chestnut
ppe_lit$taxon <- as.character(ppe_lit$taxon)
ppe_lit$taxon[ppe_lit$taxon%in%c('Fir','Alder','Larch')] <- c('Other hardwood','Other conifer','Tamarack') 
ppe_lit <- ppe_lit[order(ppe_lit$taxon),]


par_stats_stepps = data.frame(name=taxa2, mu=ppe$ppe_STEPPS_mean_stand[taxa!='Chestnut'], handle=rep('STEPPS', length(taxa2)))
par_stats_reveals = data.frame(name=taxa2, mu=ppe_lit$ppe, handle=rep('REVEALS', length(taxa2)))

dat = rbind(par_stats_stepps,par_stats_reveals)

p <- ggplot(data=dat, aes(x=reorder(name, mu), y=mu, group=handle, colour=handle)) + 
  geom_point(size=4, position=position_dodge(width=0.5)) + 
  xlab("Taxon") + ylab(parse(text='PPE')) +
  coord_flip() + theme_bw() + theme(axis.title.x=element_text(size=20), 
                                    axis.title.y=element_text(size=20), 
                                    axis.text.x=element_text(size=rel(1.3)),
                                    axis.text.y=element_text(size=rel(1.3)))

print(p)
ggsave('figures/comparison_ppe.pdf')
#-----------------------------------------------------------------------------------------------------------------
# comparison gamma and fall speeds
setwd('~/workflow_stepps_calibration/results/')
path_data <- 'data/stepps_median/'

file.names <- list.files(path_data)
files.use <- file.names[grep('Ka_Kgamma_EPs',file.names)]
file.name <- strsplit(files.use[1],'0.')[[1]][1]

# load PPEs of Ka_Kgamma
for (i in 0:11) {
  fname = paste(path_data,file.name,i,'.csv',sep='') #sprintf('%s/%s.csv', path_out, suff_fit)
  fit <- read_stan_csv(fname)
  #fit1 <- read_stan_csv(fname1)
  if (i==0)  {post = rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)}
  else{post_new <-  rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
  post <- abind(post,post_new,along=1)
  }
}

par_names <- colnames(post[,1,])
gamma.index <- grep('gamma',par_names)[1:length(taxa)]

gamma_STEPPS_mean <- colMeans(post[,1,gamma.index])
gamma_STEPPS_sd <- apply(post[,1,gamma.index],2,sd)

names(gamma_STEPPS_mean) <- names(gamma_STEPPS_sd) <- taxa
gamma_STEPPS_mean1 <- gamma_STEPPS_mean[names(gamma_STEPPS_mean)!='Chestnut']
#------------------------------------------------------------------------------------------------------------------
#load fallspeeds
setwd('~/Reveals_NEUS/')
svs = read.csv('reveals-na/data/svs_LC6K.csv', sep=',', header=TRUE, stringsAsFactors=FALSE)
#svs$taxon[svs$taxon=='Fir'] <- 'Abies'
svs[which(is.na(svs$sv)), 'sv'] = 0.01
svs_agg = aggregate(sv ~ taxon, svs, median)


svs_agg = svs_agg[svs_agg$taxon %in% taxa1, ]
svs_agg$taxon[svs_agg$taxon%in%c('Fir','Alder','Larch')] <- c('Other hardwood','Other conifer','Tamarack')
svs_agg <- svs_agg[order(svs_agg$taxon),]
#------------------------------------------------------------------------------------------------------------------
pdf('~/Reveals_NEUS/figures/gamma_fallspeed.pdf',height=5,width=5)
par(oma=c(2.2,2.2,0,0),mar=c(1,1,1,1))
plot(svs_agg$sv,gamma_STEPPS_mean[names(gamma_STEPPS_mean)%in%svs_agg$taxon],pch = 16,xlim=c(0,0.15))
mtext(side=2,text=expression(gamma*' [%]'),line=2.2,font=2)
mtext(side=1,text='Fallspeed [m/s]',line=2.2,font=2)
text(svs_agg$sv,gamma_STEPPS_mean[names(gamma_STEPPS_mean)%in%svs_agg$taxon],
     labels=as.character(names(gamma_STEPPS_mean)[names(gamma_STEPPS_mean)%in%svs_agg$taxon]),pos=4)
dev.off()
     