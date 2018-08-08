setwd('~/Reveals_NEUS/')
path_plot <- 'figures/'
path_data <- 'data/'
wd <- getwd()

if (!file.exists(path_plot)){
  dir.create(file.path(wd, path_plot))
}




#-----------------------------------------------------------------------------------------------------------------
#dispersal kernel
#-----------------------------------------------------------------------------------------------------------------
reveals.params <- read.csv(paste(path_data,'reveals_input_params_variable.csv',sep=''))
taxa <- reveals.params$species
fall.speed <- reveals.params$fallspeed
names(fall.speed) <- reveals.params$species
n <- 0.25
u <-  3
c <- 0.12
gamma <- n/2
b <- 4/sqrt(pi) * fall.speed/(n*u*c)
names(b) <- taxa
b <- as.data.frame(b)

radius <- seq(0.5,400,1)

#----------------------------------------------------------------------------------------------------------------
#plot dispersal function
#--------------------------------------------------------------------------------------------------------------------

g <- 
    (sapply(radius,function(z){
      b * gamma*z^(gamma-1) *exp(-b*z^gamma)
    }))

g <- as.data.frame(g)
g <- t(g)

sapply(taxa,function(x) plot(radius,g[,x],type='b',pch = 16,cex = 0.5,main=x,ylim=c(0,0.05)))

sapply(taxa,function(x) plot(radius,cumsum(g[,x])/sum(g[,x]),type='b',pch = 16,cex = 0.5,main=x))


#-----------------------------------------------------------------------------------------------------------------
#one dimensional integral of dispersal function
#-----------------------------------------------------------------------------------------------------------------

basin_radius <- c(0.01,0.1,0.5,1,5)
pollen_dispersal_distance <- c(50,400,700,1000)
lake_radius <- 0.5

important.taxa <- c('Beech','Hemlock','Oak','Spruce')

#pdf(paste(path_plot,'Figure_5.pdf',sep=''),height=10,width=15)
pdf(paste(path_plot,'Figure_5.pdf',sep=''),height=8,width=8)
#par(mfrow=c(3,5))
par(mfrow=c(2,2))
catchment.radius <- 
  sapply(fall.speed[names(fall.speed)%in%important.taxa],function(v){
    lapply(pollen_dispersal_distance,function(pollen_disp){
      par(mar=c(3,3,2,1))

      radius <- seq(0.5,pollen_disp,1)  
      b <- 4/sqrt(pi) * v/(n*u*c)
      integrand <- function(z)   {b * gamma*z^(gamma-1) * exp(-b*z^gamma)}# * z^2}
      integration <- sapply(radius[radius>=lake_radius],function(rad) integrate(f=integrand,
                                                                              lower=lake_radius,upper = rad)$value)
      if(pollen_disp==min(pollen_dispersal_distance)){
      plot(radius[radius>=lake_radius],integration/max(integration),pch = 16,ylim=c(0,1),
           main = paste(taxa[fall.speed==v]),type='l',lwd=2,xlim=c(0,max(pollen_dispersal_distance)))
      mtext(side=1,text='Distance from lake centre [km]',line=2.2)
      mtext(side=2, text= 'Proporiton of pollen deposited',line=2.2)
       if(v==fall.speed[[which(names(fall.speed)=='Beech')]]) {
         legend('bottomright',col=1,lty = 1:4,lwd=2,legend = paste('Dispersal distance', pollen_dispersal_distance ,'km',sep=' '))
       }
      }
      else {
        points(radius[radius>=lake_radius],integration/max(integration),pch = 16,ylim=c(0,1),
             main = paste(taxa[fall.speed==v]),type='l',col=1, lty = which(pollen_disp==pollen_dispersal_distance),lwd=2)#, 'basin_radius =',lake_radius, 'km')) 
      }
      
      integration.prop <- integration/max(integration)

      ninety <- radius[which.min(abs(integration.prop-0.9))]
      fifty <- radius[which.min(abs(integration.prop-0.5))]
      catchment.radius <- data.frame(fifty = fifty, ninety=ninety)
  })
})
dev.off()

colnames(catchment.radius) <- taxa
rownames(catchment.radius) <- c(basin_radius)


 
 
#-------------------------------------------------------------------------------------------------------------------
cdf <- cumsum(g$Hemlock)/sum(g$Hemlock)

plot(radius,cdf,type='l')

#-------------------------------------------------------------------------------------------------------------------
# #area weighting

 
sapply(taxa,function(x) { 
  
  radius <- seq(0.1,700,0.1)
  area <- pi*radius^2
  
  g <- 
    t(sapply(radius,function(z){
      b * gamma*z^(gamma-1) *exp(-b*z^gamma)
    }))
  
  g <- as.data.frame(g)
  g_weight <-g[-nrow(g),x]*diff(area)
  cdf_weight <- cumsum(g_weight)/sum(g_weight)
  plot(radius[-length(radius)],cdf_weight,type='l',main=x)
})
