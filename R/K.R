library(zipfR)
library(rioja)
setwd('~/Reveals_NEUS/')
plot.loc <- 'figures/'
#--------------------------------------------------------------------------------------------------------------]
#function to estimate species dependent depositional coefficient
Ki <- function(b, R, zmax=400000){
  ul1 <- b*(zmax-R)^(1/8)
  ul2 <- b*(zmax+R)^(1/8)
  ul3 <- b*(2*R)^(1/8)
  gamma_Ki <- Igamma(8, ul1, lower=TRUE) -  
    Igamma(8, ul2, lower=TRUE) +  
    Igamma(8, ul3, lower=TRUE)
  
  return(4*pi*R/b^8*gamma_Ki)
}

K_calc <- function(fall_speed,R,zmax,n,u,c){
  
  #calculate parameter b later on used to estimate deposition coeffcient 
  b <- 4/sqrt(pi) * fall.speed/(n*u*c)  
  
  #species specific depositional coefficient
  K_species <- 
    sapply(b,function(x){
      Ki(x[[1]],R=R,zmax = zmax)
    })
  
  #eq (5) in Sugita (2007a)
}

#----------------------------------------------------------------------------------------------------------------
#parameters

reveals.params <- read.csv('data/reveals_input_params_variable.csv')
taxa <- reveals.params$species

#fall speed
fall.speed <- reveals.params$fallspeed
names(fall.speed) <- taxa
fall.speed <- as.data.frame(t(fall.speed))

n=0.25
u=3
c=0.12

b <- 4/sqrt(pi) * fall.speed/(n*u*c)
max_dist <- c(50000,seq(100000,400000,100000))
taxa1 <- taxa[!taxa%in%c('Elm','Poplar')]

#########################################################################################################################
K <-
  sapply(max_dist,function(x){
    K_calc(fall_speed = fall.speed, R = 200,zmax = x,n=n,u=u,c=c)
  })

pdf(paste(plot.loc,'K.pdf',sep=''),height=7,width =10)
matplot(max_dist/1000,t(K),type='p',pch =1:14,xlim=c(50,500),ylab='',xlab='',col=1:14)
legend('topright',legend =taxa,pch = 1:14,col=1:14,cex=0.7)
mtext(side = 1,'m',line=2.2)
mtext(side = 2,'K/1000',line=2.2)
dev.off()



# plot(1:15,seq(1,10,length.out=15),type='n',ylab='',xlab='')
# 
# 
# k_spec_100 <- 
#   sapply(b,function(x){
#     Ki(x[[1]],R=0.5,zmax = 100)
#   })
# 
# plot_order <- names(k_spec_100[order(k_spec_100)])
# 
# plot(1:15,seq(1,10,length.out=15),type='n',ylab='',xlab='')
# 
# 
# for(i in 1:6) {
# 
#   k_spec <- 
#   sapply(b,function(x){
#     Ki(x[[1]],R=0.5,zmax = i*100)
#   })
# 
#   rel_k <- k_spec/k_spec['Tamarack']
# 
#   points(rel_k[plot_order],pch = 16,col=i)
# }
# 

