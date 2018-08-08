#----------------------------------------------------------------------------------------------------------------
#estimating volume of a three dimensional sphere (radius = 1)
# 
# basic idea: I divide the sphere into i small segments with radius ri (hence area ri^2 *pi) and height sin(ri).
# The volume of these small segments is: ri^2*pi*sin(ri). The total volume is the sum of these small segments 

x <- seq(0,1,0.01)
deg <- acos(x)


z <- sin(deg)
ri <- cos(deg)# same as x

area <- ri^2*pi
diff.area <- diff(area)

vol <- sum(diff.area*z[-length(z)]) #2.11 (is more precise with increasing number of segments)
precise.vol <- 4/6*pi #2.09

#in a first approximation these values are similar

#-----------------------------------------------------------------------------------------------------------------
# Now I use the same idea for pollen volume / contribution of pollen within a certain radius
# instead of sin(ri) I use g(ri) resulting in ri^2*pi*g(ri)

reveals.params <- read.csv('~/Reveals_NEUS/data/reveals_input_params_variable.csv')
taxa <- reveals.params$species
fall.speed <- reveals.params$fallspeed
names(fall.speed) <- taxa
fall.speed <- as.data.frame(t(fall.speed))
n <- 0.25
u <-  3
c <- 0.12
gamma <- n/2



sapply(taxa,function(x){

    b <- 4/sqrt(pi) * fall.speed[[x]]/(n*u*c)
    radius <- seq(0.01,700,0.01)

    g <- 
      sapply(radius,function(z){
        b * gamma*z^(gamma-1) *exp(-b*z^gamma)
      })

    g <- unlist(g)
    area <- radius^2*pi
    diff.area <- diff(area)
    g_weight <-g[-length(g)]*diff.area
    cdf_weight <- cumsum(g_weight)/sum(g_weight)
    plot(radius[-length(radius)],cdf_weight,type='l',ylim=c(0,1),main=x)
    mtext(side=1,text='radius',line=2.2)
    mtext(side=2, text= 'proportion of deposited pollen',line=2.2)
})


#---------------------------------------------------------------------------------------------------------------------
#plot deposition function
sapply(taxa,function(x){
  
  b <- 4/sqrt(pi) * fall.speed[[x]]/(n*u*c)
  radius <- seq(0.1,700,0.1)
  
  g <- 
    sapply(radius,function(z){
      b * gamma*z^(gamma-1)*exp(-b*z^gamma)
    })
  
  g <- unlist(g)
  
  plot(radius,g,type='b',ylim=c(0,0.5),main=x,pch = 16,cex=0.5)
  mtext(side=1,text='radius',line=2.2)
  mtext(side=2, text= 'Deposition function',line=2.2)    
})
