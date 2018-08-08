library(DISQOVER)
setwd('~/Reveals_NEUS/')

distances <- c(seq(10000,100000,10000),seq(2e+05,9e+05,1e+05))
sensitivity_region1 <- replicate(10,
sapply(distances,function(x) {
    a <- REVEALSinR(pollenFile = "data/reveals_input.csv",
                pf         = "data/reveals_input_params.csv",
                filetype   = "csv",
                dwm        = "gpm neutral",
                tBasin     = "lake",
                dBasin     = 1000, # diameter!
                regionCutoff = x,
                repeats      = 1000)
    a.mean <- a[,grep('meansim',colnames(a))]
}))  


colnames(sensitivity_region) <- distances



ten.km <- matrix(ncol=10,unlist(as.numeric(sensitivity_region1[,1,])))
apply(ten.km,1,var)

rep1 <- matrix(ncol=19,unlist(as.numeric(sensitivity_region1[,,1])))
apply(rep1,1,var)

sensitivity_region1[,c(1,19),]
sensitivity_region1[,c(1,10),]
sensitivity_region1[,c(5,10,13),]
