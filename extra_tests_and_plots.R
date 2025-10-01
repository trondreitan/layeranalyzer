library(layeranalyzer)
load("extra1b.Rdata")

# Check if second layer is still necessary      
do.temp.size.cause1.1layer=layer.analyzer(do.struct, temp.struct,
                         X.struct.1layer,
		         causal=cbind(c(1,1,3,1),c(2,1,3,1)),
                         do.maximum.likelihood = TRUE, mcmc=T,
                         maximum.likelihood.numstart = 1000, 
                         num.MCMC=1000,spacing=10,burnin=4000,
                         num.temp=1)
anova(do.temp.size.cause1.1layer, do.temp.size.cause1)


# Check if stochasticity is still found in the first layer:
do.temp.size.cause1.1det=layer.analyzer(do.struct, temp.struct,
                         X.struct.1det,
		         causal=cbind(c(1,1,3,1),c(2,1,3,1)),
                         do.maximum.likelihood = TRUE, mcmc=T,
                         maximum.likelihood.numstart = 1000, 
                         num.MCMC=1000,spacing=10,burnin=4000,
                         num.temp=1)
anova(do.temp.size.cause1.1det, do.temp.size.cause1)
