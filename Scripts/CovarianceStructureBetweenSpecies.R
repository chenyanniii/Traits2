library(ape)
library(nlme)
library(geiger)

### The covariance (correlation) structure between species is permitted to match 
### that expected under a Brownian motion process of evolution on the tree.

## define the covariance structure under Brownian motion process 
bm <- corBrownian(1, TreeAllMatrix)
bm

## a simple analysis investigating the relationship between seed mass and seed height

### Under Brownian Motion 
modelofmass1 <- gls(scaled_Mass ~ scaled_Height, data=AllMatrix_G_scaled, correlation=bm)
summary(modelofmass1)

modelofgermination11 <- gls(scaled_Germination ~ scaled_Mass, data=AllMatrix_G_scaled, correlation=bm)
summary(modelofgermination11)

modelofgermination12 <- gls(scaled_Germination ~ scaled_Height, data=AllMatrix_G_scaled, correlation=bm)
summary(modelofgermination12)

## relaxing the assumption of Brownian motion
modelofmass2 <- gls(scaled_Mass ~ scaled_Height, data=AllMatrix_G_scaled, 
                    correlation=corPagel(1,TreeAllMatrix))
summary(modelofmass2)

modelofgermination21 <- gls(scaled_Germination ~ scaled_Mass, data=AllMatrix_G_scaled, correlation=corPagel(1,TreeAllMatrix))
summary(modelofgermination21)

modelofgermination22 <- gls(scaled_Germination ~ scaled_Height, data=AllMatrix_G_scaled, correlation=corPagel(1,TreeAllMatrix))
summary(modelofgermination22)
