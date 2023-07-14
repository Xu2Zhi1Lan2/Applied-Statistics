####################################
## Example 1: Additive Model
####################################
library(mlbench)
data(Ozone) # use temperature, inversion base height
            # and inversion base temperature to predict O3 concentration

O3 = Ozone$V4
temp = Ozone$V9
ibh = Ozone$V10
ibt = Ozone$V12

## fit a linear regression
lr <- lm(O3 ~ temp + ibh + ibt)
summary(lr)

plot(resid(lr)~temp[!is.na(ibh)&!is.na(ibt)&!is.na(temp)&!is.na(O3)])
plot(resid(lr)~ibt[!is.na(ibh)&!is.na(ibt)&!is.na(temp)&!is.na(O3)])
plot(resid(lr)~ibh[!is.na(ibh)&!is.na(ibt)&!is.na(temp)&!is.na(O3)])

## fit an additive model (use gam package)
library(gam) 

gam1 <- gam(O3 ~ lo(temp) + lo(ibh) + lo(ibt)) # lo means loess smoother
summary(gam1)
R_sq1 = 1 - 4246.151/14073.83 # better fit compared with Linear regression

plot(gam1,residuals=TRUE,se=TRUE,pch=".") # transformations on predictors

## fit an additive model (use mgcv package)
library(mgcv) 

gam2 <- gam(O3 ~ s(temp) + s(ibh) + s(ibt)) # s() means smooth term
summary(gam2)
plot(gam2) # transformations on predictors

gam3 <- gam(O3 ~ s(temp,ibh) + s(ibt)) # include interaction 
summary(gam3)
plot(gam3)
vis.gam(gam3,theta=-45,color="gray")

## residual
plot(predict(gam2),residuals(gam2), xlab="Predicted", ylab="Residuals")
qqnorm (residuals (gam2), main="") # looks normal

####################################
## Example 2: Generalized Additive Model
####################################
gam.pois <- gam(O3 ~ s(temp) + s(ibh) + s(ibt), family=poisson, scale=-1) 
## scale = -1 means dispersion should be estimated rather than fixed at one
plot(gam.pois, residuals=TRUE) 



