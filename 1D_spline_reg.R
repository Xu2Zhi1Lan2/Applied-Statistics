### This file is used to learn the implementation of univariate smoothing spline regression
### loading required packages
library(mgcv)
library(HRW)

### generate spline basis
data(WarsawApts)
x <- WarsawApts$construction.date
a <- 1.01*min(x)- 0.01*max(x)
b <- 1.01*max(x)- 0.01*min(x) # [a,b] is slightly wider than min(x) to max(x)
numIntKnots <- 20 # remember number of interior knots = number of basis - 2
intKnots <- quantile(unique(x),seq(0,1,length = 
                                     (numIntKnots+2))[-c(1,(numIntKnots+2))]) # the default setting

xg <- seq(a,b,length = 1001) # grid on [a,b]
Zg <- ZOSull(xg,range.x = c(a,b),intKnots = intKnots) # evaluate the 22 basis functions with given knots on the grid xg
plot(0,type = "n",xlim = range(xg),ylim = range(Zg),bty = "l",
     xlab = "construction date (year)",ylab = "spline basis function") # create a new plot with corresponding axis
for (k in 1:ncol(Zg)) lines(xg,Zg[,k],col = k,lwd = 2) # add the splines
for (k in 1:numIntKnots) points(intKnots[k],0,pch = 18,cex = 2,col = "darkmagenta") # add the interior points

### choose lambda
y <- WarsawApts$areaPerMzloty

fitSSauto <- smooth.spline(x,y) 
fitSSauto

fitGAMauto <- gam(y~s(x,bs = "cr",k = 22)) # k = number of basis and bs = "cr" specifies cubic regression spline
summary(fitGAMauto) # automatically fit the model with optimal lambda

### choose K
fitGAMdefault <- gam(y~s(x,bs = "cr"))
gam.check(fitGAMdefault) # k-index < 1 which is not good, and indicates larger number of basis

fitGAM52 <- gam(y~s(x,bs = "cr",k=52))
gam.check(fitGAM52) # now k-index > 1

fitGAM62 <- gam(y~s(x,bs = "cr",k=62))
gam.check(fitGAM62) # and keeping increasing k provides no improvement

fitGAM27<-gam(y~s(x,bs="cr",k=27))
gam.check(fitGAM27) # Actually, 27 is enough

### residual check
stdResids <- (residuals(fitGAM27)/sqrt(summary(fitGAM27)$scale))
par(mfrow = c(1,2),mai = c(1,1,0.1,0.1))
plot(x,stdResids,col="dodgerblue",cex = 1.5,
         xlab = "construction date (year)",
         ylab = "standardized residuals",
         cex.lab = 2,cex.axis = 2,bty="l",
         ylim=c(-3,max(stdResids)))
abline(h = 0,col = "slateblue")
for (hVal in c(-3,-2,2,3)) abline(h=hVal,lty = 2,col = "slateblue")
acf(stdResids[order(x)],col="darkgreen",cex.lab = 2,cex.axis = 2,main="")

### figure 1: normality, outliers
### figure 2: no sequential dependence
