library(faraway)
data(meatspec) # 1-172: training set
               # 173-215: testing set

###########################################################################
## collinearity and overfitting
model1 <- lm(fat  ~ ., meatspec[1:172,])
summary(model1)$r.squared   # good fitting

library(car)
vif(model1) # severe multicollinearity, rule of thumb: 10

rmse <- function(x,y){sqrt(mean((x-y)^2))} # function to calculate root of mean squared error (RMSE)
rmse(model1$fit, meatspec$fat[1:172]) # in-sample fitting RMSE
rmse(predict(model1, meatspec[173:215,]), meatspec$fat[173:215]) # out-of-sample prediction RMSE

## Note: prediction RMSE is much greater than in-sample fitting RMSE, which means model1 tends to be
## an over-fitting. 

##########################################################################
## principle component regression (PCR)
library(MVA)
meatpca <- prcomp(meatspec[1:172, -101]) # the 101st is fat, so remove it when doling PCA
round(meatpca$sdev/sum(meatpca$sdev),4) # the information contained in principal components

model2 <- lm(fat ~ meatpca$x[,1:4], meatspec[1:172,]) # use first 4 PCs to fit a model
rmse(model2$fit, meatspec$fat[1:172]) # in-sample fitting RMSE

rmsmeat <- numeric(18) # information mainly contained in first 18 PCs
mm <- apply(meatspec[1:172, -101], 2, mean) # calculate mean of each variable
tx <- as.matrix(sweep(meatspec[173:215, -101], 2, mm))
for (i in 1:18) {
   nx <- tx %*% meatpca$rot[,1:i]
   mode13 <- lm(fat~ meatpca$x[,1:i], meatspec[1:172,])
   pv <- cbind(1, nx) %*% mode13$coef
   rmsmeat[i] <- rmse(pv, meatspec$fat[173:215] )
   } 
plot(rmsmeat, ylab="Test RMS")
abline(h=rmse(predict(model1, meatspec[173:215,]), meatspec$fat[173:215]), col = "red")

################################################################################
## ridge regression
library(MASS)
trainx <- as.matrix(sweep(meatspec [1:172, -101], 2, mm))
y.centered <- meatspec$fat[1:172]-mean(meatspec$fat[1:172])
gridge <- lm.ridge(y.centered~trainx, lambda=seq(0,5e-8,1e-9)) # lambda is the tuning parameter
matplot(gridge$lambda, t(gridge$coef), type="l",lty=1,
         xlab=expression(lambda), ylab=expression(hat(beta))) # traceplot of beta 

select(gridge) # use GCV score which we discussed in the lecture of LASSO CV
which.min(gridge$GCV) # the 19th
ypredg <- scale(trainx, center=FALSE,
                scale=gridge$scales) %*% gridge$coef[, 19] + mean(meatspec$fat[1:172])
rmse(ypredg,meatspec$fat[1:172]) # in-sample fitting

mmt <- apply(meatspec[173:215, -101], 2, mean) 
testx <- as.matrix(sweep(meatspec[173:215, -101], 2, mmt)) # in test sample
ytpredg <- scale(testx, center=FALSE, scale
                  =gridge$scales)%*% gridge$coef [, 19] + mean(meatspec$fat[173:215])
rmse(ytpredg, meatspec$fat[173:215]) # poor prediction



