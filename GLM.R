####################################################
### Part 1: model estimation
data(trees)
GLM1 <- glm( Volume ~ log(Height) + log(Girth), data=trees,
                  family=Gamma(link="log"))

## estimation of coefficients
coef(GLM1)

## Pearson estimator of dispersion parameter
w <- weights(GLM1, type="working")
e <- residuals(GLM1, type="working")
sum(w*e^2)/df.residual(GLM1)

summary(GLM1)$dispersion # R default reports Pearson estimator

## investigate the iteration process
GLM2 <- update(GLM1, control=glm.control(
  maxit=3, # Max of 3 iterations
  epsilon=1e-15, # Stopping criterion
  trace=TRUE) )

GLM2 <- update(GLM1, control=glm.control(
  maxit=30, # Max of 3 iterations
  epsilon=1e-15, # Stopping criterion
  trace=TRUE) )

## summary of glm object
summary(GLM1)

####################################################
### Part 2: model diagnostics

## Check independence: plot residuals against lagged residuals
plot(resid(GLM1)[2:31]~resid(GLM1)[1:30]);abline(h=0) # no special pattern
# for discrete observations, use statmod::qresid()

## Check systematic component: plot residuals against fitted value
plot(resid(GLM1)~ sqrt(fitted(GLM1)),las=1);abline(h=0) # no special pattern

## Check random component: Q-Q plot of quantile residuals
qqnorm(statmod::qresid(GLM1), las=1); qqline(statmod::qresid(GLM1)) # no obvious violation

## Outliers and Influential Observations
rs <- cbind(rD=resid(GLM1), standardrD=rstandard(GLM1),
            studentizedr=rstudent(GLM1), rQ=statmod::qresid(GLM1))
apply(abs(rs), 2, max) # compared with mean(trees$Volume)=30.17, no large residuals

im <- influence.measures(GLM1)
im$infmat # influential measures
colSums(im$is.inf) # count the number of influential observations by every measure

# if remove the most influential one by Cook
plot(im$infmat[,6], type="h", ylab="Cook's distance", las=1)
infl <- which.max(im$infmat[,6])
GLM.infl <- update(GLM1, subset=(-infl)) # refit the model
AIC(GLM.infl)>AIC(GLM1)
BIC(GLM.infl)>BIC(GLM1)

####################################################
### Part 3: testing and confidence interval

## Wald test for individual parameter
printCoefmat(coef(summary(GLM1)))
confint(GLM1,level = 0.95) # confidence interval of coefficients

## Confidence interval of fitted value
fitting <- predict(GLM1,
                se.fit=TRUE, # also return the standard error of fitted value
                type="response"  # to fit the response, without it R fits linear predictor
                )

zstar <- qt(0.975,df.residual(GLM1)) # consider 95% interval
L <- fitting$fit - zstar*fitting$se.fit
U <- fitting$fit + zstar*fitting$se.fit

plot(fitting$fit~seq(1,31,1),type = "l",col = "red",lty = 1,
     xlim = c(0,31),ylim = c(10,90),xlab = "ID", ylab = "Volume")
par(new = T)
plot(trees$Volume~seq(1,31,1),type = "l",col = "black",lty = 1,
     xlim = c(0,31),ylim = c(10,90),xlab = "",ylab = "")
par(new = T)
plot(L~seq(1,31,1),type = "l",col = "blue",lty = 2,
     xlim = c(0,31),ylim = c(10,90),xlab = "", ylab = "")
par(new = T)
plot(U~seq(1,31,1),type = "l",col = "blue",lty = 2,
     xlim = c(0,31),ylim = c(10,90),xlab = "", ylab = "")
legend("topleft",legend = c("Fitted Value","Observation","95% CI"),
       col = c("red","black","blue"),lty = c(1,1,2))

## Confidence interval for predicted value
# If we have two new observation, Girth = c(9,10), Height = c(77,81)
prediction <- predict(GLM1, 
                      newdata=data.frame(Girth = c(9,10), Height = c(77,81)), 
                      se.fit=TRUE,
                      type="response")

Lhat <- prediction$fit - zstar*prediction$se.fit
Uhat <- prediction$fit + zstar*prediction$se.fit

## ANOVA and likelihood ratio test
GLM0 <- glm( Volume ~ 1, data=trees,
             family=Gamma(link="log"))
L <- deviance(GLM0) - deviance(GLM1) # test statistic
pchisq(L, df.residual(GLM0) - df.residual(GLM1), lower.tail=FALSE) # suggest GLM1

anova(GLM1,test = "F") # how explanatory variables reduce deviance

## Score test (for example, whether to include Girth as explanatory variable)
GLM_1 <- GLM1 <- glm( Volume ~ log(Height), data=trees,
                      family=Gamma(link="log"))
z.stat <- statmod::glm.scoretest(GLM_1, trees$Girth)
p.val <- 2*pnorm(abs(z.stat), lower.tail=FALSE) # yes, should include Girth

## stepwise model selection
min.model <- glm( Volume~1, data=trees, family=Gamma(link="log"))
max.model <- glm( Volume~log(Girth) + log(Height),
                    data=trees, family=Gamma(link="log"))
m.step <- step(min.model, scope=list(lower=min.model, upper=max.model),
             direction="both") # direction can also be forward or backward

