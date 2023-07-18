####################################################
library(GLMsData)

### Example 1: Poisson Regression
data(nminer)
# $Miners: number of Noisy Miners黑头矿鸟
# $Eucs: number of eucalypt trees桉树

## Simple Poisson Regression
PR1 <-  glm(Minerab ~ Eucs, data=nminer, family=poisson)
summary(PR1)

## Diagnostics (use quantile residual because response is discrete)
qr <- statmod::qresid(PR1)
qqnorm(qr, las=1);qqline(qr) # Q-Q plot
plot(qr ~ sqrt(fitted(PR1)),las=1);abline(h=0) # residual against fitted value: check systematic component
plot(qr[2:length(qr)]~qr[1:(length(qr)-1)],las=1);abline(h=0) # residual against lagged residual: check independence
plot(cooks.distance(PR1), type="h", las=1) # detect influential observations: 17th

maxinfl <- which.max(cooks.distance(PR1))
PR2 <- update(PR1,subset=(-maxinfl)) # remove the most influential observation
c( "Original model"=coef(PR1), "Without Infl"=coef(PR2))

plot( Minerab ~ jitter(Eucs), data=nminer,
      xlab="Number of eucalypts", ylab="Number of noisy miners")
newE <- seq( 0, 35, length=100)
newM1 <- predict( PR1, newdata=data.frame(Eucs=newE), type="response")
newM2 <- predict( PR2, newdata=data.frame(Eucs=newE), type="response")
lines( newM1 ~ newE, lty=1); lines( newM2 ~ newE, lty=2) 
# observation with large Eucs contributes significantly to the model
# so need to be very careful when deal with dispersion

####################################################

### Example 2: Modeling Rates rather than counts
# When upper bound of Y is large, modeling rates is preferred.
data(danishlc) # cases of lung cancer

## descriptive statistics
hist(danishlc$Cases)

danishlc$Rate <- danishlc$Cases/danishlc$Pop*1000 # Rate per 1000: number of cases of lung cancer among 1000 people
danishlc$Age <- ordered(danishlc$Age, # Ensure age-order is preserved
                          levels=c("40-54", "55-59", "60-64", "65-69", "70-74", ">74"))

danishlc$City <- abbreviate(danishlc$City, 1) # Abbreviate city names
matplot(xtabs(Rate ~ Age+City, data=danishlc), pch=1:4, lty=1:4,
           type="b", lwd=2, col="black", axes=FALSE, ylim=c(0, 25),
           xlab="Age group", ylab="Cases/1000")
axis(side=1, at=1:6, labels=levels(danishlc$Age))
axis(side=2, las=1); box()
legend("topleft", col="black", pch=1:4, lwd=2, lty=1:4, merge=FALSE,
         legend=c("Fredericia", "Horsens", "Kolding", "Vejle") )

# suggest age may provide a positive effect

## What if model counts directly?
LC.PR0 <- glm(Cases ~ log(Pop) + City*Age,
              family=poisson, data=danishlc)
summary(LC.PR0) # nothing significant

## Poisson regression modeling rate by setting log population as offset
## Rate_i = Cases_i_i/Pop_i, E[Rate_i] = mu_i/Pop_i
## eta_i = log(mu_i/Pop_i) = log(mu_i) - log(Pop_i)
## Hence, log(mu_i) = log(Pop_i) + eta_i
LC.PR1 <- glm(Cases ~ offset(log(Pop)) + City*Age,
               family=poisson, data=danishlc)
summary(LC.PR1)

####################################################

### Example 3: Two-way contingency table

## about the data: survey on attitude on genetically modified food among low- 
## and high-income households among Australia
Counts <- c(263, 258, 151, 222)
Att <- gl(2, 2, 4, labels=c("For", "Against"))
Inc <- gl(2, 1, 4, labels=c("High", "Low"))
data.frame(Counts, Att, Inc)
# gl() is used to generate level data, 2 means generate 2 levels and 
## the second 2 means repeat them 2 times, 4 means the length of the factor
gm.table <- xtabs(Counts ~Att+Inc); gm.table

## No Marginal Totals Are Fixed
GLM1 <- glm(Counts ~ Att + Inc, family=poisson)
anova(GLM1, test="Chisq")
coef(GLM1)

# with interaction (Att and Inc are associated)
GLM2 <- glm(Counts ~ Att*Inc, family=poisson)
anova(GLM2, test="Chisq")

#####################################################

### Example 4: Sampling zero and Structural zero
# male breast cancer乳腺癌: sampling zero
# females prostate cancer前列腺癌: structural zero
# males cervical cancer宫颈癌: structural zero
data(wacancer)
cancer.table <- xtabs(Counts ~ Cancer + Gender, data = wacancer); cancer.table 

Cancer1 <- glm(Counts ~ Cancer*Gender, data=wacancer, family=poisson)
anova(Cancer1,test = "Chisq") # information of gender absorbed by gender:cancer

# remove structural zero
wc <- subset(wacancer, (Cancer!="Breast"))
wc <- subset(wc, !(Cancer=="Cervix" & Gender=="M"))
wc <- subset(wc, !(Cancer=="Prostate" & Gender=="F"))
Cancer2 <- glm(Counts ~ Cancer*Gender, data=wc, family=poisson)
anova(Cancer2)

######################################################

### Example 5: over-dispersion
data(pock)
# large counts of pock mark represent severe viral activity

## descriptive statistics
plot(Count ~ jitter(log2(Dilution)), data=pock, las=1,
      xlab="Log (base 2) of dilution", ylab="Pock mark count")
mn <- with(pock, tapply(Count, log2(Dilution), mean) ) # Group means
vr <- with(pock, tapply(Count, log2(Dilution), var) ) # Group variances
plot( log(vr) ~ log(mn), las=1,
        xlab="Group mean", ylab="Group variance") 

# Group mean much less than group variance, indicating over-dispersion

m1 <- glm( Count ~ log2(Dilution), data=pock, family=poisson )
X2 <- sum(residuals(m1, type="pearson")^2)
c(Df=df.residual(m1), Pearson.X2=X2)
qchisq(0.975,df.residual(m1)) # suggest over-dispersion

## Negative Binomial GLM
m.nb <- MASS::glm.nb(Count ~ log2(Dilution), data=pock)
m.nb <- MASS::glm.convert(m.nb)
X2.nb <- sum(residuals(m.nb, type="pearson")^2) # provide a better fit
X2.nb 

## quasi-Poisson
m.qp <- glm(Count ~ log2(Dilution), data=pock, family="quasipoisson")
X2.qp <- sum(residuals(m.qp, type="pearson")^2) # the same Pearson chisq 
X2.qp

## compare the three model
coef.mat <- rbind( coef(m1), coef(m.qp), coef(m.nb) )
rownames(coef.mat) <- c("Poisson glm", "Quasi-Poisson", "Neg bin glm")
coef.mat

## difference between Poisson and quasi-Poisson
se.m1 <- coef(summary(m1))[, "Std. Error"]
se.qp <- coef(summary(m.qp))[, "Std. Error"]
data.frame(SE.Pois=se.m1, SE.Quasi=se.qp, ratio=se.qp/se.m1)
sqrt(summary(m.qp)$dispersion)