####################################################
library(GLMsData)

### Example 1: logit link (interpretation with odds ratio)
data(turbines) # 涡轮数据: 涡轮数量，涡轮工作总时长，涡轮裂缝总数

# want to investigate how working hours contribute to the growth of fissures
GLM.logit <- glm( Fissures/Turbines ~ Hours, family=binomial(link = "logit"),
               weights=Turbines, data=turbines) # Turbines is m_i
summary(GLM.logit)

# interpret the coefficient: one more working hour will 
# increases the odds of a turbine developing fissures by exp(0.0009992)

LogOdds <- predict(GLM.logit); odds <- exp(LogOdds)
plot(LogOdds ~ turbines$Hours, type="l", las=1,
        xlim=c(0, 5000), ylim=c(-5, 1),
        ylab="Log-odds", xlab="Run-time (in hours)")
my <- turbines$Fissures; m <- turbines$Turbines
EmpiricalOdds <- (my + 0.5)/(m - my + 0.5) # To avoid log of zeros, this is Bayes estimator with beta(0.5,0.5) prior
points( log(EmpiricalOdds) ~ turbines$Hours)
plot(odds ~ turbines$Hours, las=1, xlim=c(0, 5000), ylim=c(0, 2),
          type="l", ylab="Odds", xlab="Run-time (in hours)")
points( EmpiricalOdds ~ turbines$Hours)

### Example 2: probit link (interpretation with tolerance)
GLM.probit <- glm( Fissures/Turbines ~ Hours, family=binomial(link = "probit"),
                  weights=Turbines, data=turbines) # Turbines is m_i
summary(GLM.probit)

# interpret the coefficient: more working hour means the turbines more close to the tolerance level
Beyond_tolerance <- predict(GLM.probit,type = "response")
plot(Beyond_tolerance ~ turbines$Hours, type="l", las=1,
     xlim=c(0, 5000), ylim=c(0, 1),
     ylab="Rates of Beyond Tolerance", xlab="Run-time (in hours)")
my <- turbines$Fissures; m <- turbines$Turbines
mu_MLE <- my/m 
points( mu_MLE ~ turbines$Hours)

### Example 3: log-log link
data("mammary") 
# 目的：测量某样本中干细胞的比例
# 实验：将样本中的细胞移植到小鼠身上，若观测到小鼠身上被移植区域有增生的现象，
# 则说明移植到小鼠身上的细胞中至少包含一个干细胞，记录该被移植区域结果为阳性，否则为阴性
# 数据分析：根据实验记录（Bernoulli），用log-log link Binomial GLM估计样本中干细胞的比例
y <- mammary$N.Outgrowths/mammary$N.Assays
GLM.loglog <- glm(y~offset(log(N.Cells)), family=binomial(link="cloglog"),
           weights=N.Assays, data=mammary)
summary(GLM.loglog)

hat_lambda = exp(coef(GLM.loglog)) 
frequency_stem = 1/hat_lambda 
# interpret: The mammary stem cell frequency is estimated to be about 1 in 64 cells 

# confidence interval
Estimate <- summary(GLM.loglog)$coef[, "Estimate"]
SE <- summary(GLM.loglog)$coef[, "Std. Error"]
z <- qnorm(0.05/2, lower.tail=FALSE) # known dispersion parameter, use Z-interval
CI <- c(Lower=Estimate+z*SE, Estimate=Estimate, Upper=Estimate-z*SE)
CI <- 1/exp(CI)
round(CI, digits=1)

### Example 4: over-dispersion
data(germ)
GLM1 <- glm(Germ/Total ~ Extract * Seeds, family=binomial,
              weights=Total, data=germ )
sum(resid(GLM1, type="pearson")^2)
qchisq(0.975,df.residual(GLM1)) # suggest overdispersion
GLM2 <- update(GLM1, family=quasibinomial) # use quasi-binomial instead

### Example 5:  Hauck–Donner effect (JASA 1977)
# If mu near zero or one, linear predictor and se diverge and 
# since se diverges faster than beta, Wald statistic converges to zero
# Wald test is no longer reliable but can use likelihood ratio test or score test
data(nminer)
Eucs20 <- nminer$Eucs>20
mean(Eucs20) # small success prob
m1 <- glm(Miners ~ Eucs20, data=nminer, family=binomial(link = "probit")) # The issue is more severe with logit link
printCoefmat(coef(summary(m1))) # large coef, danger! Wald test fails!
anova(m1, test="Chisq") # likelihood ratio test still works

