########################################################
### Take the lady tasting tee problem as example to illustrate different views of testing
### Let n be the number of cups
### Let t be the number of cups the lady successfully claimed
### Let theta be the probability of success
n = 10

### Fisher
## H0: theta = 1/2
t = seq(1,10,1)
pvalue = pbinom(t,10,1/2,lower.tail = F)
plot(pvalue~t,type = "l")
abline(h = 0.05,lty = 2, col = "red") # to reject H0 under high standard, the lady should at least tell 8 cups correctly 
abline(h = 0.01,lty = 2, col = "red") # to reject H0 under high standard, the lady should at least tell 9 cups correctly 

### Neyman-Pearson
## H0: theta = 1/2
## H1: theta not equals 1/2
critvalue = seq(6,10,1)
theta = seq(0.01,1,0.01)
for (k in 1:5) {
  plot(pbinom(critvalue[k],10,theta,lower.tail = F)~theta,type = "l",lty = k,xlim = c(0,1),ylim = c(0,1),ylab = "")
  par(new = T)
}
abline(v = 1/2)
abline(h = 0.05, col = "red") # assume alpha<=0.05
legend("topleft",legend = c("crit=6","crit=7","crit=8","crit=9","crit=10"),lty = c(1,2,3,4,5))
## suggest critical value = 8, if t>=8, reject H0

### Jeffreys' Bayesian (Bayes factor)
BF = numeric(10)
for (i in 1:10) {
  BF[i] = dbinom(t[i],10,1/2)/dbinom(t[i],10,t[i]/10)
}
plot(BF~t,type = "l") 
# t = 10, decisive
# t = 9, very strong
# t = 8, substantial
# t = 7, not worth more than a bare mention

### Relative Belief
RB = numeric(10)
for (i in 1:10) {
  RB[i] = dbinom(t[i],10,1/2)/(0.5*dbinom(t[i],10,1/2) + 0.5*VGAM::dbetabinom.ab(t[i],10,0.5,0.5))
}
plot(RB~t,type = "l");abline(h = 1, col = "red")
# If RB<1, reject H0
# so t>=8 supports the Lady can tell it.


