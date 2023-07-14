################################################################################
### Example 1: Gauss Process Classification
library(kernlab) # kernel method lab
data(iris)

train = sample(1:150,100,replace = F)
GPC <- gausspr(Species~.,data=iris[train,],var=2)
  # type: regression or classification
  # kernel: default rbfdot (Gaussian)
  # kpar: the list of hyper-parameters, using MMLE to estimate them
  # var: initial sigma_y^2
GPC

# class probabilities 
predict(GPC, iris[-train,], type="probabilities")

################################################################################
### Example 2: Gauss Process Regression
library(fpp2)
data(uschange)
autoplot(uschange, facets=TRUE) +
  xlab("Year") + ylab("") +
  ggtitle("Quarterly percentage changes")

y.train = uschange[2:150,1]
x.train= uschange[1:149,2:5]
GPR <- gausspr(y.train~x.train,var=0.5)
GPR

# fitting
y.fit = predict(GPR,x.train)
Time = seq(1970.25,2007.25,0.25)
plot(Time, y.train, type ="l")
lines(Time, y.fit, col="red")
