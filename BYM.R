#####################################################################
### 1989-1993 Suicides data in 32 boroughs of London
### Model 1
### y_i \sim Poisson(E_i\rho_i)
### log(rho_i) = b_0 + u_i + v_i,u_i structured effect, v_i \sim N(0,sigma_v^2) unstructured random effect
### rho_i: standardized mortality ratio
### E_i: the number of expected cases of suicides for each borough
### b_0: the average suicide rate in all the 32 boroughs

### Model 2
### y_i \sim Poisson(E_i\rho_i)
### log(rho_i) = b_0 + \beta_1x_1 + \beta_2x_2 + u_i + v_i
### x1: index of social deprivation
### x2: index of social fragmentation
#####################################################################

library(maptools)
library(spdep)
library(sf)
library(INLA) # install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(lattice)

### load data
load("LondonSuicides.RData")

london.gen = readShapePoly("LDNSuicides.shp") # read polygon shape file or say .shp file
names<- sort(london.gen$NAME)
data <- data.frame(NAME=names, y=y, E=E, x1=x1, x2=x2) # attach name to each borough

Nareas <- length(data[,1])

boroughs=london.gen
data.boroughs=attr(boroughs, "data")
order <- match(data.boroughs$NAME,data$NAME)
data <- data[order,]
ID<-seq(1,32)
data <- cbind(ID,data)

### prepare the data: specify a graph which assigns the set of neighbors for each borough
temp <- poly2nb(london.gen)
nb2INLA("LDN.graph", temp)
LDN.adj <- paste(getwd(),"/LDN.graph",sep="") # create LDN.graph under work dictionary
H <- inla.read.graph(filename="LDN.graph") # data frame to store neighboring information
image(inla.graph2matrix(H),xlab="",ylab="")

### estimate the model
formula <- y ~ 1 + f(ID, model="bym",graph=LDN.adj) # 1 specifies intercept and f() specifies u_i + v_i

# or specify other values of log tau in precision (precision = tau*neighbor correlation)
formula <- y ~ 1 + f(ID, model="bym",graph=LDN.adj,
                     scale.model=TRUE,
                     hyper=list(prec.unstruct=
                                  list(prior="loggamma",param=c(1,0.001)),
                                prec.spatial=list(
                                  prior="loggamma",param=c(1,0.001))))

mod <- inla(formula,family="poisson",data=data,E=E)

#We calculate zeta=exp(csi) where csi=upsilon + nu
m <- mod$marginals.random$ID[1:Nareas]
zeta <- lapply(m,function(x)inla.emarginal(exp,x))

#We now calculate the probability that the spatial effects zeta are above 1, 
#identifying areas with excess risk of suicides. This is equivalent to 
#calculate the probability that csi is above 0, which is easier to obtain
a=0
inlaprob<-lapply(mod$marginals.random$ID[1:Nareas], function(X){
  1-inla.pmarginal(a, X)
})

#Finally we obtain the proportion of variance explained by the spatially structured component 
#upsilon, taking the structured effect upsilon and calculating the empirical variance. 
#First we create a matrix with rows equal to the number of areas and 1000 columns. 
#Then for each area we extract 1000 values from the corresponding marginal distribution of upsilon 
#and finally we calculate the empirical variance. We also extract the expected value of the variance 
#for the unstructured component and build the spatial fractional variance.
mat.marg<-matrix(NA, nrow=Nareas, ncol=1000)
m<-mod$marginals.random$ID
for (i in 1:Nareas){
  u<-m[[Nareas+i]]
  s<-inla.rmarginal(1000, u)
  mat.marg[i,]<-s}
var.RRspatial<-mean(apply(mat.marg, 2, sd))^2
var.RRhet<-inla.emarginal(function(x) 1/x,
                          mod$marginals.hyper$"Precision for ID (iid component)")
var.RRspatial/(var.RRspatial+var.RRhet)
########################
#Now we add the covariates (deprivation - x1 and social fragmentation - x2) and repeat the steps
########################
formula.cov <- y ~ 1+ f(ID, model="bym", graph=LDN.adj) + x1 + x2
mod.cov <- inla(formula.cov,family="poisson",data=data,E=E)
mod.cov$summary.fixed
m <- mod.cov$marginals.random$ID[1:Nareas]
zeta.cov <- lapply(m,function(x)inla.emarginal(exp,x))
a=0
inlaprob.cov<-lapply(mod.cov$marginals.random$ID[1:Nareas], function(X){
  1-inla.pmarginal(a, X)
})
m<-mod.cov$marginals.random$ID
mat.marg<-matrix(NA, nrow=Nareas, ncol=1000)
for (i in 1:Nareas){
  u<-m[[Nareas+i]]
  s<-inla.rmarginal(1000, u)
  mat.marg[i,]<-s}
var.RRspatial<-mean(apply(mat.marg, 2, sd))^2
var.RRhet<-inla.emarginal(function(x) 1/x,
                          mod.cov$marginals.hyper$"Precision for ID (iid component)")
var.RRspatial/(var.RRspatial+var.RRhet)

#Finally we build the maps. First we create a dataset with all the relevant quantities 
#and classes of SMR and posterior probabilities. Then transform the continuous SMR and 
#posterior probabilities in factors, Merge the spatial polygon of London boroughs with the 
#data and map the quantities.
Spatial.results<- data.frame(NAME=data$NAME,SMR=unlist(zeta),
                             pp=unlist(inlaprob), SMR.cov = unlist(zeta.cov), pp.cov = unlist(inlaprob.cov))
SMR.cutoff<- c(0.6, 0.9, 1.0, 1.1,  1.8)
pp.cutoff <- c(0,0.2,0.8,1)
#Transform SMR and pp in factors
SMR.DM=cut(Spatial.results$SMR,breaks=SMR.cutoff,include.lowest=TRUE)
pp.DM=cut(Spatial.results$pp,breaks=pp.cutoff,include.lowest=TRUE)
SMR.COV=cut(Spatial.results$SMR.cov,breaks=SMR.cutoff,include.lowest=TRUE)
pp.COV=cut(Spatial.results$pp.cov,breaks=pp.cutoff,include.lowest=TRUE)
maps.SMR.factors <- data.frame(NAME=data$NAME,
                               SMR.DM=SMR.DM,pp.DM=pp.DM,SMR.COV=SMR.COV,pp.COV=pp.COV)
attr(boroughs, "data")=merge(data.boroughs,maps.SMR.factors,by="NAME")
trellis.par.set(axis.line=list(col=NA))
spplot(obj=boroughs, zcol= "SMR.DM", col.regions=gray(3.5:0.5/4),main="")
trellis.par.set(axis.line=list(col=NA))
spplot(obj=boroughs, zcol= "pp.DM", col.regions=gray(2.5:0.5/3),main="")
trellis.par.set(axis.line=list(col=NA))
spplot(obj=boroughs, zcol= "SMR.COV", col.regions=gray(3.5:0.5/4))
trellis.par.set(axis.line=list(col=NA))
spplot(obj=boroughs, zcol= "pp.COV", col.regions=gray(2.5:0.5/3))