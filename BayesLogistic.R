#### Generate data
DGP <- function(n,p,s){
  ## n: sample size
  ## p: dimension of parameter
  ## s: number of nonzero coefficients
  
  x = matrix(rnorm(n*p),n,p) ## generate x from iid standard normal
  beta = c(rep(sqrt(4*log(p/s)),s),rep(0,(p-s))) ## set-up coefficients
  succ_prob = 1/(1+exp(-x%*%beta)) ## calculate success probability
  y = rbinom(n,1,succ_prob) ## generate binary response
  
  return(list(y=y,x=x,beta=beta))
}

set.seed(2023)
n = 1000
p = 100
s = 10
simdata = DGP(n,p,s)

### glmnet
library(glmnet)
model.LASSO = glmnet(simdata$x,simdata$y,family = "binomial") # LASSO
coef(model.LASSO, s = s/p)
model.elasticnet = glmnet(simdata$x,simdata$y,alpha=0.5,family = "binomial") # elastic net
coef(model.elasticnet, s = s/p)

### Gelman (2008) insights
Gelman.sampling.beta <- function(y,x,beta0,ss){
  L = length(beta0)
  beta1 = numeric(L)
  accept = numeric(L)
  Z = (y-1)%*%x
  for (j in 1:L) {
    beta.star = beta0[j] + rnorm(1,0,ss[j])
    T1 = log(1+(beta0[j])^2)-log(1+beta.star^2) # prior
    C1 = sum(x[,j]*beta.star)<=-5; C2 = sum(x[,j]*beta0[j])<=-5
    if( C1 & !C2){
      T2 = Z[j]*(beta.star-beta0[j]) +sum(x[,j]*beta.star)-log(1+exp(-sum(x[,j]*beta0[j])))# likelihood
    }
    if(!C1 & C2){
      T2 = Z[j]*(beta.star-beta0[j]) - log(1+exp(-sum(x[,j]*beta.star)))-sum(x[,j]*beta0[j])# likelihood
    }
    if(C1 & C2){
      T2 = Z[j]*(beta.star-beta0[j]) +sum(x[,j]*beta.star)-sum(x[,j]*beta0[j])# likelihood
    }
    if(!C1 & !C2){    
      T2 = Z[j]*(beta.star-beta0[j]) - log(1+exp(-sum(x[,j]*beta.star)))+log(1+exp(-sum(x[,j]*beta0[j])))# likelihood
    }
    MH = min(1,exp(T1+T2))
    U = runif(1,0,1)
    if(U<=MH){beta1[j]=beta.star;accept[j] = 1}
    else{beta1[j]=beta0[j];accept[j] = 0}
  }
  return(list(beta = beta1,accept = accept))
}

gelman08 <- function(y,x,beta0,ss_beta,batch_size,max_batch,N,n_burn,n_thin){
  ## elementwise adaptive random-walk Metropolis-Hastings
  ## y: binary response
  ## x: covariates
  ## beta0: initial value
  ## ss: initial step size
  ## batch_size: batch size for tuning
  ## max_batch: maximal number of batches for tuning
  ## N: length of chain
  ## n_burn: length of burn-in samples
  ## n_thin: number of thinning
  
  L = length(beta0)
  betas = matrix(0,L,(N-n_burn)/n_thin)
  
  ## tuning procedure
  if (max_batch > 0) {
    cat("tuning begins...\n")
    for (index_batch in 1:max_batch) {
      ar_beta = rep(0, L)
      for (index_iter in 1:batch_size) {
        beta_res <- Gelman.sampling.beta(y,x,beta0,ss_beta)
        beta0 <- beta_res$beta
        ar_beta <- ar_beta + beta_res$accept
      }
      ar_beta <- ar_beta / batch_size
      for (j in 1:L) {
        if (ar_beta[j] < 0.35) ss_beta[j] <- ss_beta[j] * 0.9
        else if (ar_beta[j] > 0.45) ss_beta[j] <- ss_beta[j] * 1.1
      }
      cat("batch", index_batch, "ar_beta1:", ar_beta[1], "ss_beta1:", ss_beta[1], "ar_beta610:", ar_beta[100], "ss_beta610", ss_beta[100], "\n")
    }
    cat("tuning ends\n")
  }
  
  ## sampling procedure
  ar_beta <- rep(0, L)
  for (i in 1:N) {
    if (i %% 1000 == 0) cat(i, "\n")
    beta_res <- Gelman.sampling.beta(y,x,beta0,ss_beta)
    beta0 <- beta_res$beta
    ar_beta <- ar_beta + beta_res$accept
    
    if (i > n_burn & ((i - n_burn)%%n_thin == 0)) {
      betas[,(i-n_burn)/n_thin] = beta0
    }
  }

  ar_beta <- ar_beta / N
  
  samples = list(beta = betas, ar_beta = ar_beta)
  return(samples)
}

beta0 = rep(0.5,p)
ss_beta = rep(1,p)
batch_size = 200
max_batch = 20
N = 100000
n_burn = 50000
n_thin = 10
set.seed(2023)
model.gelman = gelman08(simdata$y,simdata$x,
                        beta0,ss_beta,
                        batch_size,max_batch,
                        N,n_burn,n_thin)

plot(model.gelman$ar_beta)
plot(model.gelman$beta[1,],type = "l")
plot(model.gelman$beta[100,],type = "l")

# plot posterior distribution
plot_post <- function(i,y,x,beta){
  Z = y%*%x
  beta0 = seq(-10,10,0.1);L = length(beta0)
  post = numeric(L)
  for (j in 1:L) {
    beta[i] = beta0[j]
    T1 = -log(1+(beta0[j])^2) # prior
    
    post[j] = exp(T1+T2)
  }
  plot(post~beta,type = "l",xlim = c(-10,10),xlab = expression(beta),ylab = "posterior kernel")
}

plot_post(1,y=simdata$y,x = simdata$x)

### Ghosh (2018) insights
library(mvtnorm)
library(truncnorm)
library(coda)
library(BayesLogit)

is.binary = function(x){
  length(unique(x)) == 2;
}

## this function is written by Professor Yingbo Li (2017)
tglm.fit = function(y, X, iter = 1e5, thin = max(1, round(iter/2e3)), burnin = 0.5, 
                    method = 'logistic', df = 1, slope.scale = 2.5, intercept.scale = 10,
                    save.latent = FALSE, center.binary = TRUE, scale.continuous = TRUE, 
                    beta.original = TRUE, track.time = TRUE, show.summary = TRUE){
  
  t.start = proc.time();
  #### normalizing X ######
  n = dim(X)[1];
  p = dim(X)[2];
  Xmean = apply(X, 2, mean);
  Xsd = apply(X, 2, sd);
  Xbinary = apply(X, 2, is.binary);
  
  ## center binary predictors
  if(center.binary == TRUE)
    X[, Xbinary] = sweep(as.matrix(X[, Xbinary]), 2, Xmean[Xbinary]);
  ## center and scale continuous predictors
  if(scale.continuous == TRUE){
    X[, !Xbinary] = sweep(as.matrix(X[, !Xbinary]), 2, Xmean[!Xbinary]);
    X[, !Xbinary] = sweep(as.matrix(X[, !Xbinary]), 2, Xsd[!Xbinary] * 2, FUN = "/");
  }
  
  ## include intercept column (all 1's) in the design matrix
  X = cbind(1, X);
  
  ## sigma: a (p + 1)-vector of prior scales
  if(length(slope.scale) == 1)
    sigma = c(intercept.scale, rep(slope.scale, p));
  if(length(slope.scale) == p)
    sigma = c(intercept.scale, slope.scale);
  if(length(slope.scale) != p && length(slope.scale) != 1)
    stop('Error: length of slope.scale should equal 1 or number of columns of X!');
  
  
  #### intial values of beta (including intercept), gamma (from their priors) ####
  beta = rnorm(p + 1, 0, sigma);
  if(df < Inf)
    gamma = 1 / rgamma(p + 1, df / 2, rate = sigma^2 * df / 2);
  if(df == Inf)
    gamma = sigma^2;
  if(method == 'logistic')
    omega = rpg(n, 1, X %*% beta);
  if(method == 'probit'){
    ## upper bounds and lower bounds of z (depend on y)
    z.lower = z.upper = rep(0, n);
    z.lower[y == 0] = -Inf;
    z.upper[y == 1] = Inf;
    z = rtruncnorm(n, a = z.lower, b = z.upper, mean = X %*% beta, sd = 1);  	
  }
  
  #### Save MCMC chains ####
  Beta = matrix(NA, nrow = iter / thin, ncol = p + 1);
  if(df < Inf)
    Gamma = matrix(NA, nrow = iter / thin, ncol = p + 1);
  if(save.latent == TRUE){
    if(method == 'logistic')
      Omega = matrix(NA, nrow = iter / thin, ncol = n);
    if(method == 'probit')
      Z = matrix(NA, nrow = iter / thin, ncol = n);
  }
  
  #### Gibbs sampler ####
  for(it in 1:iter){
    
    if(method == 'logistic'){    
      ## update beta ##
      cov.beta = solve(t(X) %*% sweep(X, 1, omega, FUN = "*") + diag(1 / gamma));
      beta = c(rmvnorm(1, cov.beta %*% (t(X) %*% (y - 0.5)), cov.beta));
      ## update omega ##
      omega = rpg(n, 1, X %*% beta);
    }
    
    if(method == 'probit'){    
      ## update beta ##
      cov.beta = solve(t(X) %*% X + diag(1 / gamma));
      beta = c(rmvnorm(1, cov.beta %*% (t(X) %*% z), cov.beta));
      ## update z ##
      z = rtruncnorm(n, a = z.lower, b = z.upper, mean = X %*% beta, sd = 1);
    }
    
    ## update gamma ##
    if(df < Inf)
      gamma = 1 / rgamma(p + 1, (df + 1) / 2, rate = (beta^2 + sigma^2 * df) / 2);
    
    ## save every thin iterations
    if(it %% thin == 0){
      record = it / thin;
      Beta[record, ] = beta;
      if(df < Inf)
        Gamma[record, ] = gamma;  	    
      if(save.latent == TRUE){
        if(method == 'logistic')
          Omega[record, ] = omega;
        if(method == 'probit')
          Z[record, ] = z;
      }
      
      
      ## show: x0% completed
      if( (it * 10) %% iter == 0 && it != iter && track.time == TRUE)
        cat(paste( (it * 100) / iter), '% completed...\n', sep = '');
      if( it == iter && track.time == TRUE){
        cat(paste( (it * 100) / iter), '% completed.\n', sep = '');
        cat('\n');
      } 
      
    }
    
  }
  
  ## track time
  if(track.time == TRUE){
    t.finish = proc.time();
    cat('Time used (in second): \n')
    print(t.finish - t.start);
    cat('\n');
  }
  
  ## transform beta back to orginal scale (before centering and scaling)
  if(center.binary == TRUE && beta.original == TRUE)
    Beta[, 1] = Beta[, 1] - as.matrix(Beta[, which(Xbinary) + 1]) %*% Xmean[Xbinary];
  if(scale.continuous == TRUE && beta.original == TRUE){
    Beta[, 1] = Beta[, 1] - as.matrix(Beta[, which(!Xbinary) + 1]) %*% (Xmean[!Xbinary] / 2 / Xsd[!Xbinary]);
    Beta[, which(!Xbinary) + 1] = as.matrix(Beta[, which(!Xbinary) + 1]) / 2 / Xsd[!Xbinary];
  }
  
  ##########################################
  ######  Posterior Inference:        ######
  ######  mean, sd, HDP               ######
  ##########################################
  
  keep = round(iter / thin * burnin): (iter / thin);
  varnames = colnames(X[, -1]);
  if(length(varnames) == 0)
    varnames = paste('X', 1:p, sep = '');
  
  inference = matrix(NA, ncol = 4, nrow = p + 1);
  rownames(inference) = c('intercept', varnames);
  colnames(inference) = c('Est', 'Std', '95HPDlower','95HPDupper');
  
  ## posterior mean and HPD
  inference[, 1] = apply(as.matrix(Beta[keep, ]), 2, mean);
  inference[, 2] = apply(as.matrix(Beta[keep, ]), 2, sd);
  inference[, 3:4] = HPDinterval(as.mcmc(Beta[keep, ]), 0.95);
  
  ## show summary: inference
  if(show.summary == TRUE){
    cat('Posterior inference: \n');
    print(round(inference, 4));
  }
  
  
  if(save.latent == FALSE && df < Inf)
    return(list(Beta = Beta, Gamma = Gamma, inference = inference));
  if(save.latent == FALSE && df == Inf)
    return(list(Beta = Beta, inference = inference));
  if(save.latent == TRUE && df < Inf){
    if(method == 'logistic')
      return(list(Beta = Beta, Gamma = Gamma, Z = Omega, inference = inference));
    if(method == 'probit')
      return(list(Beta = Beta, Gamma = Gamma, Z = Z, inference = inference));
  }
  if(save.latent == TRUE && df == Inf){
    if(method == 'logistic')
      return(list(Beta = Beta, Z = Omega, inference = inference));
    if(method == 'probit')
      return(list(Beta = Beta, Z = Z, inference = inference));
  }
  
}

model.Li <- tglm.fit(y = simdata$y, X = simdata$x) # Cauchy prior, 21 minutes to run
save(model.Li,file = "Model_Li.Rdata")
load("Model_Li.Rdata")

model.Li.normal <- tglm.fit(y = simdata$y, X = simdata$x, df = Inf) #Normal prior, recommended by Ghosh et al
save(model.Li.normal,file = "Model_Li_normal.Rdata")
load("Model_Li_normal.Rdata")


## It's not the end of the story. It there a method for high-dim data that will not cause over-estimation?
## Piironen and Vehtari 2017, regularized Horseshoe prior
