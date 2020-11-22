
load("dataex2.Rdata")

d = dataex2

log_like = function(y, mu){
  sum(y[,2]*dnorm(y[,1], mean = mu, 1.5, log = TRUE) + (1 - y[,2])*pnorm(y[,1], mean = mu, 1.5, log = TRUE))
}

require(maxLik)
mle <- maxLik(logLik = log_like, y = d, start = c(1))
summary(mle)





load("dataex4.Rdata")

EM = function(data, beta0, epsilon){
  data = data
  x = data[,1]
  y = data[,2]
  n = length(y)
  mis = which(is.na(y))
  y_obs = y[-mis]
  y_mis = y[mis]
  m = length(y_obs)
  diff = 0.1
  beta = beta0
  
  while (diff > epsilon) {
    beta.old = beta
    Q = function(pa){
      beta_0 = pa[1]
      beta_1 = pa[2]
      q_obs = sum(y_obs*(beta_0 + x[-mis]*beta_1)) - sum(log(1 + exp(beta_0 + x*beta_1)))
      q_mis = 0
      for (i in 1:(n-m)) {
        mi = mis[i]
        qt = (exp(beta.old[1] + x[mi]*beta.old[2])/(1 +  exp(beta.old[1] + x[mi]*beta.old[2])))*(beta_0 + x[mi]*beta_1)
        q_mis = q_mis + qt
        
      }
      obj = q_obs + q_mis
    }
    mle = optim(par = c(beta.old[1], beta.old[2]), fn = Q, control = list("fnscale"=-1), hessian = TRUE)
    beta = mle$par
    diff = abs(beta - beta.old)
  }
  return(beta)
}


EM(dataex4, c(0,0), 0.000001)


