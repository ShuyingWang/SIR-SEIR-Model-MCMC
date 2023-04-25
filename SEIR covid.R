### Description of the code file ###
# The first part defines the functions.
# The second part includes the illustration with COVID-19 data (Section 4.4 in the paper).
# Comment/uncomment the new/current MCMC algorithm to run comparisons.
# The last part computes the basic reproduction number by R0 package.


### 1. Functions ###

## Integrals (for calculating likelihoods)
SI_int <- function(t.x, t.y, t.z){
  t = sort(c(t.x, t.y, t.z, t0, T.end))
  n = length(t) - 1
  S = S0
  I = I0
  int = 0
  
  for (i in 1:n){
    int = int + (t[i+1] - t[i]) * S * I
    if (t[i+1] %in% t.x){
      S = S - 1
    }else if (t[i+1] %in% t.y){
      I = I + 1
    }else{
      I = I - 1
    }
  }
  return(int)
}
SI.int <- function(t.x, t.y, t.z, t1, t2){
  t = sort(c(t.x[t.x>t1 & t.x<t2], t.y[t.y>t1 & t.y<t2], t.z[t.z>t1 & t.z<t2], t1, t2))
  n = length(t) - 1
  S = S0 - sum(t.x<t1)
  I = I0 + sum(t.y<t1) - sum(t.z<t1)
  int = 0
  
  for (i in 1:n){
    int = int + (t[i+1] - t[i]) * S * I
    if (t[i+1] %in% t.x){
      S = S - 1
    }else if (t[i+1] %in% t.y){
      I = I + 1
    }else{
      I = I - 1
    }
  }
  return(int)
}
I_int <- function(t.y, t.z){
  t = sort(c(t.y, t.z, t0, T.end))
  n = length(t) - 1
  I = I0
  int = 0
  
  for (i in 1:n){
    int = int + (t[i+1] - t[i]) * I
    if (t[i+1] %in% t.y){
      I = I + 1
    }else{
      I = I - 1
    }
  }
  return(int)
}
I.int <- function(t.y, t.z, t1, t2){
  t = sort(c(t.y[t.y>t1 & t.y<t2], t.z[t.z>t1 & t.z<t2], t1, t2))
  n = length(t) - 1
  I = I0 + sum(t.y<t1) - sum(t.z<t1)
  int = 0
  
  for (i in 1:n){
    int = int + (t[i+1] - t[i]) * I
    if (t[i+1] %in% t.y){
      I = I + 1
    }else{
      I = I - 1
    }
  }
  return(int)
}
E_int <- function(t.x, t.y){
  t = sort(c(t.x, t.y, t0, T.end))
  n = length(t) - 1
  E = E0
  int = 0
  
  for (i in 1:n){
    int = int + (t[i+1] - t[i]) * E
    if (t[i+1] %in% t.x){
      E = E + 1
    }else{
      E = E - 1
    }
  }
  return(int)
}
E.int <- function(t.x, t.y, t1, t2){
  t = sort(c(t.x[t.x>t1 & t.x<t2], t.y[t.y>t1 & t.y<t2], t1, t2))
  n = length(t) - 1
  E = E0 + sum(t.x<t1) - sum(t.y<t1)
  int = 0
  
  for (i in 1:n){
    int = int + (t[i+1] - t[i]) * E
    if (t[i+1] %in% t.x){
      E = E + 1
    }else{
      E = E - 1
    }
  }
  return(int)
}


## Products (for calculating likelihoods)
Log.SI_prod <- function(t.x, t.y, t.z){
  times = sort(c(t.x, t.y, t.z))
  S = S0
  I = I0
  lp = 0
  
  for (t in times){
    if (t %in% t.x){
      lp = lp + log(S) + log(I)
      S = S - 1
    }else if (t %in% t.y){
      I = I + 1
    }else{
      I = I - 1
    }
  }
  return(lp)
}
Log.I_prod <- function(t.y, t.z){
  times = sort(c(t.y, t.z))
  I = I0
  lp = 0
  
  for (t in times){
    if (t %in% t.z){
      lp = lp + log(I)
      I = I - 1
    }else{
      I = I + 1
    }
  }
  return(lp)
}
Log.E_prod <- function(t.x, t.y){
  times = sort(c(t.x, t.y))
  E = E0
  lp = 0
  
  for (t in times){
    if (t %in% t.y){
      lp = lp + log(E)
      E = E - 1
    }else{
      E = E + 1
    }
  }
  return(lp)
}
Log.SI.prod <- function(t.x, t.y, t.z, t1, t2){
  times = sort(c(t.x[t.x>t1 & t.x<t2], t.y[t.y>t1 & t.y<t2], t.z[t.z>t1 & t.z<t2]))
  S = S0 - sum(t.x<t1)
  I = I0 + sum(t.y<t1) - sum(t.z<t1)
  lp = 0
  
  for (t in times){
    if (t %in% t.x){
      lp = lp + log(S) + log(I)
      S = S - 1
    }else if (t %in% t.y){
      I = I + 1
    }else{
      I = I - 1
    }
  }
  return(lp)
}
Log.I.prod <- function(t.y, t.z, t1, t2){
  times = sort(c(t.y[t.y>t1 & t.y<t2], t.z[t.z>t1 & t.z<t2]))
  I = I0 + sum(t.y<t1) - sum(t.z<t1)
  lp = 0
  
  for (t in times){
    if (t %in% t.z){
      lp = lp + log(I)
      I = I - 1
    }else{
      I = I + 1
    }
  }
  return(lp)
}
Log.E.prod <- function(t.x, t.y, t1, t2){
  times = sort(c(t.x[t.x>t1 & t.x<t2], t.y[t.y>t1 & t.y<t2]))
  E = E0 + sum(t.x<t1) - sum(t.y<t1)
  lp = 0
  
  for (t in times){
    if (t %in% t.y){
      lp = lp + log(E)
      E = E - 1
    }else{
      E = E + 1
    }
  }
  return(lp)
}
Indicator <- function(t.x, t.y, t.z){
  times = sort(c(t.x, t.y, t.z))
  S = S0
  E = E0
  I = I0
  ind = 1
  
  for (t in times){
    if (t %in% t.x){
      if (S<=0 | I<=0){
        ind = 0
      }
      S = S - 1
      E = E + 1
    }else if (t %in% t.y){
      if (E <= 0){
        ind = 0
      }
      E = E - 1
      I = I + 1
    }else{
      I = I - 1
    }
  }
  return(ind)
}


## New proposal for the latent process
New_proposal.x <- function(beta, nx, t.y, t.z, t1, t2){
  S = S0 - nx
  I = I0 + sum(t.y<t1) - sum(t.z<t1)
  t = t1
  t.x = c()
  t.y = c(t.y[t.y>t1 & t.y<t2], t2)
  t.z = c(t.z[t.z>t1 & t.z<t2], t2)
  i = 1
  j = 1
  
  while (t < t2){
    if (I > 0){
      t_new = t + rexp(1, beta/N * S * I)
      
      if (t_new < t.y[i] & t_new < t.z[j]){
        t = t_new
        t.x = append(t.x, t)
        S = S - 1
      }else if (t.y[i] < t.z[j]){
        t = t.y[i]
        I = I + 1
        i = i + 1
      }else{
        t = t.z[j]
        I = I - 1
        j = j + 1
      }
      
    }else{
      t = t.y[i]
      I = I + 1
      i = i + 1
    }
  }
  return(t.x)
}





### 2. Illustration with Covid Data ###
N = 285000
E0 = 15
I0 = 3
S0 = N - I0 - E0
t0 = 0
T.end = 35
k = 16
beta1 = 0.1
beta2 = 0.1
alpha = 1/6
gamma = 1/10

infection = c(3,3,3,3,3,4,6,9,11,21,30,35,42,50,63,74,82,94,101,116,128,134,147,164,171,184,186,191,196,203,207,217,221,233,245)
recovery = c(0,0,0,0,0,0,0,0,0,0,0,0,1,4,4,8,10,11,11,12,15,16,16,17,24,25,40,52,62,70,72,72,80,97,121)

infection = infection - c(3, infection[1:34])
recovery = recovery - c(0, recovery[1:34])

t.y = c()
for (i in 1:T.end){
  if (infection[i] > 0){
    t.y = append(t.y, runif(infection[i], i-1, i))
  }
}
t.y = sort(t.y)
n.y = length(t.y)

t.z = c()
for (i in 1:T.end){
  if (recovery[i] > 0){
    t.z = append(t.z, runif(recovery[i], i-1, i))
  }
}
t.z = sort(t.z)
n.z = length(t.z)

t.x.1 = New_proposal.x(0.3, 0, t.y, t.z, t0, k)
n.x.1 = length(t.x.1)
t.x.2 = New_proposal.x(0.1, n.x.1, t.y, t.z, k, T.end)
n.x.2 = length(t.x.2)
t.x = c(t.x.1, t.x.2)
n.x = n.x.1 + n.x.2

### MCMC ###
M = 5000
sample.nx = c()
sample.beta1 = c()
sample.beta2 = c()
sample.alpha = c()
sample.gamma = c()
m = 10
acc = 0
for (j in 1:M){
  ## Update parameters
  beta1 = rgamma(1, n.x.1+5, SI.int(t.x.1, t.y, t.z, t0, k)/N+50)
  beta2 = rgamma(1, n.x.2+5, SI.int(t.x, t.y, t.z, k, T.end)/N+50)
  alpha = rgamma(1, n.y+100, E_int(t.x, t.y)+600)
  gamma = rgamma(1, n.z+50, I_int(t.y, t.z)+500)
  sample.beta1[j] = beta1
  sample.beta2[j] = beta2
  sample.alpha[j] = alpha
  sample.gamma[j] = gamma



  ## update latent Y
  t.y = c(t0, t.y, T.end)
  for (i in 2:(n.y+1)){
    a = max(t.y[i-1], floor(t.y[i]))
    b = min(t.y[i+1], ceiling(t.y[i]))
    t.new = runif(1, a, b)
    t.y.new = t.y
    t.y.new[i] = t.new

    lp_new = Log.SI.prod(t.x, t.y.new, t.z, a, b) - ifelse(t.new<k, beta1, beta2)/N * SI.int(t.x, t.y.new, t.z, a, b) +
      Log.E.prod(t.x, t.y.new, a, b) - alpha * E.int(t.x, t.y.new, a, b) +
      Log.I.prod(t.y.new, t.z, a, b) - gamma * I.int(t.y.new, t.z, a, b)

    lp_old = Log.SI.prod(t.x, t.y, t.z, a, b) - ifelse(t.new<k, beta1, beta2)/N * SI.int(t.x, t.y, t.z, a, b) +
      Log.E.prod(t.x, t.y, a, b) - alpha * E.int(t.x, t.y, a, b) +
      Log.I.prod(t.y, t.z, a, b) - gamma * I.int(t.y, t.z, a, b)

    lp = lp_new - lp_old
    u = runif(1)
    if(is.na(lp)){
      lp = -Inf
    }
    if(log(u) < lp){
      t.y = t.y.new
    }
  }
  t.y = t.y[2:(n.y+1)]


  ## update latent Z
  t.z = c(t0, t.z, T.end)
  for (i in 2:(n.z+1)){
    a = max(t.z[i-1], floor(t.z[i]))
    b = min(t.z[i+1], ceiling(t.z[i]))
    t.new = runif(1, a, b)
    t.z.new = t.z
    t.z.new[i] = t.new

    lp_new = Log.SI.prod(t.x, t.y, t.z.new, a, b) - ifelse(t.new<k, beta1, beta2)/N * SI.int(t.x, t.y, t.z.new, a, b) +
      Log.I.prod(t.y, t.z.new, a, b) - gamma * I.int(t.y, t.z.new, a, b)

    lp_old = Log.SI.prod(t.x, t.y, t.z, a, b) - ifelse(t.new<k, beta1, beta2)/N * SI.int(t.x, t.y, t.z, a, b) +
      Log.I.prod(t.y, t.z, a, b) - gamma * I.int(t.y, t.z, a, b)

    lp = lp_new - lp_old
    u = runif(1)
    if(is.na(lp)){
      lp = -Inf
    }
    if(log(u) < lp){
      t.z = t.z.new
    }
  }
  t.z = t.z[2:(n.z+1)]


  
  # Update latent X by new MCMC algorithm
  # valid_proposal = FALSE
  # while (valid_proposal == FALSE){
  #   tx1.new = New_proposal.x(beta1, 0, t.y, t.z, t0, k)
  #   tx.new = c(tx1.new, t.x.2)
  #   valid_proposal = Indicator(tx.new, t.y, t.z)
  # }
  # log_p.new = Log.SI.prod(tx.new, t.y, t.z, k, T.end) - beta2/N * SI.int(tx.new, t.y, t.z, k, T.end) +
  #   Log.E_prod(tx.new, t.y) - alpha * E_int(tx.new, t.y)
  # 
  # log_p.old = Log.SI.prod(t.x, t.y, t.z, k, T.end) - beta2/N * SI.int(t.x, t.y, t.z, k, T.end) +
  #   Log.E_prod(t.x, t.y) - alpha * E_int(t.x, t.y)
  # 
  # log_p = log_p.new - log_p.old
  # u = runif(1)
  # if(log(u) < log_p){
  #   t.x.1 = tx1.new
  #   t.x = tx.new
  #   n.x.1 = length(t.x.1)
  #   n.x = n.x.1 + n.x.2
  # }
  # 
  # valid_proposal = FALSE
  # while (valid_proposal == FALSE){
  #   tx2.new = New_proposal.x(beta2, n.x.1, t.y, t.z, k, T.end)
  #   tx.new = c(t.x.1, tx2.new)
  #   valid_proposal = Indicator(tx.new, t.y, t.z)
  # }
  # 
  # log_p.new = Log.E.prod(tx.new, t.y, k, T.end) - alpha * E.int(tx.new, t.y, k, T.end)
  # log_p.old = Log.E.prod(t.x, t.y, k, T.end) - alpha * E.int(t.x, t.y, k, T.end)
  # 
  # log_p = log_p.new - log_p.old
  # u = runif(1)
  # if(log(u) < log_p){
  #   t.x.2 = tx2.new
  #   t.x = tx.new
  #   n.x.2 = length(t.x.2)
  #   n.x = n.x.1 + n.x.2
  # }
  # sample.nx[j] = n.x
  
  
  
  ## Update latent X by current MCMC algorithm
  tx.new = t.x
  log_a = 0
  for (i in 1:m){
    u = runif(1, 0, 1)
    if (u < 1/3){
      # Add one time point
      tx.new = sort(c(runif(1, t0, T.end), tx.new))
      log_a = log_a + log(T.end-t0) - log(length(tx.new))
    }else{
      if (u > 2/3){
        # Move one time point
        tx.new = tx.new[-sample(1:length(tx.new), 1)]
        tx.new = sort(c(runif(1, t0, T.end), tx.new))
      }else{
        # Remove one time point
        tx.new = tx.new[-sample(1:length(tx.new), 1)]
        log_a = log_a + log(length(tx.new)+1) - log(T.end-t0)
      }
    }
  }
  if (Indicator(tx.new, t.y, t.z)){
    tx1.new = tx.new[tx.new<k]
    nx1.new = length(tx1.new)
    tx2.new = tx.new[tx.new>k]
    nx2.new = length(tx2.new)
    nx.new = length(tx.new)
    
    log_p.new = Log.E_prod(tx.new, t.y) - alpha * E_int(tx.new, t.y) +
      nx1.new*log(beta1/N) + nx2.new*log(beta2/N) + Log.SI_prod(tx.new, t.y, t.z) -
      beta1/N * SI.int(tx1.new, t.y, t.z, t0, k) - beta2/N * SI.int(tx.new, t.y, t.z, k, T.end)
    
    log_p.old = Log.E_prod(t.x, t.y) - alpha * E_int(t.x, t.y) +
      n.x.1*log(beta1/N) + n.x.2*log(beta2/N) + Log.SI_prod(t.x, t.y, t.z) -
      beta1/N * SI.int(t.x.1, t.y, t.z, t0, k) - beta2/N * SI.int(t.x, t.y, t.z, k, T.end)
    
    log_a = log_a + log_p.new - log_p.old
    u = runif(1)
    if(log(u) < log_a){
      t.x = tx.new
      t.x.1 = tx1.new
      n.x.1 = nx1.new
      t.x.2 = tx2.new
      n.x.2 = nx2.new
      n.x = nx.new
      acc = acc + 1
    }
  }
  sample.nx[j] = n.x
  print(j)
}
acc/M

Samples.new = list(nx=sample.nx, beta1=sample.beta1, beta2=sample.beta2, alpha=sample.alpha, gamma=sample.gamma)
save(Samples.new, file='SEIR_new.RData')
Samples.current = list(nx=sample.nx, beta1=sample.beta1, beta2=sample.beta2, alpha=sample.alpha, gamma=sample.gamma)
save(Samples.current, file='SEIR_current.RData')


par(mfrow=c(2, 1))
plot(Samples.new$nx, type = 'l',
     main = expression(paste('Trace plot of ', n[x], ' with new MCMC algorithm')),
     ylab = expression(n[x]))
plot(Samples.current$nx, type = 'l',
     main = expression(paste('Trace plot of ', n[x], ' with current MCMC algorithm')),
     ylab = expression(n[x]))

par(mfrow=c(2, 1))
plot(Samples.new$beta1, type = 'l', main = expression(paste('Trace plot of ', beta[1], ' with new MCMC algorithm')),
     ylab = expression(beta[1]))
plot(Samples.current$beta1, type = 'l', main = expression(paste('Trace plot of ', beta[1], ' with current MCMC algorithm')),
     ylab = expression(beta[1]))

par(mfrow=c(2, 1))
plot(Samples.new$beta2, type = 'l', main = expression(paste('Trace plot of ', beta[2], ' with new MCMC algorithm')),
     ylab = expression(beta[2]))
plot(Samples.current$beta2, type = 'l', main = expression(paste('Trace plot of ', beta[2], ' with current MCMC algorithm')),
     ylab = expression(beta[2]))

par(mfrow=c(2, 1))
plot(Samples.new$alpha, type = 'l', main = expression(paste('Trace plot of ', alpha, ' with new MCMC algorithm')),
     ylab = expression(alpha))
plot(Samples.current$alpha, type = 'l', main = expression(paste('Trace plot of ', alpha, ' with current MCMC algorithm')),
     ylab = expression(alpha))

par(mfrow=c(2, 1))
plot(Samples.new$gamma, type = 'l', main = expression(paste('Trace plot of ', gamma, ' with new MCMC algorithm')),
     ylab = expression(gamma))
plot(Samples.current$gamma, type = 'l', main = expression(paste('Trace plot of ', gamma, ' with current MCMC algorithm')),
     ylab = expression(gamma))

mean(Samples.new$nx[200:M])
mean(Samples.new$beta1[200:M])
mean(Samples.new$beta2[200:M])
mean(Samples.new$alpha[200:M])
mean(Samples.new$gamma[200:M])

mean(Samples.current$nx[1000:M])
mean(Samples.current$beta1[1000:M])
mean(Samples.current$beta2[1000:M])
mean(Samples.current$alpha[1000:M])
mean(Samples.current$gamma[1000:M])

### 3. R0 ###
library(R0) 
incidences = c(1,0,0,1,1,1,2,3,2,10,9,5,7,8,13,11,8,12,7,15,12,6,13,17,7,13,2,5,5,7,4,10,4,12,12)

mGT = generation.time ("gamma", c(11, 5))

# before change point
estR0 = estimate.R(incidences, GT=mGT, begin=1, end=16, methods=c("EG", "ML"), 
                   pop.size=285000, nsim=100)
attributes(estR0)
estR0

# after change point
estR0 = estimate.R(incidences, GT=mGT, begin=17, end=35, methods=c("EG", "ML"), 
                   pop.size=285000, nsim=100)
attributes(estR0)
estR0