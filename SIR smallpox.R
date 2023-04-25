### Description of the code file ###
# The first part defines the functions.
# The second part includes the illustration with the smallpox dataset (Section 4.3 in the paper).

library(mcmcse)
library(ReIns)

### 1. Functions ###

# Likelihood calculation
SI_int <- function(t.x, t.y, n.x, n.y, t0, T.end){
  t = sort(c(t.x, t.y, t0, T.end))
  n = n.x + n.y + 1
  S = S0
  I = I0
  int = 0
  
  for (i in 1:n){
    int = int + (t[i+1] - t[i]) * S * I
    if (t[i+1] %in% t.x){
      S = S - 1
      I = I + 1
    }else{
      I = I - 1
    }
  }
  return(int)
}
I_int <- function(t.x, t.y, n.x, n.y, t0, T.end){
  t = sort(c(t.x, t.y, t0, T.end))
  n = n.x + n.y + 1
  I = I0
  int = 0
  
  for (i in 1:n){
    int = int + (t[i+1] - t[i]) * I
    if (t[i+1] %in% t.x){
      I = I + 1
    }else{
      I = I - 1
    }
  }
  return(int)
}
Log.SI_prod <- function(t.x, t.y){
  times = sort(c(t.x, t.y))
  S = S0
  I = I0
  lp = 0
  
  for (t in times){
    if (t %in% t.x){
      lp = lp + log(S) + log(I)
      S = S - 1
      I = I + 1
    }else{
      I = I - 1
    }
  }
  return(lp)
}
Log.I_prod <- function(t.x, t.y){
  times = sort(c(t.x, t.y))
  I = I0
  lp = 0
  
  for (t in times){
    if (t %in% t.y){
      lp = lp + log(I)
      I = I - 1
    }else{
      I = I + 1
    }
  }
  return(lp)
}
Indicator <- function(t.x, t.y){
  times = sort(c(t.x, t.y))
  S = S0
  I = I0
  ind = 1
  
  for (t in times){
    if (t %in% t.x){
      if (S<=0 | I<=0){
        ind = 0
      }
      S = S - 1
      I = I + 1
    }else{
      if (I <= 0){
        ind = 0
      }
      I = I - 1
    }
  }
  return(ind)
}


# Conversions
X_to_u <- function(t.x, t.y, t0, T.end, beta){
  t.x = c(t.x, T.end)
  t = sort(c(t.x, t.y, t0))
  S = S0
  I = I0
  U = c()
  u = 0
  
  for (i in 1:(length(t) - 1)){
    u = u + (t[i+1] - t[i]) * S * I
    if (t[i+1] %in% t.x){
      S = S - 1
      I = I + 1
      U = c(U, u)
      u = 0
    }else{
      I = I - 1
    }
  }
  U = U * beta / N
  n = length(U)
  U[n] = U[n] + rexp(1, 1)
  return(U)
}
u_to_u.new <- function(u, r){
  n = length(u)
  u.new = ifelse(rbinom(n,1,r), rexp(n,1), u)
  return(u.new)
}
u_to_X <- function(u, t.y, t0, T.end, beta){
  t.y = c(t.y, T.end)
  n = length(u)
  S = S0
  I = I0
  t.x = c()
  t = t0
  j = 1
  i = 1
  
  while (t < T.end & S > 0 & I > 0){
    
    if ((t + u[i] / (beta*S*I/N)) < t.y[j]){
      t = t + u[i] / (beta*S*I/N)
      t.x = c(t.x, t)
      S = S - 1
      I = I + 1
      i = i + 1
      if (i > n){
        u[i] = rexp(1, 1)
      }
    }else{
      u[i] = u[i] - (t.y[j] - t) * (beta*S*I/N)
      t = t.y[j]
      I = I - 1
      j = j + 1
    }
  }
  return(t.x)
}


# New MCMC algorithm 
New_proposal.x <- function(t.y, n.y, beta, t0, T.end){
  ty = c(t.y, T.end)
  S = S0
  I = I0
  t = t0
  tx.new = c()
  i = 1
  while (t < T.end & S > 0 & I > 0){
    t_new = t + rexp(1, beta/N * S * I)
    if (I==1 & i<n.y & t_new>ty[i]){
      S = S0
      I = I0
      t = t0
      tx.new = c()
      i = 1
    }else if (t_new < ty[i]){
      t = t_new
      tx.new = c(tx.new, t)
      I = I + 1
      S = S - 1
    }else{
      t = ty[i]
      I = I - 1
      i = i + 1
    }
  }
  return(tx.new)
}
New_MH_adjust.x <- function(r, t.x, n.x, t.y, n.y, beta, gamma, t0, T.end){
  
  # update initial infection time
  t0 = t.x[1]-rexp(1, gamma + beta*S0/N + 0.1)
  
  # update beta and gamma
  beta = rgamma(1, n.x+10, SI_int(t.x, t.y, n.x, n.y, t0, T.end)/N+100)
  gamma = rgamma(1, n.y+10, I_int(t.x, t.y, n.x, n.y, t0, T.end)+100)
  
  # update X
  u = X_to_u(t.x, t.y, t0, T.end, beta)
  u.new = u_to_u.new(u, r)
  tx.new = u_to_X(u.new, t.y, t0, T.end, beta)
  nx.new = length(tx.new)
  acc = 0
  if (Indicator(tx.new, t.y)){
    log_p = Log.I_prod(t.x, t.y) - gamma * I_int(t.x, t.y, n.x, n.y, t0, T.end)
    log_p.new = Log.I_prod(tx.new, t.y) - gamma * I_int(tx.new, t.y, nx.new, n.y, t0, T.end)
    log_a = log_p.new-log_p
    acc = 0
    if (log(runif(1)) < log_a){
      t.x = tx.new
      n.x = nx.new
      acc = 1
    }
  }
  return(list("t0"=t0, "beta"=beta, "gamma"=gamma, "tx"=t.x, "nx"=n.x, "acc"=acc))
}

# Current MCMC algorithm (Random Walk Metropolis)
Current_MH.x <- function(k, t.x, n.x, t.y, n.y, beta, gamma, t0, T.end){
  
  # update initial infection time
  t0 = t.x[1]-rexp(1, gamma + beta*S0/N + 0.1)
  
  # update beta and gamma
  beta = rgamma(1, n.x+10, SI_int(t.x, t.y, n.x, n.y, t0, T.end)/N+100)
  gamma = rgamma(1, n.y+10, I_int(t.x, t.y, n.x, n.y, t0, T.end)+100)
  
  # update X
  tx.new = t.x
  log_a = 0
  for (i in 1:k){
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
  nx.new = length(tx.new)
  acc = 0
  if (Indicator(tx.new, t.y)){
    log_p = n.x * log(beta/N) + Log.SI_prod(t.x, t.y) - beta/N * SI_int(t.x, t.y, n.x, n.y, t0, T.end) +
      Log.I_prod(t.x, t.y) - gamma * I_int(t.x, t.y, n.x, n.y, t0, T.end)
    log_p.new = nx.new * log(beta/N) + Log.SI_prod(tx.new, t.y) - beta/N * SI_int(tx.new, t.y, nx.new, n.y, t0, T.end) +
      Log.I_prod(tx.new, t.y) - gamma * I_int(tx.new, t.y, nx.new, n.y, t0, T.end)
    log_a = log_a + log_p.new - log_p
    if (log(runif(1)) < log_a){
      t.x = tx.new
      n.x = nx.new
      acc = 1
    }
  }
  return(list("t0"=t0, "beta"=beta, "gamma"=gamma, "tx"=t.x, "nx"=n.x, "acc"=acc))
}

# Simulation-based Bayesian inference for epidemic models
Simulation_MH.x <- function(t.x, n.x, t.y, n.y, beta, gamma, t0, T.end){
  
  # update initial infection time
  t0 = t.x[1]-rexp(1, gamma + beta*S0/N + 0.1)
  
  # update beta and gamma
  beta = rgamma(1, n.x+10, SI_int(t.x, t.y, n.x, n.y, t0, T.end)/N+100)
  gamma = rgamma(1, n.y+10, I_int(t.x, t.y, n.x, n.y, t0, T.end)+100)
  
  # update X
  S = S0
  I = I0
  t = t0
  tx.new = c()
  lp.new = 0
  for (i in 0:(T.end-1)){
    ny = sum(t.y==i)
    nx = ny+1 - I
    if (nx > 0){
      for (j in 1:nx){
        lp.new = lp.new + log(1 - exp(-beta*S*I/N * (i-t)))
        t = t + rtexp(1, beta/N * S * I, i-t)
        tx.new = c(tx.new, t)
        I = I + 1
        S = S - 1
      }
    }
    if (S > 0){
      t_new = t + rexp(1, beta/N * S * I)
      while (t_new < i){
        t = t_new
        tx.new = c(tx.new, t)
        I = I + 1
        S = S - 1
        if (S == 0){
          break
        }
        t_new = t + rexp(1, beta/N * S * I) 
      }
    }
    t = i
    I = I - ny
  }
  if (S > 0){
    t_new = t + rexp(1, beta/N * S * I)
    while (t_new < T.end){
      t = t_new
      tx.new = c(tx.new, t)
      I = I + 1
      S = S - 1
      if (S == 0){
        break
      }
      t_new = t + rexp(1, beta/N * S * I)
    }
  }
  nx.new = length(tx.new)
  log_p.new = lp.new + Log.I_prod(tx.new, t.y) - gamma * I_int(tx.new, t.y, nx.new, n.y, t0, T.end)
  
  S = S0
  I = I0
  t = t0
  k = 1
  tx = c(t.x, T.end)
  lp = 0
  for (i in 0:(T.end-1)){
    ny = sum(t.y==i)
    nx = ny+1 - I
    if (nx > 0){
      for (j in 1:nx){
        lp = lp + log(1 - exp(-beta*S*I/N * (i-t)))
        t = tx[k]
        k = k + 1
        I = I + 1
        S = S - 1
      }
    }
    while (tx[k] < i){
      t = tx[k]
      k = k + 1
      I = I + 1
      S = S - 1
    }
    t = i
    I = I - ny
  }
  log_p = lp + Log.I_prod(t.x, t.y) - gamma * I_int(t.x, t.y, n.x, n.y, t0, T.end)
  acc = 0
  if (log(runif(1)) < (log_p.new-log_p)){
    t.x = tx.new
    n.x = nx.new
    acc = 1
  }
  return(list("t0"=t0, "beta"=beta, "gamma"=gamma, "tx"=t.x, "nx"=n.x, "acc"=acc))
}




### 2. Illustration with smallpox data ###

# Data entry
I0 = 1
N = 120
S0 = 119
t.y = c(0, 13, 7, 2, 3, 0, 0, 1, 4, 5, 3, 2, 0, 2, 0, 5, 3, 1, 4, 0, 1, 0, 1, 1, 2, 0, 1, 2, 3, 0, 5, 5)
t.y = cumsum(t.y)
T.end = 76
n.y = length(t.y)



# New MCMC
set.seed(3)
B = 10000
gamma = 0.1
beta = 0.1
t.x = t.y - rexp(n.y, 1)
t0 = t.x[1]-rexp(1, gamma + beta)
t.x = sort(c(t.x, runif(10, t0, T.end)))
n.x = length(t.x)
sample_nx.new = rep(NA, B)
sample_t0.new = rep(NA, B)
sample_beta.new = rep(NA, B)
sample_gamma.new = rep(NA, B)
acc.new = 0
r = 0.05
for (i in 1:B){
  Sample = New_MH_adjust.x(r, t.x, n.x, t.y, n.y, beta, gamma, t0, T.end)
  t.x = Sample$tx
  n.x = Sample$nx
  t0 = Sample$t0
  beta = Sample$beta
  gamma = Sample$gamma
  acc.new = acc.new + Sample$acc
  sample_nx.new[i] = n.x
  sample_t0.new[i] = t0
  sample_beta.new[i] = beta
  sample_gamma.new[i] = gamma
  print(i)
}
acc.new/B #0.2784



# Random walk MCMC
set.seed(3)
B = 10000
gamma = 0.1
beta = 0.1
t.x = t.y - rexp(n.y, 1)
t0 = t.x[1]-rexp(1, gamma + beta)
t.x = sort(c(t.x, runif(10, t0, T.end)))
n.x = length(t.x)
sample_nx.current = rep(NA, B)
sample_t0.current = rep(NA, B)
sample_beta.current = rep(NA, B)
sample_gamma.current = rep(NA, B)
acc.current = 0
k = 5
for (i in 1:B){
  Sample = Current_MH.x(k, t.x, n.x, t.y, n.y, beta, gamma, t0, T.end)
  t.x = Sample$tx
  n.x = Sample$nx
  t0 = Sample$t0
  beta = Sample$beta
  gamma = Sample$gamma
  acc.current = acc.current + Sample$acc
  sample_nx.current[i] = n.x
  sample_t0.current[i] = t0
  sample_beta.current[i] = beta
  sample_gamma.current[i] = gamma
  print(i)
}
acc.current/B # 0.2556



# Simulation-based MCMC
set.seed(3)
B = 10000
gamma = 0.1
beta = 0.1
t.x = t.y - rexp(n.y, 1)
t0 = t.x[1]-rexp(1, gamma + beta)
t.x = sort(c(t.x, runif(10, t0, T.end)))
n.x = length(t.x)
sample_nx.simulate = rep(NA, B)
sample_t0.simulate = rep(NA, B)
sample_beta.simulate = rep(NA, B)
sample_gamma.simulate = rep(NA, B)
acc.simulate = 0
for (i in 1:B){
  Sample = Simulation_MH.x(t.x, n.x, t.y, n.y, beta, gamma, t0, T.end)
  t.x = Sample$tx
  n.x = Sample$nx
  t0 = Sample$t0
  beta = Sample$beta
  gamma = Sample$gamma
  acc.simulate = acc.simulate + Sample$acc
  sample_nx.simulate[i] = n.x
  sample_t0.simulate[i] = t0
  sample_beta.simulate[i] = beta
  sample_gamma.simulate[i] = gamma
  print(i)
}
acc.simulate/B # 0.004




Samples.new = list(nx=sample_nx.new, t0=sample_t0.new, beta=sample_beta.new, gamma=sample_gamma.new)
save(Samples.new, file='new.RData')
Samples.current = list(nx=sample_nx.current, t0=sample_t0.current, beta=sample_beta.current, gamma=sample_gamma.current)
save(Samples.current, file='current.RData')
Samples.simulate = list(nx=sample_nx.simulate, t0=sample_t0.simulate, beta=sample_beta.simulate, gamma=sample_gamma.simulate)
save(Samples.simulate, file='simulate.RData')

## Plots
par(mfrow=c(3, 1))
plot(Samples.new$nx, type = 'l',
     main = expression(paste('Trace plot of ', n[x], ' with new MCMC algorithm')),
     ylab = expression(n[x]))
plot(Samples.current$nx, type = 'l',
     main = expression(paste('Trace plot of ', n[x], ' with current MCMC algorithm')),
     ylab = expression(n[x]))
plot(Samples.simulate$nx, type = 'l',
     main = expression(paste('Trace plot of ', n[x], ' with simulation-based MCMC algorithm')),
     ylab = expression(n[x]))

par(mfrow=c(3, 1))
plot(Samples.new$t0, type = 'l', main = expression(paste('Trace plot of ', t[0]^x, ' with new MCMC algorithm')),
     ylab = expression(t[0]^x))
plot(Samples.current$t0, type = 'l', main = expression(paste('Trace plot of ', t[0]^x, ' with current MCMC algorithm')),
     ylab = expression(t[0]^x))
plot(Samples.simulate$t0, type = 'l', main = expression(paste('Trace plot of ', t[0]^x, ' with simulation-based MCMC algorithm')),
     ylab = expression(t[0]^x))

par(mfrow=c(3, 1))
plot(Samples.new$beta, type = 'l', main = expression(paste('Trace plot of ', beta, ' with new MCMC algorithm')),
     ylab = expression(beta))
plot(Samples.current$beta, type = 'l', main = expression(paste('Trace plot of ', beta, ' with current MCMC algorithm')),
     ylab = expression(beta))
plot(Samples.simulate$beta, type = 'l', main = expression(paste('Trace plot of ', beta, ' with simulation-based MCMC algorithm')),
     ylab = expression(beta))

par(mfrow=c(3, 1))
plot(Samples.new$gamma, type = 'l', main = expression(paste('Trace plot of ', gamma, ' with new MCMC algorithm')),
     ylab = expression(gamma))
plot(Samples.current$gamma, type = 'l', main = expression(paste('Trace plot of ', gamma, ' with current MCMC algorithm')),
     ylab = expression(gamma))
plot(Samples.simulate$gamma, type = 'l', main = expression(paste('Trace plot of ', gamma, ' with simulation-based MCMC algorithm')),
     ylab = expression(gamma))

par(mfrow=c(3, 1))
hist(Samples.new$beta, freq = FALSE, main = expression(paste('Histogram of ', beta, ' with new MCMC algorithm')),
     xlab = expression(beta), xlim=c(0.04, 0.2))
hist(Samples.current$beta, freq = FALSE, main = expression(paste('Histogram of ', beta, ' with current MCMC algorithm')),
     xlab = expression(beta), xlim=c(0.04, 0.2))
hist(Samples.simulate$beta, freq = FALSE, main = expression(paste('Histogram of ', beta, ' with simulation-based MCMC algorithm')),
     xlab = expression(beta), xlim=c(0.04, 0.2))

par(mfrow=c(3, 1))
hist(Samples.new$gamma, freq = FALSE, main = expression(paste('Histogram of ', gamma, ' with new MCMC algorithm')),
     xlab = expression(gamma), xlim=c(0.03, 0.16))
hist(Samples.current$gamma, freq = FALSE, main = expression(paste('Histogram of ', gamma, ' with current MCMC algorithm')),
     xlab = expression(gamma), xlim=c(0.03, 0.16))
hist(Samples.simulate$gamma, freq = FALSE, main = expression(paste('Histogram of ', gamma, ' with simulation-based MCMC algorithm')),
     xlab = expression(gamma), xlim=c(0.03, 0.16))

mean(Samples.new$beta) #0.110
mean(Samples.new$gamma) #0.086
var(Samples.new$beta) #0.00035
var(Samples.new$gamma) #0.00030
ess(Samples.new$beta) #1566.731
ess(Samples.new$gamma) #784.1569

mean(Samples.current$beta) #0.105
mean(Samples.current$gamma) #0.081
var(Samples.current$beta) #0.00046
var(Samples.current$gamma) #0.00042
ess(Samples.current$beta) #588.9186
ess(Samples.current$gamma) #390.2475

mean(Samples.simulate$beta) #
mean(Samples.simulate$gamma) #
var(Samples.simulate$beta) #
var(Samples.simulate$gamma) #
ess(Samples.simulate$beta) #
ess(Samples.simulate$gamma) #
