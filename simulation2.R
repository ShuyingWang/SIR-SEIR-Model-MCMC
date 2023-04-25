### Description of the code file ###
# The first part defines the functions.
# The second part includes the second illustration with simulated data (Section 4.2 in the paper).
# The simulated data can be generated with the seeds set.


### 1. Functions ###

# Data simulation
Sim_SIR <- function(N, I0, beta, gamma, T.end){
  # initial numbers
  S = N - I0
  I = I0
  
  # recording time;
  t.x = c()
  t.y = c()
  t = 0
  
  while (I > 0){
    
    # time to next infection
    if (S > 0){
      T.to.I = rexp(1, (beta/N)*I*S)
    }
    else{
      T.to.I = Inf
    }
    
    # time to next removal
    T.to.R = rexp(1, gamma*I)
    
    # Check time
    t = t + min(T.to.I, T.to.R)
    if (t > T.end){
      break
    }
    
    # Update time and states
    if (T.to.I < T.to.R){
      
      # infection happens
      I = I + 1
      S = S - 1
      t.x = append(t.x, t)
    }else{
      
      # removal happens
      I = I - 1
      t.y = append(t.y, t)
    }
  }
  return(list('x'=t.x, 'y'=t.y))
}

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
  u.new = ifelse(rbinom(n, 1, r), rexp(n, 1), u)
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
  # update beta and gamma
  beta = rgamma(1, n.x+0.1, SI_int(t.x, t.y, n.x, n.y, t0, T.end)/N+0.1)
  gamma = rgamma(1, n.y+0.1, I_int(t.x, t.y, n.x, n.y, t0, T.end)+0.1)
  
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
  return(list("beta"=beta, "gamma"=gamma, "tx"=t.x, "nx"=n.x, "acc"=acc))
}

# Current MCMC algorithm (Random Walk Metropolis)
Current_MH.x <- function(k, t.x, n.x, t.y, n.y, beta, gamma, t0, T.end){
  # update beta and gamma
  beta = rgamma(1, n.x+0.1, SI_int(t.x, t.y, n.x, n.y, t0, T.end)/N+0.1)
  gamma = rgamma(1, n.y+0.1, I_int(t.x, t.y, n.x, n.y, t0, T.end)+0.1)
  
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
  return(list("beta"=beta, "gamma"=gamma, "tx"=t.x, "nx"=n.x, "acc"=acc))
}





### 2. Illustration using simulated data ###

library(mcmcse)
set.seed(3)

# Data simulation
I0 = 100
N = 5000
S0 = N - I0
beta = 0.3
gamma = 0.1
t0 = 0
T.end = 30

t = Sim_SIR(N, I0, beta, gamma, T.end)
t.x = t$x
t.y = t$y
nx.true = length(t.x)
n.y = length(t.y)
nx.true
n.y





# New MCMC
set.seed(3)
B = 20000
beta = 0.3
gamma = 0.1
t.x = New_proposal.x(t.y, n.y, beta, t0, T.end)
n.x = length(t.x)
sample_nx.new = rep(NA, B)
sample_beta.new = rep(NA, B)
sample_gamma.new = rep(NA, B)
acc.new = 0
r = 0.1
for (i in 1:B){
  Sample = New_MH_adjust.x(r, t.x, n.x, t.y, n.y, beta, gamma, t0, T.end)
  t.x = Sample$tx
  n.x = Sample$nx
  beta = Sample$beta
  gamma = Sample$gamma
  acc.new = acc.new + Sample$acc
  sample_nx.new[i] = n.x
  sample_beta.new[i] = beta
  sample_gamma.new[i] = gamma
  print(i)
}
acc.new/B # 0.25



# Random walk MCMC
set.seed(3)
B = 20000
beta = 0.3
gamma = 0.1
t.x = New_proposal.x(t.y, n.y, beta, t0, T.end)
n.x = length(t.x)
sample_nx.current = rep(NA, B)
sample_beta.current = rep(NA, B)
sample_gamma.current = rep(NA, B)
acc.current = 0
k = 20
for (i in 1:B){
  Sample = Current_MH.x(k, t.x, n.x, t.y, n.y, beta, gamma, t0, T.end)
  t.x = Sample$tx
  n.x = Sample$nx
  beta = Sample$beta
  gamma = Sample$gamma
  acc.current = acc.current + Sample$acc
  sample_nx.current[i] = n.x
  sample_beta.current[i] = beta
  sample_gamma.current[i] = gamma
  print(i)
}
acc.current/B # 0.22



Samples.new = list(nx=sample_nx.new, beta=sample_beta.new, gamma=sample_gamma.new, t.x)
save(Samples.new, file='New.RData')
Samples.current = list(nx=sample_nx.current, beta=sample_beta.current, gamma=sample_gamma.current, t.x)
save(Samples.current, file='Current.RData')


# Plots 
par(mfrow=c(2, 1))
plot(Samples.new$nx, type = 'l',
     main = expression(paste('Trace plot of ', n[x], ' with new MCMC algorithm')),
     ylab = expression(n[x]))
plot(Samples.current$ny, type = 'l',
     main = expression(paste('Trace plot of ', n[x], ' with current MCMC algorithm')),
     ylab = expression(n[x]))

par(mfrow=c(2, 1))
plot(Samples.new$beta, type = 'l',
     main = expression(paste('Trace plot of ', beta, ' with new MCMC algorithm')),
     ylab = expression(beta))
plot(Samples.current$beta, type = 'l',
     main = expression(paste('Trace plot of ', beta, ' with current MCMC algorithm')),
     ylab = expression(beta))

par(mfrow=c(2, 1))
plot(Samples.new$gamma, type = 'l',
     main = expression(paste('Trace plot of ', gamma, ' with new MCMC algorithm')),
     ylab = expression(gamma))
plot(Samples.current$gamma, type = 'l',
     main = expression(paste('Trace plot of ', gamma, ' with current MCMC algorithm')),
     ylab = expression(gamma))

par(mfrow=c(2, 1))
hist(Samples.new$nx[1000:B], freq=FALSE, breaks=10, xlim=c(3500, 4600),
     main = expression(paste('Histogram of ', n[x], ' with new MCMC')),
     xlab = expression(n[x]), ylim=c(0, 0.004))
hist(Samples.current$nx[1000:B], freq=FALSE, breaks=10, xlim=c(3500, 4600),
     main = expression(paste('Histogram of ', n[x], ' with current MCMC')),
     xlab = expression(n[x]), ylim=c(0, 0.004))

par(mfrow=c(2, 1))
hist(Samples.new$beta[1000:B], freq=FALSE,
     main = expression(paste('Histogram of ', beta, ' with new MCMC')),
     xlab = expression(beta), xlim=c(0.28, 0.34), ylim=c(0, 55))
hist(Samples.current$beta[1000:B], freq=FALSE,
     main = expression(paste('Histogram of ', beta, ' with current MCMC')),
     xlab = expression(beta), xlim=c(0.28, 0.34), ylim=c(0, 55))

par(mfrow=c(2, 1))
hist(Samples.new$gamma[1000:B], freq=FALSE,
     main = expression(paste('Histogram of ', gamma, ' with new MCMC')),
     xlab = expression(gamma), xlim=c(0.05, 0.15), ylim=c(0, 28))
hist(Samples.current$gamma[1000:B], freq=FALSE,
     main = expression(paste('Histogram of ', gamma, ' with current MCMC')),
     xlab = expression(gamma), xlim=c(0.05, 0.15), ylim=c(0, 28))

# Posterior means
mean(Samples.new$nx[1000:B]) # 4208
mean(Samples.current$ny[1000:B]) # 3788

mean(Samples.new$beta[1000:B]) # 0.307
mean(Samples.current$beta[1000:B]) # 0.309

mean(Samples.new$gamma[1000:B]) # 0.091
mean(Samples.current$gamma[1000:B]) # 0.122


