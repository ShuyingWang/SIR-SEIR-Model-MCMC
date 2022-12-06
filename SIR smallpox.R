### Description of the code file ###
# The first part defines the functions.
# The second part includes the illustration with the smallpox dataset (Section 4.3 in the paper).
# Comment/uncomment the new/current MCMC algorithms to run comparisons.


### 1. Functions ###

# Integrals (for calculating likelihoods)

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


# Products (for calculating likelihoods)

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
X.to.u <- function(t.x, t.y, t0, T.end, beta){
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
  return(U)
}

u.to.X <- function(u, t.y, t0, T.end, beta){
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





### 2. Illustration with smallpox data ###

library(mcmcse)

## Data entry
I0 = 1
N = 120
S0 = 119
t.y = c(0, 13, 7, 2, 3, 0, 0, 1, 4, 5, 3, 2, 0, 2, 0, 5, 3, 1, 4, 0, 1, 0, 1, 1, 2, 0, 1, 2, 3, 0, 5, 5)
t.y = cumsum(t.y)
T.end = 76
n.y = length(t.y)

## MCMC
M = 5000
gamma = 0.09
beta = 0.1
t.x = t.y - rexp(n.y, 1)
t0 = t.x[1]-rexp(1, gamma + beta)
t.x = sort(c(t.x, runif(10, t0, T.end)))
n.x = length(t.x)
sample.nx = c()
sample.t0 = c()
sample.beta = c()
sample.gamma = c()
acc = 0

ptm <- proc.time()
for (j in 1:M){
  
  # update initial infection time
  t0 = t.x[1]-rexp(1, gamma + beta*S0/N + 0.1)
  
  # update beta and gamma
  beta = rgamma(1, n.x + 10, SI_int(t.x, t.y, n.x, n.y, t0, T.end)/N + 100)
  gamma = rgamma(1, n.y + 10, I_int(t.x, t.y, n.x, n.y, t0, T.end) + 100)

  
  ## update latent infection time by new MCMC algorithm
  # u = X.to.u(t.x, t.y, t0, T.end, beta)
  # k = ceiling((n.x + 1) * 0.05)
  # 
  # valid_proposal = FALSE
  # while (valid_proposal == FALSE){
  #   u.new = u
  #   u.new[sample(1:(n.x+1), k, replace = FALSE)] = rexp(k, 1)
  #   tx.new = u.to.X(u.new, t.y, t0, T.end, beta)
  #   nx.new = length(tx.new)
  #   valid_proposal = Indicator(tx.new, t.y)
  # }
  # log_p.new = Log.I_prod(tx.new, t.y) - gamma * I_int(tx.new, t.y, nx.new, n.y, t0, T.end)
  # log_p.old = Log.I_prod(t.x, t.y) - gamma * I_int(t.x, t.y, n.x, n.y, t0, T.end)
  # log_p = log_p.new - log_p.old
  # u = runif(1)
  # if(log(u) < log_p){
  #   t.x = tx.new
  #   n.x = nx.new
  #   acc = acc + 1
  # }
  
  
  # update latent infection time by current MCMC algorithm (Random Walk Metropolis)
  valid_proposal = FALSE
  while (valid_proposal == FALSE){
    u = runif(1)
    if (u < 1/3){
      # Add one time point
      nx.new = n.x + 1
      t.new = runif(1, t0, T.end)
      tx.new = sort(c(t.new, t.x))
      move_type = 1
    }else{
      if (u > 2/3){
        # Move one time point
        nx.new = n.x
        t.new = runif(1, t0, T.end)
        tx.new = t.x[-sample(1:n.x, 1)]
        tx.new = sort(c(t.new, tx.new))
        move_type = 2
      }else{
        # Remove one time point
        nx.new = n.x - 1
        tx.new = t.x[-sample(1:n.x, 1)]
        move_type = 3
      }
    }
    valid_proposal = Indicator(tx.new, t.y)
  }
  log_p.new = nx.new * log(beta/N) + Log.SI_prod(tx.new, t.y) - beta/N * SI_int(tx.new, t.y, nx.new, n.y, t0, T.end) +
    Log.I_prod(tx.new, t.y) - gamma * I_int(tx.new, t.y, nx.new, n.y, t0, T.end)
  log_p.old = n.x * log(beta/N) + Log.SI_prod(t.x, t.y) - beta/N * SI_int(t.x, t.y, n.x, n.y, t0, T.end) +
    Log.I_prod(t.x, t.y) - gamma * I_int(t.x, t.y, n.x, n.y, t0, T.end)
  if (move_type == 1){
    # Add one time point
    log_p = log_p.new + log(T.end-t0) - (log_p.old + log(n.x+1))
  }else{
    if (move_type == 3){
      # Remove one time point
      log_p = log_p.new + log(n.x) - (log_p.old + log(T.end-t0))
    }else{
      # Move one time point
      log_p = log_p.new - log_p.old
    }
  }
  u = runif(1)
  if(log(u) < log_p){
    t.x = tx.new
    n.x = nx.new
    acc = acc + 1
  }
  
  sample.nx[j] = n.x
  sample.t0[j] = t0
  sample.beta[j] = beta
  sample.gamma[j] = gamma
}
proc.time() - ptm
acc/M

acc.new = acc/M #0.204
acc.current = acc/M #0.6398

Samples.new = list(nx=sample.nx, t0=sample.t0, beta=sample.beta, gamma=sample.gamma)
save(Samples.new, file='new.RData')
Samples.current = list(nx=sample.nx, t0=sample.t0, beta=sample.beta, gamma=sample.gamma)
save(Samples.current, file='current.RData')

## Plots
par(mfrow=c(2, 1))
plot(Samples.new$nx, type = 'l',
     main = expression(paste('Trace plot of ', n[x], ' with new MCMC algorithm')),
     ylab = expression(n[x]))
plot(Samples.current$nx, type = 'l',
     main = expression(paste('Trace plot of ', n[x], ' with current MCMC algorithm')),
     ylab = expression(n[x]))

par(mfrow=c(2, 1))
plot(Samples.new$t0, type = 'l', main = expression(paste('Trace plot of ', t[0]^x, ' with new MCMC algorithm')),
     ylab = expression(t[0]^x))
plot(Samples.current$t0, type = 'l', main = expression(paste('Trace plot of ', t[0]^x, ' with current MCMC algorithm')),
     ylab = expression(t[0]^x))

par(mfrow=c(2, 1))
plot(Samples.new$beta, type = 'l', main = expression(paste('Trace plot of ', beta, ' with new MCMC algorithm')),
     ylab = expression(beta))
plot(Samples.current$beta, type = 'l', main = expression(paste('Trace plot of ', beta, ' with current MCMC algorithm')),
     ylab = expression(beta))

par(mfrow=c(2, 1))
plot(Samples.new$gamma, type = 'l', main = expression(paste('Trace plot of ', gamma, ' with new MCMC algorithm')),
     ylab = expression(gamma))
plot(Samples.current$gamma, type = 'l', main = expression(paste('Trace plot of ', gamma, ' with current MCMC algorithm')),
     ylab = expression(gamma))

par(mfrow=c(2, 1))
hist(Samples.new$beta, freq = FALSE, main = expression(paste('Histogram of ', beta, ' with new MCMC algorithm')),
     xlab = expression(beta), xlim=c(0.04, 0.2))
hist(Samples.current$beta, freq = FALSE, main = expression(paste('Histogram of ', beta, ' with current MCMC algorithm')),
     xlab = expression(beta), xlim=c(0.04, 0.2))

par(mfrow=c(2, 1))
hist(Samples.new$gamma, freq = FALSE, main = expression(paste('Histogram of ', gamma, ' with new MCMC algorithm')),
     xlab = expression(gamma), xlim=c(0.03, 0.16))
hist(Samples.current$gamma, freq = FALSE, main = expression(paste('Histogram of ', gamma, ' with current MCMC algorithm')),
     xlab = expression(gamma), xlim=c(0.03, 0.16))

mean(Samples.new$beta) #0.110
mean(Samples.new$gamma) #0.084
var(Samples.new$beta) #0.00032
var(Samples.new$gamma) #0.00028
ess(Samples.new$beta) #1588.848
ess(Samples.new$gamma) #621.0135

mean(Samples.current$beta) #0.101
mean(Samples.current$gamma) #0.077
var(Samples.current$beta) #0.00045
var(Samples.current$gamma) #0.00045
ess(Samples.current$beta) #346.0594
ess(Samples.current$gamma) #174.275
