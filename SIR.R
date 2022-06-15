### Description of the code file ###
# The first part defines the functions.
# The second part includes the first three illustrations presented in the paper.
# The simulated data used in illustration 1 and illustration 2 can be generated with the seeds set.
# In each illustration, comment the new MCMC algorithm and uncomment the current MCMC algorithm to run comparisons.


### Define functions ###

## 1. Data simulation

Sim_SIR <- function(N, I0, beta, gamma){
  # initial numbers
  S = N - I0
  I = I0
  
  # recording time;
  times = c()
  
  # type of event (1 = infection, 2 = removal)
  type = c()
  
  while (I > 0 & length(type) < 10000){
    
    # time to next infection
    if (S > 0){
      T.to.I = rexp(1, (beta/N)*I*S)
    }
    else{
      T.to.I = Inf
    }
    
    # time to next removal
    T.to.R = rexp(1, gamma*I);
    
    # check which of the two events happens first
    T.to.next.event = min(T.to.I, T.to.R)
    type.next.event = which.min(c(T.to.I, T.to.R))
    
    if (type.next.event == 1){
      
      # infection happens first
      I = I + 1
      S = S - 1
      type = append(type, 1)
      times = append(times, T.to.I)
    }
    else{
      
      # removal happens first
      I = I - 1
      type = append(type, 2)
      times = append(times, T.to.R)
    }
  }
  res = list("t"=times, "type"=type)
  return(res)
}



## 2. Metropolis Hasting proposals

Generate_x <- function(t.y, t0, T.end){
  t.y = c(t.y, T.end)
  S = S0
  I = I0
  t = t0
  t.x = c()
  i = 1
  
  while (t < T.end & S > 0 & I > 0){
    t_new = t + rexp(1, beta/N * S * I)
    if (t_new < t.y[i]){
      t = t_new
      t.x = c(t.x, t)
      I = I + 1
      S = S - 1
    }else{
      t = t.y[i]
      I = I - 1
      i = i + 1
    }
  }
  return(t.x)
}


Generate_y <- function(t.x, t0, T.end){
  t.x = c(t.x, T.end)
  I = I0
  t = t0
  t.y = c()
  i = 1
  
  while (t < T.end & I > 0){
    t_new = t + rexp(1, gamma * I)
    if (t_new < t.x[i]){
      t = t_new
      t.y = c(t.y, t)
      I = I - 1
    }else{
      t = t.x[i]
      I = I + 1
      i = i + 1
    }
  }
  return(t.y)
}



## 3. Integrals (for calculating likelihoods)

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



## 4. Products (for calculating likelihoods)

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



## 4. Conversions
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





### Illustrations ###

## 1. First simple illustration using simulated data

library(mcmcse)
set.seed(6)

# Data simulation
I0 = 100
N = 1000000
S0 = N - I0
beta = 0.2
gamma = 0.2

res = Sim_SIR(N, I0, beta, gamma)
time_points = cumsum(res$t)
type = res$type

t.x = time_points[type==1]
t.y = time_points[type==2]

T.end = 10

t.x = t.x[t.x < T.end]
nx.true = length(t.x)
t.y = t.y[t.y < T.end]
n.y = length(t.y)

# MCMC
M = 1000
n.x = sample(100:300, 1)
t.x = runif(n.x, 0, T.end)
sample.nx = c()
acc = 0

ptm <- proc.time()
for (j in 1:M){
  
  sample.nx[j] = n.x
  
  ## New MCMC algorithm
  t.x.new = Generate_x(t.y, 0, T.end)
  n.x.new = length(t.x.new)
  lp_new = Log.I_prod(t.x.new, t.y) - gamma * I_int(t.x.new, t.y, n.x.new, n.y, 0, T.end)
  lp_old = Log.I_prod(t.x, t.y) - gamma * I_int(t.x, t.y, n.x, n.y, 0, T.end)
  lp = lp_new - lp_old
  u = runif(1)
  if(is.na(lp)){
    lp = -Inf
  }
  if(log(u) < lp){
    t.x = t.x.new
    n.x = n.x.new
    acc = acc + 1
  }
  
  
  ## Current MCMC algorithm
  # u = runif(1, 0, 1)
  # if (u < 1/3){
  #   # Add one time point
  #   n.x.new = n.x + 1
  #   t.new = runif(1, 0, T.end)
  #   t.x.new = sort(c(t.new, t.x))
  #   move_type = 1
  # }else{
  #   if (u > 2/3){
  #     # Move one time point
  #     n.x.new = n.x
  #     t.new = runif(1, 0, T.end)
  #     t.x.new = t.x[-sample(1:n.x, 1)]
  #     t.x.new = sort(c(t.new, t.x.new))
  #     move_type = 2
  #   }else{
  #     # Remove one time point
  #     n.x.new = n.x - 1
  #     t.x.new = t.x[-sample(1:n.x, 1)]
  #     move_type = 3
  #   }
  # }
  # lp_new = n.x.new * log(beta/N) + Log.SI_prod(t.x.new, t.y) - beta/N * SI_int(t.x.new, t.y, n.x.new, n.y, 0, T.end) +
  #   Log.I_prod(t.x.new, t.y) - gamma * I_int(t.x.new, t.y, n.x.new, n.y, 0, T.end)
  # 
  # lp_old = n.x * log(beta/N) + Log.SI_prod(t.x, t.y) - beta/N * SI_int(t.x, t.y, n.x, n.y, 0, T.end) +
  #   Log.I_prod(t.x, t.y) - gamma * I_int(t.x, t.y, n.x, n.y, 0, T.end)
  # 
  # if (move_type == 1){
  #   # Add one time point
  #   lp = lp_new + log(T.end-0) - (lp_old + log(n.x+1))
  # }else{
  #   if (move_type == 3){
  #     # Remove one time point
  #     lp = lp_new + log(n.x) - (lp_old + log(T.end-0))
  #   }else{
  #     # Move one time point
  #     lp = lp_new - lp_old
  #   }
  # }
  # u = runif(1)
  # if(is.na(lp)){
  #   lp = -Inf
  # }
  # if(log(u) < lp){
  #   t.x = t.x.new
  #   n.x = n.x.new
  #   acc = acc + 1
  # }
}
proc.time() - ptm
ess(sample.nx)
acc/M


## Plots for new MCMC algorithm
plot(sample.nx, type = 'l',
     main = expression(paste('Trace plot of ', n[x], ' with new MCMC algorithm')),
     ylab = expression(n[x]))
acf(sample.nx, main = expression(paste('Auto-correlation of ', n[x], ' with new MCMC algorithm')))


## Plots for current MCMC algorithm
# plot(sample.nx, type = 'l',
#       main = expression(paste('Trace plot of ', n[x], ' with current MCMC algorithm')),
#       ylab = expression(n[x]))
# acf(sample.nx, main = expression(paste('Auto-correlation of ', n[x], ' with current MCMC algorithm')))





## 2. Second illustration using simulated data

library(mcmcse)
set.seed(2)

# Data simulation
I0 = 100
N = 1000000
S0 = N - I0
beta = 0.25
gamma = 0.15

res = Sim_SIR(N, I0, beta, gamma)
time_points = cumsum(res$t)
type = res$type

t.x = time_points[type==1]
t.y = time_points[type==2]

T.end = 10

t.x = t.x[t.x < T.end]
n.x = length(t.x)
t.y = t.y[t.y < T.end]
ny.true = length(t.y)


# MCMC 
M = 3000
n.y = sample(200:450, 1)
t.y = runif(n.y, 0, T.end)
sample.ny = c()
sample.beta = c()
sample.gamma = c()
acc = 0

ptm <- proc.time()
for (j in 1:M){
  
  sample.ny[j] = n.y
  
  # update beta and gamma
  beta = rgamma(1, n.x+0.1, SI_int(t.x, t.y, n.x, n.y, 0, T.end)/N+0.1)
  gamma = rgamma(1, n.y+0.1, I_int(t.x, t.y, n.x, n.y, 0, T.end)+0.1)
  sample.beta[j] = beta
  sample.gamma[j] = gamma
  

  ## update latent recovery times by new MCMC algorithm
  t.y.new = Generate_y(t.x, 0, T.end)
  n.y.new = length(t.y.new)
  lp_new = Log.SI_prod(t.x, t.y.new) - beta/N * SI_int(t.x, t.y.new, n.x, n.y.new, 0, T.end)
  lp_old = Log.SI_prod(t.x, t.y) - beta/N * SI_int(t.x, t.y, n.x, n.y, 0, T.end)
  lp = lp_new - lp_old
  u = runif(1)
  if(is.na(lp)){
    lp = -Inf
  }
  if(log(u) < lp){
    t.y = t.y.new
    n.y = n.y.new
    acc = acc + 1
  }
  
  
  ## update latent recovery times by current MCMC algorithm
  # u = runif(1, 0, 1)
  # if (u < 1/3){
  #   # Add one time point
  #   n.y.new = n.y + 1
  #   t.new = runif(1, 0, T.end)
  #   t.y.new = sort(c(t.new, t.y))
  #   move_type = 1
  # }else{
  #   if (u > 2/3){
  #     # Move one time point
  #     n.y.new = n.y
  #     t.new = runif(1, 0, T.end)
  #     t.y.new = t.y[-sample(1:n.y, 1)]
  #     t.y.new = sort(c(t.new, t.y.new))
  #     move_type = 2
  #   }else{
  #     # Remove one time point
  #     n.y.new = n.y - 1
  #     t.y.new = t.y[-sample(1:n.y, 1)]
  #     move_type = 3
  #   }
  # }
  # lp_new = Log.SI_prod(t.x, t.y.new) - beta/N * SI_int(t.x, t.y.new, n.x, n.y.new, 0, T.end) +
  #   n.y.new * log(gamma) + Log.I_prod(t.x, t.y.new) - gamma * I_int(t.x, t.y.new, n.x, n.y.new, 0, T.end)
  # lp_old = Log.SI_prod(t.x, t.y) - beta/N * SI_int(t.x, t.y, n.x, n.y, 0, T.end) +
  #   n.y * log(gamma) + Log.I_prod(t.x, t.y) - gamma * I_int(t.x, t.y, n.x, n.y, 0, T.end)
  # if (move_type == 1){
  #   # Add one time point
  #   lp = lp_new + log(T.end-0) - (lp_old + log(n.y+1))
  # }else{
  #   if (move_type == 3){
  #     # Remove one time point
  #     lp = lp_new + log(n.y) - (lp_old + log(T.end-0))
  #   }else{
  #     # Move one time point
  #     lp = lp_new - lp_old
  #   }
  # }
  # u = runif(1)
  # if(is.na(lp)){
  #   lp = -Inf
  # }
  # if(log(u) < lp){
  #   t.y = t.y.new
  #   n.y = n.y.new
  #   acc = acc + 1
  # }
}
proc.time() - ptm

## Plots for new MCMC algorithm
plot(sample.beta, type = 'l', main = expression(paste('Trace plot of ', beta, ' with new MCMC algorithm')),
     ylab = expression(beta))
plot(sample.gamma, type = 'l', main = expression(paste('Trace plot of ', gamma, ' with new MCMC algorithm')),
     ylab = expression(gamma))
mean(sample.beta[300:M])
mean(sample.gamma[300:M])


## Plots for current MCMC algorithm
# plot(sample.beta, type = 'l', main = expression(paste('Trace plot of ', beta, ' with current MCMC algorithm')),
#      ylab = expression(beta))
# plot(sample.gamma, type = 'l', main = expression(paste('Trace plot of ', gamma, ' with current MCMC algorithm')),
#      ylab = expression(gamma))
# mean(sample.beta[8000:M])
# mean(sample.gamma[8000:M])





## 3. Third illustration using real data

library(mcmcse)

# Data entry
I0 = 1
N = 120
S0 = 119
t.y = c(0, 13, 7, 2, 3, 0, 0, 1, 4, 5, 3, 2, 0, 2, 0, 5, 3, 1, 4, 0, 1, 0, 1, 1, 2, 0, 1, 2, 3, 0, 5, 5)
t.y = cumsum(t.y)
T.end = 76
n.y = length(t.y)


# MCMC
M = 20000
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
  u = X.to.u(t.x, t.y, t0, T.end, beta)
  k = ceiling((n.x + 1) * 0.03)
  u.new = u
  upi = sample(1:(n.x+1), k, replace = FALSE)
  u.new[upi] = rexp(k, 1)
  t.x.new = u.to.X(u.new, t.y, t0, T.end, beta)
  n.x.new = length(t.x.new)
  lp_new = Log.I_prod(t.x.new, t.y) - gamma * I_int(t.x.new, t.y, n.x.new, n.y, t0, T.end)
  lp_old = Log.I_prod(t.x, t.y) - gamma * I_int(t.x, t.y, n.x, n.y, t0, T.end)
  lp = lp_new - lp_old
  u = runif(1)
  if(is.nan(lp)){
    lp = -Inf
  }
  if(log(u) < lp){
    t.x = t.x.new
    n.x = n.x.new
    acc = acc + 1
  }
  
  
  ## update latent infection time by current MCMC algorithm
  # u = runif(1)
  # if (u < 1/3){
  #   # Add one time point
  #   n.x.new = n.x + 1
  #   t.new = runif(1, t0, T.end)
  #   t.x.new = sort(c(t.new, t.x))
  #   move_type = 1
  # }else{
  #   if (u > 2/3){
  #     # Move one time point
  #     n.x.new = n.x
  #     t.new = runif(1, t0, T.end)
  #     t.x.new = t.x[-sample(1:n.x, 1)]
  #     t.x.new = sort(c(t.new, t.x.new))
  #     move_type = 2
  #   }else{
  #     # Remove one time point
  #     n.x.new = n.x - 1
  #     t.x.new = t.x[-sample(1:n.x, 1)]
  #     move_type = 3
  #   }
  # }
  # lp_new = n.x.new * log(beta/N) + Log.SI_prod(t.x.new, t.y) - beta/N * SI_int(t.x.new, t.y, n.x.new, n.y, t0, T.end) +
  #   Log.I_prod(t.x.new, t.y) - gamma * I_int(t.x.new, t.y, n.x.new, n.y, t0, T.end)
  # lp_old = n.x * log(beta/N) + Log.SI_prod(t.x, t.y) - beta/N * SI_int(t.x, t.y, n.x, n.y, t0, T.end) +
  #   Log.I_prod(t.x, t.y) - gamma * I_int(t.x, t.y, n.x, n.y, t0, T.end)
  # if (move_type == 1){
  #   # Add one time point
  #   lp = lp_new + log(T.end-t0) - (lp_old + log(n.x+1))
  # }else{
  #   if (move_type == 3){
  #     # Remove one time point
  #     lp = lp_new + log(n.x) - (lp_old + log(T.end-t0))
  #   }else{
  #     # Move one time point
  #     lp = lp_new - lp_old
  #   }
  # }
  # u = runif(1)
  # if(is.na(lp)){
  #   lp = -Inf
  # }
  # if(log(u) < lp){
  #   t.x = t.x.new
  #   n.x = n.x.new
  #   acc = acc + 1
  # }
  
  sample.nx[j] = n.x
  sample.t0[j] = t0
  sample.beta[j] = beta
  sample.gamma[j] = gamma
}
proc.time() - ptm
acc/M

## Plots for new MCMC algorithm
plot(sample.nx, type = 'l', main = expression(paste('Trace plot of ', n[x], ' with new MCMC algorithm')),
     ylab = expression(n[x]))
plot(sample.t0, type = 'l', main = expression(paste('Trace plot of ', t[0]^x, ' with new MCMC algorithm')),
     ylab = expression(t[0]^x))
plot(sample.beta, type = 'l', main = expression(paste('Trace plot of ', beta, ' with new MCMC algorithm')),
     ylab = expression(beta))
plot(sample.gamma, type = 'l', main = expression(paste('Trace plot of ', gamma, ' with new MCMC algorithm')),
     ylab = expression(gamma))
hist(sample.beta, freq = FALSE, main = expression(paste('Histogram of ', beta, ' with new MCMC algorithm')),
     xlab = expression(beta))
hist(sample.gamma, freq = FALSE, main = expression(paste('Histogram of ', gamma, ' with new MCMC algorithm')),
     xlab = expression(gamma))
mean(sample.beta)
mean(sample.gamma)
var(sample.beta)
var(sample.gamma)
ess(sample.beta)
ess(sample.gamma)


## Plots for current MCMC algorithm
# plot(sample.nx, type = 'l', main = expression(paste('Trace plot of ', n[x], ' with current MCMC algorithm')),
#      ylab = expression(n[x]))
# plot(sample.t0, type = 'l', main = expression(paste('Trace plot of ', t[0]^x, ' with current MCMC algorithm')),
#      ylab = expression(t[0]^x))
# plot(sample.beta, type = 'l', main = expression(paste('Trace plot of ', beta, ' with current MCMC algorithm')),
#      ylab = expression(beta))
# plot(sample.gamma, type = 'l', main = expression(paste('Trace plot of ', gamma, ' with current MCMC algorithm')),
#      ylab = expression(gamma))
# hist(sample.beta, freq = FALSE, main = expression(paste('Histogram of ', beta, ' with current MCMC algorithm')),
#      xlab = expression(beta))
# hist(sample.gamma, freq = FALSE, main = expression(paste('Histogram of ', gamma, ' with current MCMC algorithm')),
#      xlab = expression(gamma))
# mean(sample.beta)
# mean(sample.gamma)
# var(sample.beta)
# var(sample.gamma)
# ess(sample.beta)
# ess(sample.gamma)


