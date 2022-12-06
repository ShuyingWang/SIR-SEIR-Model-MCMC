### Description of the code file ###
# The first part defines the functions.
# The second part includes the second illustration with simulated data (Section 4.2 in the paper).
# The simulated data can be generated with the seeds set.
# Comment/uncomment the new/current MCMC algorithms to run comparisons.


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


# New MCMC algorithm Proposal

New_proposal.y <- function(t.x, n.x, gamma, t0, T.end){
  tx = c(t.x, T.end)
  I = I0
  t = t0
  ty.new = c()
  i = 1
  
  while (t < T.end & I > 0){
    t_new = t + rexp(1, gamma * I)
    if (I==1 & i<n.x & t_new<tx[i]){
      I = I0
      t = t0
      ty.new = c()
      i = 1
    }else if (t_new < tx[i]){
      t = t_new
      ty.new = c(ty.new, t)
      I = I - 1
    }else{
      t = tx[i]
      I = I + 1
      i = i + 1
    }
  }
  return(ty.new)
}


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



### 2. Illustration using simulated data ###

library(mcmcse)
set.seed(2)

# Data simulation
I0 = 100
N = 1000000
S0 = N - I0
beta = 0.25
gamma = 0.15
t0 = 0
T.end = 10

t = Sim_SIR(N, I0, beta, gamma, T.end)
t.x = t$x
t.y = t$y
n.x = length(t.x)
ny.true = length(t.y)

# MCMC 
M = 10000
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
  # ty.new = New_proposal.y(t.x, n.x, gamma, t0, T.end)
  # ny.new = length(ty.new)
  # log_p.new = Log.SI_prod(t.x, ty.new) - beta/N * SI_int(t.x, ty.new, n.x, ny.new, t0, T.end)
  # log_p.old = Log.SI_prod(t.x, t.y) - beta/N * SI_int(t.x, t.y, n.x, n.y, t0, T.end)
  # log_p = log_p.new - log_p.old
  # u = runif(1)
  # if(log(u) < log_p){
  #   t.y = ty.new
  #   n.y = ny.new
  #   acc = acc + 1
  # }
  
  
  # update latent recovery times by current MCMC algorithm (Random Walk Metropolis)
  valid_proposal = FALSE
  while (valid_proposal == FALSE){
    u = runif(1, 0, 1)
    if (u < 1/3){
      # Add one time point
      ny.new = n.y + 1
      t.new = runif(1, t0, T.end)
      ty.new = sort(c(t.new, t.y))
      move_type = 1
    }else{
      if (u > 2/3){
        # Move one time point
        ny.new = n.y
        t.new = runif(1, t0, T.end)
        ty.new = t.y[-sample(1:n.y, 1)]
        ty.new = sort(c(t.new, ty.new))
        move_type = 2
      }else{
        # Remove one time point
        ny.new = n.y - 1
        ty.new = t.y[-sample(1:n.y, 1)]
        move_type = 3
      }
    }
    valid_proposal = Indicator(t.x, ty.new)
  }
  
  log_p.new = Log.SI_prod(t.x, ty.new) - beta/N * SI_int(t.x, ty.new, n.x, ny.new, t0, T.end) +
    ny.new * log(gamma) + Log.I_prod(t.x, ty.new) - gamma * I_int(t.x, ty.new, n.x, ny.new, t0, T.end)
  log_p.old = Log.SI_prod(t.x, t.y) - beta/N * SI_int(t.x, t.y, n.x, n.y, t0, T.end) +
    n.y * log(gamma) + Log.I_prod(t.x, t.y) - gamma * I_int(t.x, t.y, n.x, n.y, t0, T.end)
  if (move_type == 1){
    # Add one time point
    log_p = log_p.new + log(T.end-t0) - (log_p.old + log(n.y+1))
  }else{
    if (move_type == 3){
      # Remove one time point
      log_p = log_p.new + log(n.y) - (log_p.old + log(T.end-t0))
    }else{
      # Move one time point
      log_p = log_p.new - log_p.old
    }
  }
  u = runif(1)
  if(log(u) < log_p){
    t.y = ty.new
    n.y = ny.new
    acc = acc + 1
  }
}
proc.time() - ptm

acc.new = acc/M #0.5466
acc.current = acc/M #0.8776

Samples.new = list(ny=sample.ny, beta=sample.beta, gamma=sample.gamma)
save(Samples.new, file='New.RData')
Samples.current = list(ny=sample.ny, beta=sample.beta, gamma=sample.gamma)
save(Samples.current, file='Current.RData')


## Plots 
par(mfrow=c(2, 1))
plot(Samples.new$ny, type = 'l',
     main = expression(paste('Trace plot of ', n[y], ' with new MCMC algorithm')),
     ylab = expression(n[y]), ylim=c(100, 400))
plot(Samples.current$ny, type = 'l',
     main = expression(paste('Trace plot of ', n[y], ' with current MCMC algorithm')),
     ylab = expression(n[y]), ylim=c(100, 400))

par(mfrow=c(2, 1))
plot(Samples.new$beta, type = 'l',
     main = expression(paste('Trace plot of ', beta, ' with new MCMC algorithm')),
     ylab = expression(beta), ylim=c(0.15, 0.7))
plot(Samples.current$beta, type = 'l',
     main = expression(paste('Trace plot of ', beta, ' with current MCMC algorithm')),
     ylab = expression(beta), ylim=c(0.15, 0.7))

par(mfrow=c(2, 1))
plot(Samples.new$gamma, type = 'l',
     main = expression(paste('Trace plot of ', gamma, ' with new MCMC algorithm')),
     ylab = expression(gamma), ylim=c(0, 0.8))
plot(Samples.current$gamma, type = 'l',
     main = expression(paste('Trace plot of ', gamma, ' with current MCMC algorithm')),
     ylab = expression(gamma), ylim=c(0, 0.8))

par(mfrow=c(2, 1))
hist(Samples.new$ny[1000:M], freq=FALSE,
     main = expression(paste('Hist of ', n[y], ' with new MCMC')),
     xlab = expression(n[y]), xlim=c(100, 370), ylim=c(0, 0.025))
hist(Samples.current$ny[1000:M], freq=FALSE,
     main = expression(paste('Hist of ', n[y], ' with current MCMC')),
     xlab = expression(n[y]), xlim=c(100, 370), ylim=c(0, 0.025))

par(mfrow=c(2, 1))
hist(Samples.new$beta[1000:M], freq=FALSE,
     main = expression(paste('Histogram of ', beta, ' with new MCMC')),
     xlab = expression(beta), xlim=c(0.15, 0.55), ylim=c(0, 11))
hist(Samples.current$beta[1000:M], freq=FALSE,
     main = expression(paste('Histogram of ', beta, ' with current MCMC')),
     xlab = expression(beta), xlim=c(0.15, 0.55), ylim=c(0, 11))

par(mfrow=c(2, 1))
hist(Samples.new$gamma[1000:M], freq=FALSE,
     main = expression(paste('Histogram of ', gamma, ' with new MCMC')),
     xlab = expression(gamma), xlim=c(0.03, 0.55), ylim=c(0, 9))
hist(Samples.current$gamma[1000:M], freq=FALSE,
     main = expression(paste('Histogram of ', gamma, ' with current MCMC')),
     xlab = expression(gamma), xlim=c(0.03, 0.55), ylim=c(0, 9))

## Posterior means
mean(Samples.new$ny[1000:M]) #240.25
mean(Samples.current$ny[1000:M]) #329.58

mean(Samples.new$beta[1000:M]) #0.2397
mean(Samples.current$beta[1000:M]) #0.3525

mean(Samples.new$gamma[1000:M]) #0.1675
mean(Samples.current$gamma[1000:M]) #0.3314


