### Description of the code file ###
# The first part defines the functions.
# The second part includes the first illustration with simulated data (Section 4.1 in the paper).
# The simulated data can be generated with the seeds set.

library(mcmcse)
library(survival)
library(nbpMatching)

### 1. Functions ###
# Data Simulation
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
New_MH.x <- function(t.x, n.x, t.y, n.y, beta, gamma, t0, T.end, log_p){
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
  nx.new = length(tx.new)
  log_p.new = Log.I_prod(tx.new, t.y) - gamma * I_int(tx.new, t.y, nx.new, n.y, t0, T.end)
  acc = 0
  if (log(runif(1)) < (log_p.new-log_p)){
    t.x = tx.new
    n.x = nx.new
    log_p = log_p.new
    acc = 1
  }
  return(list("tx"=t.x, "nx"=n.x, "lp"=log_p, "acc"=acc))
}
New_MH_adjust.x <- function(r, t.x, n.x, t.y, n.y, beta, gamma, t0, T.end, log_p){
  u = X_to_u(t.x, t.y, t0, T.end, beta)
  u.new = u_to_u.new(u, r)
  tx.new = u_to_X(u.new, t.y, t0, T.end, beta)
  nx.new = length(tx.new)
  acc = 0
  if (Indicator(tx.new, t.y)){
    log_p.new = Log.I_prod(tx.new, t.y) - gamma * I_int(tx.new, t.y, nx.new, n.y, t0, T.end)
    log_a = log_p.new-log_p
    acc = 0
    if (log(runif(1)) < log_a){
      t.x = tx.new
      n.x = nx.new
      log_p = log_p.new
      acc = 1
    }
  }
  return(list("tx"=t.x, "nx"=n.x, "lp"=log_p, "acc"=acc))
}

# Current MCMC algorithm (Random Walk Metropolis)
Current_MH.x <- function(k, t.x, n.x, t.y, n.y, beta, gamma, t0, T.end, log_p){
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
    log_p.new = nx.new * log(beta/N) + Log.SI_prod(tx.new, t.y) - beta/N * SI_int(tx.new, t.y, nx.new, n.y, t0, T.end) +
      Log.I_prod(tx.new, t.y) - gamma * I_int(tx.new, t.y, nx.new, n.y, t0, T.end)
    log_a = log_a + log_p.new - log_p
    if (log(runif(1)) < log_a){
      t.x = tx.new
      n.x = nx.new
      log_p = log_p.new
      acc = 1
    }
  }
  return(list("tx"=t.x, "nx"=n.x, "lp"=log_p, "acc"=acc))
}

# Simulation-based Bayesian inference for epidemic models
Simulation_proposal.x <- function(t.y, n.y, beta, t0, T.end){
  ty = c(t.y, T.end)
  S = S0
  I = I0
  t = t0
  tx.new = c()
  i = 1
  while (t < T.end & S > 0 & I > 0){
    t_new = t + rexp(1, beta/N * S * I)
    if (I == 1 & i < n.y){
      while (t_new > ty[i]){
        t_new = t + rexp(1, beta/N * S * I)
      }
      t = t_new
      tx.new = c(tx.new, t)
      I = I + 1
      S = S - 1
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
Simulation_MH.x <- function(t.x, n.x, t.y, n.y, beta, gamma, t0, T.end, log_p){
  ty = c(t.y, T.end)
  S = S0
  I = I0
  t = t0
  tx.new = c()
  i = 1
  lp.new = 0
  while (t < T.end & S > 0 & I > 0){
    t_new = t + rexp(1, beta/N * S * I)
    if (I == 1 & i < n.y){
      lp.new = lp.new + log(1 - exp(-beta*S*I/N * (ty[i]-t)))
      while (t_new > ty[i]){
        t_new = t + rexp(1, beta/N * S * I)
      }
      t = t_new
      tx.new = c(tx.new, t)
      I = I + 1
      S = S - 1
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
  nx.new = length(tx.new)
  log_p.new = lp.new + Log.I_prod(tx.new, t.y) - gamma * I_int(tx.new, t.y, nx.new, n.y, t0, T.end)
  acc = 0
  if (log(runif(1)) < (log_p.new-log_p)){
    t.x = tx.new
    n.x = nx.new
    log_p = log_p.new
    acc = 1
  }
  return(list("tx"=t.x, "nx"=n.x, "lp"=log_p, "acc"=acc))
}

# Sample distance
crossmatchtest <- function(z, D){
  if ( !isSymmetric(D) )
  {
    stop("Invalid distance matrix: your distance matrix is not symmetric")
    return(NA)
  }
  if ( sum(D < 0 ) > 0 )
  {
    stop("Invalid distance matrix: your distance matrix includes negative values")
    return(NA)
  }
  plainmatrix <- 100000*D/max(as.vector(D))
  diag(plainmatrix) <- 0
  mdm <- distancematrix(plainmatrix)
  nzero <- sum(mdm==0) - length(diag(mdm))
  if (nzero/((dim(mdm)[1])^2-dim(mdm)[1])>.95)
  {
    warning("Your distance matrix has some very large relative distances such that more than 95 percent of distances were rounded to zero")
  }
  res <- nonbimatch(mdm)
  mt <- pmin(as.numeric(res$matches$Group1.Row),as.numeric(res$matches$Group2.Row))
  if ( length(z) < length(mt) ) ##if the number of observations is odd remove observation that paired with ghost
  {
    mt[ mt==mt[length(mt)] ] <- 0
    mt <- mt[1:(length(mt)-1)]
  }
  z0 <- z[mt>0]
  mt0 <- factor(mt[mt>0])
  tab <- table(factor(z0),mt0)
  a1 <- sum(tab[1,]==1)
  bigN <- length(z0)
  n <- sum(z0)
  if (bigN<340) ##if the number of observations is below 340 compute the exact null distribution
  {
    dist <- crossmatchdist(bigN,n)
    pval <- dist[5,dist[2,]==a1]
  }
  else
  {
    pval<-NA
  }
  m <- bigN-n
  Ea1 <- (n*m/(bigN-1))
  Va1 <- 2*n*(n-1)*m*(m-1)/((bigN-3)*(bigN-1)*(bigN-1))
  dev <- (a1-Ea1)/sqrt(Va1)
  approx <- pnorm(dev)
  list (a1=a1,Ea1=Ea1,Va1=Va1,dev=dev,pval=pval,approxpval=approx)
}
Distance_matrix <- function(x, y){
  data = c(x, y)
  n = length(data)
  D = matrix(nrow = n, ncol = n)
  for (i in 1:n){
    for (j in 1:n){
      D[i, j] = abs(data[i] - data[j])
    }
  }
  return(D)
}





### 2. Illustrations with simulated data ###

## Setting
# Small
N = 1000
I0 = 10
S0 = N - I0
beta = 1.5
gamma = 1
t0 = 0
T.end = 3

# Moderate
N = 1000
I0 = 20
S0 = N - I0
beta = 2
gamma = 1
t0 = 0
T.end = 3

# Large
N = 1000
I0 = 50
S0 = N - I0
beta = 2.5
gamma = 1
t0 = 0
T.end = 5





## Simulate data
set.seed(3)
t = Sim_SIR(N, I0, beta, gamma, T.end)
t.x = t$x
t.y = t$y
nx.true = length(t.x)
n.y = length(t.y)
nx.true
n.y





## Compare proposal distributions
B = 10000
Sample_proposal.new = rep(NA, B)
Sample_proposal.simulation = rep(NA, B)
for (i in 1:B){
  Sample_proposal.new[i] = length(New_proposal.x(t.y, n.y, beta, 0, T.end))
  Sample_proposal.simulation[i] = length(Simulation_proposal.x(t.y, n.y, beta, 0, T.end))
  print(i)
}
par(mfrow = c(2, 1))
hist(Sample_proposal.new, xlab=expression(n[x]), freq=FALSE, breaks=20, main='New proposal')
hist(Sample_proposal.simulation, xlab=expression(n[x]), freq=FALSE, breaks=20, main='Simulation-based proposal')





## Compare MCMC samples
# New MCMC
set.seed(3)
B = 10000
n.x = n.y
t.x = sort(t.y - rexp(n.y, gamma))
t.x = ifelse(t.x<0, 0, t.x)
log_p = Log.I_prod(t.x, t.y) - gamma * I_int(t.x, t.y, n.x, n.y, t0, T.end)
sample.new = rep(NA, B)
acc.new = 0
# r = 0.1
# r = 0.03
for (i in 1:B){
  # Sample = New_MH.x(t.x, n.x, t.y, n.y, beta, gamma, t0, T.end, log_p)
  Sample = New_MH_adjust.x(r, t.x, n.x, t.y, n.y, beta, gamma, t0, T.end, log_p)
  t.x = Sample$tx
  n.x = Sample$nx
  log_p = Sample$lp
  acc.new = acc.new + Sample$acc
  sample.new[i] = n.x
  print(i)
}
acc.new/B # 0.1197, 0.1613, 0.3032
ess(sample.new) # 885.83, 549.16, 467.84
mean(sample.new) # 115.5, 418.9, 824.1
plot(sample.new, type='l', ylab=expression(n[x]), main='New MCMC algorithm')


# Simulation based MCMC
set.seed(3)
n.x = n.y
t.x = sort(t.y - rexp(n.y, gamma))
t.x = ifelse(t.x<0, 0, t.x)
log_p = Log.I_prod(t.x, t.y) - gamma * I_int(t.x, t.y, n.x, n.y, t0, T.end)
sample.simulation = rep(NA, B)
acc.simulation = 0
for (i in 1:B){
  Sample = Simulation_MH.x(t.x, n.x, t.y, n.y, beta, gamma, t0, T.end, log_p)
  t.x = Sample$tx
  n.x = Sample$nx
  log_p = Sample$lp
  acc.simulation = acc.simulation + Sample$acc
  sample.simulation[i] = n.x
  print(i)
}
acc.simulation/B # 0.0447, 0.0425, 0.0377
ess(sample.simulation) # 362.06, 199.18, 190.41
mean(sample.simulation) # 115.0, 419.5, 824.1
plot(sample.simulation, type='l', ylab=expression(n[x]), main='Simulation-based MCMC algorithm')


# Random walk MCMC
set.seed(3)
n.x = n.y
t.x = sort(t.y - rexp(n.y, gamma))
t.x = ifelse(t.x<0, 0, t.x)
log_p = n.x * log(beta/N) + Log.SI_prod(t.x, t.y) - beta/N * SI_int(t.x, t.y, n.x, n.y, t0, T.end) +
  Log.I_prod(t.x, t.y) - gamma * I_int(t.x, t.y, n.x, n.y, t0, T.end)
sample.current = rep(NA, B)
acc.current = 0
k = 20
# k = 20
# k = 5
for (i in 1:B){
  Sample = Current_MH.x(k, t.x, n.x, t.y, n.y, beta, gamma, t0, T.end, log_p)
  t.x = Sample$tx
  n.x = Sample$nx
  log_p = Sample$lp
  acc.current = acc.current + Sample$acc
  sample.current[i] = n.x
  print(i)
}
acc.current/B # 0.3114, 0.2289, 0.3137
ess(sample.current) # 63.95, 18.12, 22.58
mean(sample.current[1000:10000]) # 114.8, 411.1, 824.3
plot(sample.current, type='l', ylab=expression(n[x]), main='Current MCMC algorithm')





## Results
Samples.small = list(Sample_proposal.new, Sample_proposal.simulation, sample.new, sample.simulation, sample.current)
save(Samples.small, file='small.RData')

Samples.moderate = list(Sample_proposal.new, Sample_proposal.simulation, sample.new, sample.simulation, sample.current)
save(Samples.moderate, file='moderate.RData')

Samples.large = list(Sample_proposal.new, Sample_proposal.simulation, sample.new, sample.simulation, sample.current)
save(Samples.large, file='large.RData')


# Acceptance rates
acc.new/B # 0.1197, 0.1613, 0.3032
acc.simulation/B # 0.0447, 0.0425, 0.0377
acc.current/B # 0.3114, 0.2289, 0.3137


# Effective sample sizes
ess(sample.new) # 885.83, 549.16, 467.84
ess(sample.simulation) # 362.06, 199.18, 190.41
ess(sample.current) # 63.95, 18.12, 22.58


# Plots
par(mfrow = c(3, 1))
hist(Sample_proposal.new, xlab=expression(n[x]), freq=FALSE, breaks=20, main='New proposal', ylim=c(0, 0.05), xlim=c(765, 940))
hist(Sample_proposal.simulation, xlab=expression(n[x]), freq=FALSE, breaks=20, main='Simulation-based proposal', ylim=c(0, 0.05), xlim=c(765, 940))
hist(sample.new, breaks=20, xlab=expression(n[x]), freq=FALSE, ylim=c(0, 0.05), main='Target distribution', xlim=c(765, 940))

par(mfrow = c(3, 1))
plot(sample.new, type='l', ylab=expression(n[x]), main='New MCMC algorithm')
plot(sample.simulation, type='l', ylab=expression(n[x]), main='Simulation-based MCMC algorithm')
plot(sample.current, type='l', ylab=expression(n[x]), main='Current MCMC algorithm')

par(mfrow = c(3, 1))
acf(sample.new,  main='New MCMC algorithm')
acf(sample.simulation,  main='Simulation-based MCMC algorithm')
acf(sample.current, main='Current MCMC algorithm')


# Sample distance between proposal and target by cross match
label = c(rep(0, 2000), rep(1, 2000))
sample.target = sample.new[5*(1:2000)]
crossmatch.new = crossmatchtest(label, Distance_matrix(Sample_proposal.new[1:2000], sample.target))
1 - crossmatch.new$a1 / crossmatch.new$Ea1 #0.518, 0.592, 0.536
crossmatch.simulation = crossmatchtest(label, Distance_matrix(Sample_proposal.simulation[1:2000], sample.target))
1 - crossmatch.simulation$a1 / crossmatch.simulation$Ea1 #0.794, 0.788, 0.736





