### Description of the code file ###
# The first part defines the functions.
# The second part includes the illustration with COVID-19 data.
# Comment the new MCMC algorithm and uncomment the current MCMC algorithm to run comparisons.


### Define functions ###
## Integrals 
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




## Products
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




## Generate
Generate.x <- function(beta, nx, t.y, t.z, t1, t2){
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





### Covid Data ###
N = 285000
E0 = 15
I0 = 3
S0 = N - I0 - E0
t0 = 0
T.end = 35
k = 15
beta1 = 0.39
beta2 = 0.075
alpha = 0.15
gamma = 0.05

Y = c()
Z = c()

for (i in 1:1000){
  data = Sim_SEIR(N, E0, I0, beta1, alpha, gamma)
  time_points = cumsum(data$t)
  type = data$type
  
  t.x.1 = time_points[type==1]
  t.y.1 = time_points[type==2]
  t.z.1 = time_points[type==3]
  t.x.1 = t.x.1[t.x.1 < k]
  t.y.1 = t.y.1[t.y.1 < k]
  t.z.1 = t.z.1[t.z.1 < k]
  n.x.1 = length(t.x.1)
  n.y.1 = length(t.y.1)
  n.z.1 = length(t.z.1)
  
  E = E0 + n.x.1 - n.y.1
  I = I0 + n.y.1 - n.z.1
  
  data = Sim_SEIR(N, E, I, beta2, alpha, gamma)
  time_points = cumsum(data$t)
  type = data$type
  
  t.x.2 = time_points[type==1]
  t.y.2 = time_points[type==2]
  t.z.2 = time_points[type==3]
  
  t.x.2 = t.x.2[t.x.2 < (T.end - k)] + k
  n.x.2 = length(t.x.2)
  t.y.2 = t.y.2[t.y.2 < (T.end - k)] + k
  n.y.2 = length(t.y.2)
  t.z.2 = t.z.2[t.z.2 < (T.end - k)] + k
  n.z.2 = length(t.z.2)
  
  t.x = c(t.x.1, t.x.2)
  t.y = c(t.y.1, t.y.2)
  t.z = c(t.z.1, t.z.2)
  nx.true = n.x.1 + n.x.2
  n.x = nx.true
  Y = append(Y, n.y.1 + n.y.2)
  Z = append(Z, n.z.1 + n.z.2)
}

mean(Z)
mean(Y)
hist(Y)
hist(Z)




library(R0) # Loading package

data = infection[16:35]
mGT<-generation.time("gamma", c(2, 0.5))
estR0<-estimate.R(data, mGT, begin=1, end=20, methods=c("EG", "ML", "TD", "AR", "SB"), 
                  pop.size=N, nsim=200)

attributes(estR0)
estR0





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

t.x.1 = Generate.x(0.3, 0, t.y, t.z, t0, k)
n.x.1 = length(t.x.1)
t.x.2 = Generate.x(0.15, n.x.1, t.y, t.z, k, T.end)
n.x.2 = length(t.x.2)
t.x = c(t.x.1, t.x.2)
n.x = n.x.1 + n.x.2

### MCMC ###
M = 10000
sample.nx = c()
sample.beta1 = c()
sample.beta2 = c()
sample.alpha = c()
sample.gamma = c()
acc1 = 0
acc2 = 0
acc3 = 0
acc4 = 0
ptm <- proc.time()
for (j in 1:M){
  ## Update parameters
  beta1 = rgamma(1, n.x.1+20, SI.int(t.x.1, t.y, t.z, t0, k)/N+100)
  beta2 = rgamma(1, n.x.2+20, SI.int(t.x, t.y, t.z, k, T.end)/N+100)
  alpha = rgamma(1, n.y+100, E_int(t.x, t.y)+600)
  gamma = rgamma(1, n.z+100, I_int(t.y, t.z)+1000)
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
      acc3 = acc3 + 1
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
      acc4 = acc4 + 1
    }
  }
  t.z = t.z[2:(n.z+1)]



  
  
  ## Update latent X by new MCMC algorithm
  # tx1.new = Generate.x(beta1, 0, t.y, t.z, t0, k)
  # t.x.new = c(tx1.new, t.x.2)
  # 
  # lp_new = Log.SI.prod(t.x.new, t.y, t.z, k, T.end) - beta2/N * SI.int(t.x.new, t.y, t.z, k, T.end) +
  #   Log.E_prod(t.x.new, t.y) - alpha * E_int(t.x.new, t.y)
  # 
  # lp_old = Log.SI.prod(t.x, t.y, t.z, k, T.end) - beta2/N * SI.int(t.x, t.y, t.z, k, T.end) +
  #   Log.E_prod(t.x, t.y) - alpha * E_int(t.x, t.y)
  # 
  # lp = lp_new - lp_old
  # u = runif(1)
  # if(is.na(lp)){
  #   lp = -Inf
  # }
  # if(log(u) < lp){
  #   t.x.1 = tx1.new
  #   t.x = t.x.new
  #   n.x.1 = length(t.x.1)
  #   n.x = n.x.1 + n.x.2
  #   acc1 = acc1 + 1
  # }
  # 
  # tx2.new = Generate.x(beta2, n.x.1, t.y, t.z, k, T.end)
  # t.x.new = c(t.x.1, tx2.new)
  # 
  # lp_new = Log.E.prod(t.x.new, t.y, k, T.end) - alpha * E.int(t.x.new, t.y, k, T.end)
  # lp_old = Log.E.prod(t.x, t.y, k, T.end) - alpha * E.int(t.x, t.y, k, T.end)
  # 
  # lp = lp_new - lp_old
  # u = runif(1)
  # if(is.na(lp)){
  #   lp = -Inf
  # }
  # if(log(u) < lp){
  #   t.x.2 = tx2.new
  #   t.x = t.x.new
  #   n.x.2 = length(t.x.2)
  #   n.x = n.x.1 + n.x.2
  #   acc2 = acc2 + 1
  # }
  # 
  # sample.nx[j] = n.x
  
  
  
  
  
  ## Update latent X by current MCMC algorithm
  u = runif(1, 0, 1)
  if (u < 1/3){
    # Add one time point
    nx.new = n.x + 1
    t.new = runif(1, t0, T.end)
    t.x.new = sort(c(t.new, t.x))
    move_type = 1
  }else{
    if (u > 2/3){
      # Move one time point
      nx.new = n.x
      t.new = runif(1, t0, T.end)
      t.x.new = t.x[-sample(1:n.x, 1)]
      t.x.new = sort(c(t.new, t.x.new))
      move_type = 2
    }else{
      # Remove one time point
      nx.new = n.x - 1
      t.x.new = t.x[-sample(1:n.x, 1)]
      move_type = 3
    }
  }
  tx1.new = t.x.new[t.x.new<k]
  nx1.new = length(tx1.new)
  tx2.new = t.x.new[t.x.new>k]
  nx2.new = length(tx2.new)

  lp_new = Log.E_prod(t.x.new, t.y) - alpha * E_int(t.x.new, t.y) +
    nx1.new*log(beta1/N) + nx2.new*log(beta2/N) + Log.SI_prod(t.x.new, t.y, t.z) -
    beta1/N * SI.int(tx1.new, t.y, t.z, t0, k) - beta2/N * SI.int(t.x.new, t.y, t.z, k, T.end)

  lp_old = Log.E_prod(t.x, t.y) - alpha * E_int(t.x, t.y) +
    n.x.1*log(beta1/N) + n.x.2*log(beta2/N) + Log.SI_prod(t.x, t.y, t.z) -
    beta1/N * SI.int(t.x.1, t.y, t.z, t0, k) - beta2/N * SI.int(t.x, t.y, t.z, k, T.end)

  
  if (move_type == 1){
    # Add one time point
    lp = lp_new + log(T.end-t0) - (lp_old + log(n.x+1))
  }else{
    if (move_type == 3){
      # Remove one time point
      lp = lp_new + log(n.x) - (lp_old + log(T.end-t0))
    }else{
      # Move one time point
      lp = lp_new - lp_old
    }
  }
  u = runif(1)
  if(is.na(lp)){
    lp = -Inf
  }
  if(log(u) < lp){
    t.x = t.x.new
    t.x.1 = tx1.new
    n.x.1 = nx1.new
    t.x.2 = tx2.new
    n.x.2 = nx2.new
    n.x = nx.new
    acc1 = acc1 + 1
  }

  sample.nx[j] = n.x
}
proc.time() - ptm
acc1/M
acc2/M

plot(sample.nx, type = 'l',
     main = expression(paste('Trace plot of ', n[x], ' with current MCMC algorithm')),
     ylab = expression(n[x]))

plot(sample.beta1, type = 'l', main = expression(paste('Trace plot of ', beta[1], ' with current MCMC algorithm')),
     ylab = expression(beta[1]))
plot(sample.beta2, type = 'l', main = expression(paste('Trace plot of ', beta[2], ' with current MCMC algorithm')),
     ylab = expression(beta[2]))
plot(sample.alpha, type = 'l', main = expression(paste('Trace plot of ', alpha, ' with current MCMC algorithm')),
     ylab = expression(alpha))
plot(sample.gamma, type = 'l', main = expression(paste('Trace plot of ', gamma, ' with current MCMC algorithm')),
     ylab = expression(gamma))
mean(sample.beta1[5000:M])
mean(sample.beta2[5000:M])
mean(sample.alpha[5000:M])
mean(sample.gamma[5000:M])

