# Bayesian Gaussian regression on the radiata pine

chosen.model <- 1 # Model 1 or 2
nb.runs <- 10 # Number of independent runs 
choice <- 1 # choice=1, independent runs on each phase to measure the variability of each increment 
            # choice=2, independent runs on the entire algorithm to measure the variability of log(normalizing constant)

# Read the data

RadiataPine = read.table("RadiataPine.txt", sep = " ", header = TRUE)

y = RadiataPine$y

n = length(y)

X1 = cbind(rep(1, n), RadiataPine$x - mean(RadiataPine$x))

X2 = cbind(rep(1, n), RadiataPine$z - mean(RadiataPine$z))

d = ncol(X1)

tau <- 10 ^ (-5)

mu0 = c(3000, 185)
Q0 <- diag(x = c(0.06, 6.00))

# Choose the model

if (chosen.model == 1) {
  X1 = X1
} else {
  X1 = X2
}

# Definition of U, gradient U, hessian U

log.pi = function(theta) {
  term1 <- norm(y - X1 %*% theta, type = "2") ^ 2
  term2 <- t(theta - mu0) %*% Q0 %*% (theta - mu0)
  return(as.numeric(-(tau / 2) * (term1 + term2)))
}

gradient.log.pi = function(theta) {
  term1 <- (-1) * as.vector(t(X1) %*% y)
  term2 <- as.vector(t(X1) %*% X1 %*% theta)
  term3 <- as.vector(Q0 %*% (theta - mu0))
  return(-tau * (term1 + term2 + term3))
}

hessian.log.pi = function(theta) {
  term1 <- t(X1) %*% X1 + Q0
  return(-tau * term1)
}

# Exact Value Expected

temp1 <-
  norm(y, type = "2") ^ 2 + t(mu0) %*% Q0 %*% mu0 - t(t(X1) %*% y + Q0 %*% mu0) %*%
  solve(t(X1) %*% X1 + Q0) %*% (t(X1) %*% y + Q0 %*% mu0)
temp1 <- as.numeric(temp1)
temp2 <-
  (-d / 2) * log(2 * pi) - 0.5 * log(det(diag(d) + solve(Q0) %*% t(X1) %*%
                                           X1)) + (d / 2) * log(tau)
lnZ.Exact <- (-tau / 2) * temp1 + temp2

# Preliminary calculations

A <- t(X1) %*% X1 + Q0
B <- Q0 %*% mu0 + t(X1) %*% y
argmin <- solve(A, B)
log.pi.max <- log.pi(argmin)

vp = eigen(t(X1) %*% X1 + Q0, symmetric = T, only.values = T)

L = tau * max(vp$values)
m = tau * min(vp$values)
d = ncol(X1)

epsilon = 0.1
mu = 0.1

sigma2.i <- sigma02 <- (2 * log(1 + epsilon / 3)) / (d * (L - m))

lnZ = (d / 2) * log(2 * pi * sigma02) - (d / 2) * log(1 + sigma02 * m) +
  log.pi.max + d * log(tau / (2 * pi)) + 0.5 * log(det(Q0))



# Count the number of phases

sigma2.i <- sigma02
compt <- 1
sigma2.tab <- NULL

while (sigma2.i < (2 * d + 7) / m) {
  sigma2.tab[compt] = sigma2.i
  m.i = m + 1 / sigma2.i
  a.i = m.i / (4 * (d + 4))
  sigma2.i = (1 / sigma2.i - 2 * a.i) ^ (-1)
  compt <- compt + 1
}

sigma2.tab <- sigma2.tab[-length(sigma2.tab)]

# Save for analysis

analysis.tab <- matrix(nrow = length(sigma2.tab), ncol = 5)
analysis.tab[, 1] <- 1:nrow(analysis.tab)
colnames(analysis.tab) <-
  c("iterations", "sigma2", "mean", "standard deviation", "lnZ")

# Save the parameters for each phase

parameters.tab <- matrix(nrow = length(sigma2.tab), ncol = 6)
parameters.tab[, 1] <- 1:nrow(parameters.tab)
colnames(parameters.tab) <-
  c("iterations", "sigma2", "N.i", "n.i", "gam.i", "a.i")

for (k in 1:length(sigma2.tab)) {
  sigma2.i <- sigma2.tab[k]
  m.i = m + 1 / sigma2.i
  L.i = L + 1 / sigma2.i
  kappa.i = (2 * m.i * L.i) / (m.i + L.i)
  if (k < length(sigma2.tab)) {
    a.i = m.i / (4 * (d + 4))
    gam.i = (kappa.i * sigma2.i * m.i) / (d * L.i ^ 2)
    n.i = m.i ^ (1 / 2) / (kappa.i ^ 2 * sigma2.i ^ (1 / 2) * gam.i)
    N.i = (kappa.i * gam.i) ^ (-1)
  }
  else {
    a.i = 1 / (2 * sigma2.i)
    gam.i = kappa.i / (L.i ^ 2)
    N.i = n.i = (kappa.i * gam.i) ^ (-1)
  }
  parameters.tab[k, 2] <- sigma2.i
  parameters.tab[k, 3] <- N.i
  parameters.tab[k, 4] <- n.i
  parameters.tab[k, 5] <- gam.i
  parameters.tab[k, 6] <- a.i
}

# Adjustment of constants

parameters.tab.adjusted <- parameters.tab
parameters.tab.adjusted[, 3] = 10 ^ 3 * parameters.tab.adjusted[, 3]
parameters.tab.adjusted[, 4] = 10 ^ 4 * parameters.tab.adjusted[, 4]
parameters.tab.adjusted[, 5] = 10 ^ (-2) * parameters.tab.adjusted[, 5]

# Save a useful calculation for gradient U
term1 <- (-1) * as.vector(t(X1) %*% y)

# Algorithm

if (choice == 1) {
  
  # nb.runs independent runs for each phase to measure the variability of each increment to log(Z)
  
  for (k in 1:length(sigma2.tab)) {
    a.i = parameters.tab.adjusted[k, 6]
    sigma2.i = parameters.tab.adjusted[k, 2]
    N.i = as.integer(parameters.tab.adjusted[k, 3])
    n.i = as.integer(parameters.tab.adjusted[k, 4])
    gam.i = parameters.tab.adjusted[k, 5]
    
    chain <- rep(NA, nb.runs)
    
    for (l in 1:nb.runs) {
      piHat.i = 0
      theta = rep(0, d)
      
      for (j in 1:N.i) {
        X = theta + argmin
        term2 <- as.vector(t(X1) %*% X1 %*% X)
        term3 <- as.vector(Q0 %*% (X - mu0))
        gradient.U <- tau * (term1 + term2 + term3)
        gradient.U.i = theta / sigma2.i + gradient.U
        theta = theta - gam.i * gradient.U.i + sqrt(2 * gam.i) * rnorm(d)
      }
      
      for (j in (1 + N.i):(n.i + N.i)) {
        X = theta + argmin
        term2 <- as.vector(t(X1) %*% X1 %*% X)
        term3 <- as.vector(Q0 %*% (X - mu0))
        gradient.U <- tau * (term1 + term2 + term3)
        gradient.U.i = theta / sigma2.i + gradient.U
        theta = theta - gam.i * gradient.U.i + sqrt(2 * gam.i) * rnorm(d)
        piHat.i = piHat.i + exp(a.i * norm(theta, type = "2") ^ 2)
      }
      
      chain[l] = log(piHat.i / n.i)
      print(paste("END phase ", l))
      
    }
    
    mean.chain <- mean(chain)
    std.chain <- sd(chain)
    
    lnZ = lnZ + mean.chain
    print(paste("lnZ = ", lnZ))
    print(paste("mean increment = ", mean.chain))
    print(paste("standard deviation = ", std.chain))
    
    analysis.tab[k, 2] = sigma2.i
    analysis.tab[k, 3] = mean.chain
    analysis.tab[k, 4] = std.chain
    analysis.tab[k, 5] = lnZ
    
    print(paste("iteration = ", k))
    
  }
} else {
  
  # nb.runs independent runs of the entire algorithm to measure the variability of log(Z)
  
  lnZ.tab <- rep(NA, nb.runs)
  
  for (l in 1:nb.runs) {
    lnZ = (d / 2) * log(2 * pi * sigma02) - (d / 2) * log(1 + sigma02 * m) +
      log.pi.max + d * log(tau / (2 * pi)) + 0.5 * log(det(Q0))
    
    for (k in 1:length(sigma2.tab)) {
      a.i = parameters.tab.adjusted[k, 6]
      sigma2.i = parameters.tab.adjusted[k, 2]
      N.i = as.integer(parameters.tab.adjusted[k, 3])
      n.i = as.integer(parameters.tab.adjusted[k, 4])
      gam.i = parameters.tab.adjusted[k, 5]
      
      piHat.i = 0
      theta = rep(0, d)
      
      for (j in 1:N.i) {
        X = theta + argmin
        term2 <- as.vector(t(X1) %*% X1 %*% X)
        term3 <- as.vector(Q0 %*% (X - mu0))
        gradient.U <- tau * (term1 + term2 + term3)
        gradient.U.i = theta / sigma2.i + gradient.U
        theta = theta - gam.i * gradient.U.i + sqrt(2 * gam.i) * rnorm(d)
      }
      
      for (j in (1 + N.i):(n.i + N.i)) {
        X = theta + argmin
        term2 <- as.vector(t(X1) %*% X1 %*% X)
        term3 <- as.vector(Q0 %*% (X - mu0))
        gradient.U <- tau * (term1 + term2 + term3)
        gradient.U.i = theta / sigma2.i + gradient.U
        theta = theta - gam.i * gradient.U.i + sqrt(2 * gam.i) * rnorm(d)
        piHat.i = piHat.i + exp(a.i * norm(theta, type = "2") ^ 2)
      }
      
      increment = log(piHat.i / n.i)
      
      
      
      
      lnZ = lnZ + increment
      print(paste("lnZ = ", lnZ))
      print(paste("increment = ", increment))
      
      print(paste("iteration = ", k))
    }
    
    lnZ.tab[l] = lnZ
    print(paste("END phase ", l))
    
  }
}

save.image(file = paste(
  "save-reg-gaus-model",
  chosen.model,
  "-",
  choice,
  ".RData",
  sep = ""
))
