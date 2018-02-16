# Multivariate Gaussian in dimension 10, 25, 50

d <- 10 # Dimension
nb.runs <- 10 # Number of independent runs
nb.jumps <- 1 # Number of jumped phases to accelerate the algorithm, set to 1 in the simulations.
choice <- 2 # choice=1, independent runs on each phase to measure the variability of each increment 
            # choice=2, independent runs on the entire algorithm to measure the variability of log(normalizing constant)

L = 2
m = 1
invcov = c(2, rep(1, d - 1))

lnZ.Exact <- (d / 2) * log(2 * pi) - 0.5 * log(2) # Exact Value Expected

epsilon = 0.1
mu = 0.1

sigma2.i <- sigma02 <- (2 * log(1 + epsilon / 3)) / (d * (L - m))

lnZ = (d / 2) * log(2 * pi * sigma02) - (d / 2) * log(1 + sigma02 * m)

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

# Parameters of the algorithm

parameters.tab <- matrix(nrow = length(sigma2.tab), ncol = 3)
parameters.tab[, 1] <- 1:nrow(parameters.tab)
colnames(parameters.tab) <-
  c("iterations", "sigma2", "a.i")

for (k in 1:length(sigma2.tab)) {
  sigma2.i <- sigma2.tab[k]
  m.i = m + 1 / sigma2.i
  L.i = L + 1 / sigma2.i
  if (k < length(sigma2.tab)) {
    a.i = m.i / (4 * (d + 4))
  }
  else {
    a.i = 1 / (2 * sigma2.i)
  }
  parameters.tab[k, 2] <- sigma2.i
  parameters.tab[k, 3] <- a.i
}

n = length(sigma2.tab)
q = n %/% nb.jumps

# Algorithm

if (choice == 1) {
  
  # nb.runs independent runs for each phase to measure the variability of each increment to log(Z)
  
  for (k in 1:q) {
    sigma2.i = parameters.tab[k * nb.jumps, 2]
    if (k < q) {
      sigma2.iplus1 = parameters.tab[(k + 1) * nb.jumps, 2]
    } else {
      sigma2.iplus1 = Inf
    }
    m.i = m + 1/sigma2.i
    L.i = L + 1/sigma2.i
    
    N.i = 10 ^ 4
    gam.i = 10^(-2) * (m.i + L.i) ^ (-1)
    n.i = 10 ^ 5
    
    if (k < q) {
      a.i = 0.5 * (sigma2.i ^ (-1) - sigma2.iplus1 ^ (-1))
    } else{
      a.i = 2 * parameters.tab[nrow(parameters.tab) - 1, "a.i"]
    }
    
    chain <- rep(NA, nb.runs)
    
    for (l in 1:nb.runs) {
      piHat.i = 0
      theta = rep(0, d)
      
      for (j in 1:N.i) {
        gradient.U <- invcov * theta
        gradient.U.i = theta / sigma2.i + gradient.U
        theta = theta - gam.i * gradient.U.i + sqrt(2 * gam.i) * rnorm(d)
      }
      
      for (j in (1 + N.i):(n.i + N.i)) {
        gradient.U <- invcov * theta
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
    lnZ = (d / 2) * log(2 * pi * sigma02) - (d / 2) * log(1 + sigma02 * m)
    
    for (k in 1:q) {
      sigma2.i = parameters.tab[k * nb.jumps, 2]
      if (k < q) {
        sigma2.iplus1 = parameters.tab[(k + 1) * nb.jumps, 2]
      } else {
        sigma2.iplus1 = Inf
      }
      m.i = m + 1/sigma2.i
      L.i = L + 1/sigma2.i
      
      N.i = 10 ^ 4
      gam.i = 10^(-2) * (m.i + L.i) ^ (-1)
      n.i = 10 ^ 5
      
      if (k < q) {
        a.i = 0.5 * (sigma2.i ^ (-1) - sigma2.iplus1 ^ (-1))
      } else{
        a.i = 2 * parameters.tab[nrow(parameters.tab) - 1, "a.i"]
      }
      
      piHat.i = 0
      theta = rep(0, d)
      
      for (j in 1:N.i) {
        gradient.U <- invcov * theta
        gradient.U.i = theta / sigma2.i + gradient.U
        theta = theta - gam.i * gradient.U.i + sqrt(2 * gam.i) * rnorm(d)
      }
      
      for (j in (1 + N.i):(n.i + N.i)) {
        gradient.U <- invcov * theta
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

save.image(file = paste("save-gaus-d", d, "-", choice, ".RData", sep = ""))