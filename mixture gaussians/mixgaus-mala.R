# mixture of Gaussian distributions

nb.runs <- 10 # Number of independent runs
nb.jumps <- 5 # Number of jumped phases to accelerate the algorithm, set to 1 in the simulations.
choice <- 2 # choice=1, independent runs on each phase to measure the variability of each increment 
            # choice=2, independent runs on the entire algorithm to measure the variability of log(normalizing constant)

nb.mix <- 4 # Number of mixtures

# Read the data

y <- read.table("smcsampler_section4_data.txt")
y.max <- max(y)
y.min <- min(y)
xi <- (y.min + y.max) / 2
kappa <- (y.max - y.min) ^ (-2)
alpha <- 2
beta <- 2 * kappa
lambda <- alpha * beta
n <- nrow(y)

# potential U and its gradient

U <- function(mu) {
  temp <-
    -(n / 2) * log(lambda / (2 * pi)) - n * log(1 / nb.mix) - sum(log(apply(y, 1, function(x)
      sum(exp(-(lambda / 2) * (rep(x, nb.mix) - mu) ^ 2 ))))) 
    - (nb.mix/2) * log(kappa / (2 * pi)) + (kappa / 2) * sum((mu - rep(xi, nb.mix)) ^ 2)
  return(temp)
}

gradU <- function(mu) {
  temp <- apply(y, 1, function(x) lambda*(mu-rep(x,nb.mix))*exp(-(lambda/2)*(rep(x,nb.mix)-mu)^2) / sum(exp(-(lambda/2)*(rep(x,nb.mix)-mu)^2)))
  return(rowSums(temp) + kappa*(mu-rep(xi,nb.mix)))
}

# ml <- optim(rep(0,nb.mix), U, gr = gradU, method = "BFGS")
# argmin = ml$par
argmin = rep(1.76562, nb.mix)
Umin = U(argmin)

L = 1
m = kappa
d = nb.mix

epsilon = 0.1
mu = 0.1

sigma2.i <- sigma02 <- (2 * log(1 + epsilon / 3)) / (d * (L - m))

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

parameters.tab <- matrix(nrow = length(sigma2.tab), ncol = 3)
parameters.tab[, 1] <- 1:nrow(parameters.tab)
colnames(parameters.tab) <-
  c("iterations", "sigma2", "gam.i")

for (k in 1:length(sigma2.tab)) {
  sigma2.i <- sigma2.tab[k]
  m.i = m + 1 / sigma2.i
  L.i = L + 1 / sigma2.i
  kappa.i = (2 * m.i * L.i) / (m.i + L.i)
  if (k < length(sigma2.tab)) {
    gam.i = (kappa.i * sigma2.i * m.i) / (d * L.i ^ 2)
  }
  else {
    gam.i = kappa.i / (L.i ^ 2)
  }
  parameters.tab[k, 2] <- sigma2.i
  parameters.tab[k, 3] <- gam.i
}

parameters.tab.adjusted <- parameters.tab
parameters.tab.adjusted[, 3] = 10^(-1) * parameters.tab.adjusted[, 3]

lnZ <-
  (d / 2) * log(2 * pi * sigma02) - (d / 2) * log(1 + sigma02 * m) - Umin

q = length(sigma2.tab) %/% nb.jumps

# Algorithm

if (choice == 1) {
  for (k in 1:q) {
    sigma2.i = parameters.tab.adjusted[k * nb.jumps, 2]
    if (k < q) {
      sigma2.iplus1 = parameters.tab.adjusted[(k + 1) * nb.jumps, 2]
    } else {
      sigma2.iplus1 = Inf
    }
    gam.i = parameters.tab.adjusted[k * nb.jumps, 3]
    
    N.i = 10^4
    n.i = 10^6
    
    a.i = 0.5 * (sigma2.i ^ (-1) - sigma2.iplus1 ^ (-1))
    
    chain <- rep(NA, nb.runs)
    
    for (l in 1:nb.runs) {
      piHat.i = 0
      theta = rep(0, d)
      
      for (j in 1:N.i) {
        gradient.U <- gradU(theta + argmin) + theta / sigma2.i
        theta.prop = theta - gam.i * gradient.U + sqrt(2 * gam.i) * rnorm(d)
        gradient.U.prop <-
          gradU(theta.prop + argmin) + theta.prop / sigma2.i
        ratio <-
          -U(theta.prop + argmin) - (2 * sigma2.i) ^ (-1) * sum(theta.prop ^ 2) + U(theta + argmin) + (2 * sigma2.i) ^ (-1) * sum(theta ^ 2) + (4 * gam.i) ^ (-1) * (sum((theta.prop - theta + gam.i * gradient.U) ^ 2) - sum((theta - theta.prop + gam.i * gradient.U.prop) ^ 2))
        if (log(runif(1)) < ratio) {
          theta = theta.prop
        }
        
      }
      
      for (j in (1 + N.i):(n.i + N.i)) {
        gradient.U <- gradU(theta + argmin) + theta / sigma2.i
        theta.prop = theta - gam.i * gradient.U + sqrt(2 * gam.i) * rnorm(d)
        gradient.U.prop <-
          gradU(theta.prop + argmin) + theta.prop / sigma2.i
        ratio <-
          -U(theta.prop + argmin) - (2 * sigma2.i) ^ (-1) * sum(theta.prop ^ 2) + U(theta + argmin) + (2 * sigma2.i) ^ (-1) * sum(theta ^ 2) + (4 * gam.i) ^ (-1) * (sum((theta.prop - theta + gam.i * gradient.U) ^ 2) - sum((theta - theta.prop + gam.i * gradient.U.prop) ^ 2))
        if (log(runif(1)) < ratio) {
          theta = theta.prop
        }
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
  lnZ.tab <- rep(NA, nb.runs)
  
  for (l in 1:nb.runs) {
    lnZ = (d / 2) * log(2 * pi * sigma02) - (d / 2) * log(1 + sigma02 * m) - Umin
    
    for (k in 1:q) {
      sigma2.i = parameters.tab.adjusted[k * nb.jumps, 2]
      if (k < q) {
        sigma2.iplus1 = parameters.tab.adjusted[(k + 1) * nb.jumps, 2]
      } else {
        sigma2.iplus1 = Inf
      }
      gam.i = parameters.tab.adjusted[k * nb.jumps, 3]
      
      N.i = 10^4
      n.i = 10^6
      
      a.i = 0.5 * (sigma2.i ^ (-1) - sigma2.iplus1 ^ (-1))
      
      piHat.i = 0
      theta = rep(0, d)
      
      for (j in 1:N.i) {
        gradient.U <- gradU(theta + argmin) + theta / sigma2.i
        theta.prop = theta - gam.i * gradient.U + sqrt(2 * gam.i) * rnorm(d)
        gradient.U.prop <-
          gradU(theta.prop + argmin) + theta.prop / sigma2.i
        ratio <-
          -U(theta.prop + argmin) - (2 * sigma2.i) ^ (-1) * sum(theta.prop ^ 2) + U(theta + argmin) + (2 * sigma2.i) ^ (-1) * sum(theta ^ 2) + (4 * gam.i) ^ (-1) * (sum((theta.prop - theta + gam.i * gradient.U) ^ 2) - sum((theta - theta.prop + gam.i * gradient.U.prop) ^ 2))
        if (log(runif(1)) < ratio) {
          theta = theta.prop
        }
        
      }
      
      for (j in (1 + N.i):(n.i + N.i)) {
        gradient.U <- gradU(theta + argmin) + theta / sigma2.i
        theta.prop = theta - gam.i * gradient.U + sqrt(2 * gam.i) * rnorm(d)
        gradient.U.prop <-
          gradU(theta.prop + argmin) + theta.prop / sigma2.i
        ratio <-
          -U(theta.prop + argmin) - (2 * sigma2.i) ^ (-1) * sum(theta.prop ^ 2) + U(theta + argmin) + (2 * sigma2.i) ^ (-1) * sum(theta ^ 2) + (4 * gam.i) ^ (-1) * (sum((theta.prop - theta + gam.i * gradient.U) ^ 2) - sum((theta - theta.prop + gam.i * gradient.U.prop) ^ 2))
        if (log(runif(1)) < ratio) {
          theta = theta.prop
        }
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
  "save-mix-mala-",
  choice,
  ".RData",
  sep = ""
))
