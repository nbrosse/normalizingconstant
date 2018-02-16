# Logistic regression on the Pima Indians dataset

chosen.model <- 1 # Model 1 or 2
nb.runs <- 10 # Number of independent runs 
nb.jumps <- 5 # Number of jumped phases to accelerate the algorithm
choice <- 1 # choice=1, independent runs on each phase to measure the variability of each increment 
            # choice=2, independent runs on the entire algorithm to measure the variability of log(normalizing constant)

# Read the data

A = read.table("pima.txt", sep = ",")
y = A[, 1]
X = cbind(rep(1, nrow(A)), A[2:ncol(A)])

model1 = c(1, 2, 3, 6, 7)
model2 = c(1, 2, 3, 6, 7, 8)

X1 = X[, model1]
X1M = as.matrix(X1)

X2 = X[, model2]
X2M = as.matrix(X2)

tau = 0.01

if (chosen.model == 1) {
  X1M = X1M
} else {
  X1M = X2M
}

# Library maxLik to compute the minimum and argmin of U

# library(maxLik)
#
# log.pi = function(beta) {
#   term1 <- as.numeric(t(beta) %*% t(X1M) %*% y)
#   term2 <- sum(log(1 + exp(X1M %*% beta)))
#   term3 <- (tau / 2) * norm(beta, type = "2") ^ 2
#   term4 <- (ncol(X1M) / 2) * log(tau / (2 * pi))
#   return(term1 - term2 - term3 + term4)
# }
#
# gradient.log.pi = function(beta) {
#   term1 <- as.vector(t(X1M) %*% y)
#   term2 <- as.vector(t(1 / (1 + exp(-X1M %*% beta))) %*% X1M)
#   term3 <- tau * beta
#   return(term1 - term2 - term3)
# }
#
# ml <-
#   maxLik(logLik = log.pi,
#          start = rep(0, ncol(X1M)),
#          grad = gradient.log.pi)
#
# ml.save<-ml[1:2]
# save(ml.save, file = paste("maxLik-log-tau001-M", chosen.model, ".RData", sep=""))

load(file = paste("maxLik-log-tau001-M", chosen.model, ".RData", sep = ""))
ml <- ml.save

argmin <- ml$estimate

vp = eigen(t(X1M) %*% X1M, symmetric = T, only.values = T)

L = max(vp$values) / 4 + tau
m = tau
d = ncol(X1M)

epsilon = 0.1
mu = 0.1

sigma2.i <- sigma02 <- (2 * log(1 + epsilon / 3)) / (d * (L - m))

# Count the number of phases

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

lnZ <-
  (d / 2) * log(2 * pi * sigma02) - (d / 2) * log(1 + sigma02 * m) +
  ml$maximum

term1 <- as.vector(t(X1M) %*% y) # useful term for the computation of gradient U

n = length(sigma2.tab)
q = n %/% nb.jumps

if (choice == 1) {
  for (k in 1:q) {
    sigma2.i = sigma2.tab[k * nb.jumps]
    if (k < q) {
      sigma2.iplus1 = sigma2.tab[(k + 1) * nb.jumps]
    } else {
      sigma2.iplus1 = Inf
    }
    m.i = m + 1 / sigma2.i
    L.i = L + 1 / sigma2.i
    
    N.i = 10 ^ 4
    if (k <= 30) {
      gam.i = 10 ^ (-2) * (m.i + L.i) ^ (-1)
      n.i = 10 ^ 6
    }
    else {
      gam.i = 10 ^ (-1) * (m.i + L.i) ^ (-1)
      n.i = 10 ^ 5
    }
    
    a.i = 0.5 * (sigma2.i ^ (-1) - sigma2.iplus1 ^ (-1))
    
    chain <- rep(NA, nb.runs)
    
    for (l in 1:nb.runs) {
      piHat.i = 0
      theta = rep(0, d)
      
      for (j in 1:N.i) {
        X = theta + argmin
        term2 <- as.vector(t(1 / (1 + exp(-X1M %*% X))) %*% X1M)
        term3 <- tau * X
        gradient.U <- (-term1) + term2 + term3
        gradient.U.i = theta / sigma2.i + gradient.U
        theta = theta - gam.i * gradient.U.i + sqrt(2 * gam.i) * rnorm(d)
      }
      
      for (j in (1 + N.i):(n.i + N.i)) {
        X = theta + argmin
        term2 <- as.vector(t(1 / (1 + exp(-X1M %*% X))) %*% X1M)
        term3 <- tau * X
        gradient.U <- (-term1) + term2 + term3
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
  lnZ.tab <- rep(NA, nb.runs)
  
  for (l in 1:nb.runs) {
    lnZ = (d / 2) * log(2 * pi * sigma02) - (d / 2) * log(1 + sigma02 * m) +
      ml$maximum
    
    for (k in 1:q) {
      sigma2.i = sigma2.tab[k * nb.jumps]
      if (k < q) {
        sigma2.iplus1 = sigma2.tab[(k + 1) * nb.jumps]
      } else {
        sigma2.iplus1 = Inf
      }
      m.i = m + 1 / sigma2.i
      L.i = L + 1 / sigma2.i
      
      N.i = 10 ^ 4
      if (k <= 30) {
        gam.i = 10 ^ (-2) * (m.i + L.i) ^ (-1)
        n.i = 10 ^ 6
      }
      else {
        gam.i = 10 ^ (-1) * (m.i + L.i) ^ (-1)
        n.i = 10 ^ 5
      }
      
      a.i = 0.5 * (sigma2.i ^ (-1) - sigma2.iplus1 ^ (-1))
      
      piHat.i = 0
      theta = rep(0, d)
      
      for (j in 1:N.i) {
        X = theta + argmin
        term2 <- as.vector(t(1 / (1 + exp(-X1M %*% X))) %*% X1M)
        term3 <- tau * X
        gradient.U <- (-term1) + term2 + term3
        gradient.U.i = theta / sigma2.i + gradient.U
        theta = theta - gam.i * gradient.U.i + sqrt(2 * gam.i) * rnorm(d)
      }
      
      for (j in (1 + N.i):(n.i + N.i)) {
        X = theta + argmin
        term2 <- as.vector(t(1 / (1 + exp(-X1M %*% X))) %*% X1M)
        term3 <- tau * X
        gradient.U <- (-term1) + term2 + term3
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
  "save-reg-log-model",
  chosen.model,
  "-",
  choice,
  ".RData",
  sep = ""
))
