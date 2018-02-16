# Logistic regression

load("save-log-Wyse.RData") 
# data coming from the execution of the code of Jason Wyse
# see the file "code Jason Wyse" and run "RunPimaModified.R"

tab.M1 <- as.numeric(result.M1[[1]])
tab.M2 <- as.numeric(result.M2[[1]])

for(i in 2:10){
  tab.M1 <- rbind(tab.M1, as.numeric(result.M1[[i]]))
  tab.M2 <- rbind(tab.M2, as.numeric(result.M2[[i]]))
}

tab.M1 <- tab.M1[,-ncol(tab.M1)]
tab.M2 <- tab.M2[,-ncol(tab.M2)]

load("save-reg-log-model1-2.Rdata")
M1 <- lnZ.tab
load("save-reg-log-model2-2.Rdata")
M2 <- lnZ.tab

tab.M1 <- cbind(tab.M1, M1)
tab.M2 <- cbind(tab.M2, M2)

tab.M1 <- tab.M1[,-3]
tab.M2 <- tab.M2[,-3]

noms <- c("L", "L-MAP", "C", "AIS", "PP", "AV")

colnames(tab.M2) <- colnames(tab.M1) <- noms

par(mfrow=c(2,1))

boxplot(tab.M1)
title("Model 1")
boxplot(tab.M2)
title("Model 2")

# Gaussian regression

par(mfrow=c(1,1))

load("save-reg-gaus-model1-1.RData")
mean.increment <- analysis.tab[, "mean"]
# sig <- analysis.tab[, "sigma2"]
mean.plus <- mean.increment + analysis.tab[, "standard deviation"]
mean.less <- mean.increment - analysis.tab[, "standard deviation"]
y <- cbind(mean.increment, mean.plus, mean.less)

matplot(y, type = "l", xlab = "phase i", ylab = expression(hat(pi)[i] (g[i])), col = c(1,2,4))

# Logistic regression

par(mfrow=c(1,1))

load("save-reg-log-model1-1.RData")
mean.increment <- analysis.tab[1:49, "mean"]
# sig <- analysis.tab[, "sigma2"]
mean.plus <- mean.increment + analysis.tab[1:49, "standard deviation"]
mean.less <- mean.increment - analysis.tab[1:49, "standard deviation"]
y <- cbind(mean.increment, mean.plus, mean.less)

matplot(y, type = "l", xlab = "phase i", ylab = expression(hat(pi)[i] (g[i])), col = c(1,2,4))

# Gaussian regression

load("save-reg-gaus-model1-2.Rdata")
M1 <- lnZ.tab
lnZ.Exact.1 <- lnZ.Exact
load("save-reg-gaus-model2-2.Rdata")
M2 <- lnZ.tab
lnZ.Exact.2 <- lnZ.Exact

gaus <- cbind(M1, M2)

par(mfrow=c(1,2))

boxplot(M1)
points(x = 1, y=lnZ.Exact.1, col = 'red', cex = 1.2, pch = 21, bg = 'red')
title("Model 1")

boxplot(M2)
points(x = 1, y=lnZ.Exact.2, col = 'red', cex = 1.2, pch = 21, bg = 'red')
title("Model 2")

