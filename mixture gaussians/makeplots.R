# Mixture of gaussians

load("save-mix-2.Rdata")
M1 <- lnZ.tab
load("save-mix-mala-2.Rdata")
lnZ.mala <- mean(lnZ.tab)

boxplot(M1)
points(x = 1, y=lnZ.mala, col = 'red', cex = 1.2, pch = 21, bg = 'red')
