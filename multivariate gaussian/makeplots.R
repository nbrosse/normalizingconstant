# Boxplots for the multivariate Gaussian in dimension 10, 25, 50

d=10
load(file = "save-gaus-d10-2.RData")
lnZ.tab.10 <- lnZ.tab
lnZ.Exact.10 <- lnZ.Exact

d=25
load(file = "save-gaus-d25-2.RData")
lnZ.tab.25 <- lnZ.tab
lnZ.Exact.25 <- lnZ.Exact

d=50
load(file = "save-gaus-d50-2.RData")
lnZ.tab.50 <- lnZ.tab
lnZ.Exact.50 <- lnZ.Exact

par(mfrow=c(1,3))

boxplot(lnZ.tab.10)
points(x = 1, y=lnZ.Exact.10, col = 'red', cex = 1.2, pch = 21, bg = 'red')
title("Dimension 10")

boxplot(lnZ.tab.25)
points(x = 1, y=lnZ.Exact.25, col = 'red', cex = 1.2, pch = 21, bg = 'red')
title("Dimension 25")

boxplot(lnZ.tab.50)
points(x = 1, y=lnZ.Exact.50, col = 'red', cex = 1.2, pch = 21, bg = 'red')
title("Dimension 50")
