
#written by Jason Wyse: last modified 2/9/2011

#This is just a front programme to do the computations presented in the paper.
#The main lines of interest is the call to evidence.obj(). The arguments to this
#function are the response (y), design matrix (X), a vector of zeros giving the prior
#mean on the regression parameters and a diagonal matrix giving the prior precisions
#of the regression parameters.

#source functions
source("ModelEvidenceLogistic.R")

#read data
A = read.table("pima.txt",sep=",")
y = A[,1]
X = cbind(rep(1,nrow(A)),A[2:ncol(A)])

model1 = c(1,2,3,6,7)
model2 = c(1,2,3,6,7,8)

X1 = X[,model1]

X2 = X[,model2]


E1 = evidence.obj(y,X1,rep(0,ncol(X1)),0.01*diag(1,ncol(X1)))
E2 = evidence.obj(y,X2,rep(0,ncol(X2)),0.01*diag(1,ncol(X2)))

for(j in 1:5){

#E = evidence.obj(y,X,rep(0,ncol(X)),0.01*diag(1,ncol(X)))
cat("\n****Running Model 1****\n")
A = compare.log.evidences(E1)

names = A$method.name
Result1 = A$method.log.evidence 

cat("\n****Running Model 2****\n")
B = compare.log.evidences(E2)
Result2 = B$method.log.evidence

cat("\nSummary of Bayes Factors for different approaches for computing the log evidence of a Logistic regression model.")

BayesFactor = matrix(nrow = length(names)+1,ncol = 2)
BayesFactor[2:nrow(BayesFactor),1] = names
BayesFactor[2:nrow(BayesFactor),2] = exp(-as.numeric(Result2) + as.numeric(Result1))
rownames(BayesFactor) = rep(" ",nrow(BayesFactor))
colnames(BayesFactor) = c(" "," ")
BayesFactor[1,] = c("Method","Bayes Factor")
BayesFactor = as.table(BayesFactor)
cat("\n\n\n\n")
print.table(BayesFactor,justify = "right",width = 1)

}