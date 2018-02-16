#### RUN EXPERIMENTS FOR THE DIFFERENT METHODS ####

source("ModelEvidence.R")

run.experiment = function(evidence.obj,n.methods,av.times){

	A = numeric(n.methods)

	n.its = 505000
	burn.in = 101000
	
	t1 = Sys.time()
	A[1] = evidence.obj$log.laplace.evidence.map(Its = n.its,BurnIn = burn.in)
	t2 = Sys.time()
	av.times[1] = av.times[1]+as.numeric(t2-t1)
	t1 = Sys.time()
	A[2] = evidence.obj$log.harmonic.mean.estimator.evidence(Its = n.its,BurnIn = burn.in)
	t2 = Sys.time()
	av.times[2] = av.times[2]+as.numeric(t2-t1)
	t1 = Sys.time()
	A[3] = evidence.obj$log.chib.evidence(Its = 205000,BurnIn = 55000)
	t2 = Sys.time()
	av.times[3] = av.times[3]+as.numeric(t2-t1)
	t1 = Sys.time()
	A[4] = evidence.obj$log.annealed.importance.sampling.evidence(sample.size = 1000,temps = seq(1,0,by=-0.01)^5)
	t2 = Sys.time()
	av.times[4] = av.times[4]+as.numeric(t2-t1)
	t1 = Sys.time()
	A[5] = evidence.obj$log.nested.sampling.evidence(N = 1000,tolerance = 1/sqrt(505000))
	t2 = Sys.time()
	av.times[5] = av.times[5]+as.numeric(t2-t1)
	t1 = Sys.time()
	A[6] = evidence.obj$log.power.posterior.evidence(Its = 5000,BurnIn = 1000,temps = seq(0,1,by=0.01)^5)
	t2 = Sys.time()
	av.times[6] = av.times[6]+as.numeric(t2-t1)

	return(A)

}



#### READ DATA AND LOAD MODELS ####

RadiataPine = read.table("RadiataPine.txt",sep = " ",header = TRUE)

y = RadiataPine$y

n = length(y)

X1 = cbind(rep(1,n),RadiataPine$x-mean(RadiataPine$x))

X2 = cbind(rep(1,n),RadiataPine$z - mean(RadiataPine$z))

mu0 = c(3000,185)
Lambda0 = diag(1,2)
Lambda0[1,1] = 0.06
Lambda0[2,2] = 6.00

E1 = evidence.obj(y,X1,mu0,Lambda0,3,2*300^2)
E2 = evidence.obj(y,X2,mu0,Lambda0,3,2*300^2)

#### CARRY OUT EXPERIMENTS ####
#
n.times = 18

n.methods = 6

store.log.evidences = matrix(nrow = n.times,ncol = n.methods)

av.times = numeric(n.methods)

for(i in 1:n.times){

	J = run.experiment(E1,n.methods,av.times)

	store.log.evidences[i,] = J
	
	cat("\nDone to experiment ",i,"\n")
}

write.table(store.log.evidences,"RadiataModel1logevidences.txt",sep=",",row.names=FALSE,col.names=FALSE,eol="\n")





n.times = 18

n.methods = 6

store.log.evidences = matrix(nrow = n.times,ncol = n.methods)

av.times = numeric(n.methods)

for(i in 1:n.times){

	J = run.experiment(E2,n.methods,av.times)

	store.log.evidences[i,] = J
	
	cat("\nDone to experiment ",i,"\n")
	
	write.table(store.log.evidences[i,],append=TRUE,"RadiataModel2logevidences.txt",sep=",",row.names=FALSE,col.names=FALSE,eol="\n")

}

#write.table(store.log.evidences,append=TRUE,"RadiataModel2Rerun.txt",sep=",",row.names=FALSE,col.names=FALSE,eol="\n")

#### 	SUMMARY ####

cat("\nMean evidence \n")
apply(store.log.evidences,2,mean)
cat("\nStandard Deviation evidence \n")
apply(store.log.evidences,2,sd)




