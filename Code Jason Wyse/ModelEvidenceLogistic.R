#written by Jason Wyse. Last modified 28/11/11.

#functions with parameters specified and the data to compute log of posterior etc. 
#	for Logistic Regression models, note that Its should always be greater than BurnIn.


	compare.log.evidences = function(evidence.obj,Its = 25000, BurnIn = 5000,power.posterior.Its = 25000,power.posterior.BurnIn = 5000){
		
		options(warn = -1,digits=5)
		
		NumTests = 7
		
		result = numeric(NumTests)
		times = numeric(NumTests)
		
		t1 = Sys.time()
		result[1] = evidence.obj$log.laplace.evidence()
		t2 = Sys.time()

		cat("\nDone Laplace")

		times[1] = as.numeric(as.difftime(t2-t1))
		t1 = Sys.time()
		result[2] = evidence.obj$log.laplace.map.evidence(Its=250000,BurnIn=50000)
		t2 = Sys.time()
		times[2] = as.numeric(as.difftime(t2-t1))
		
		cat("\nDone Laplace at MAP")

		t1 = Sys.time()
		result[3] = evidence.obj$log.harmonic.mean.estimator.evidence(Its=250000,BurnIn=50000)
		t2 = Sys.time()

		cat("\nDone Harmonic mean")


		times[3] = as.numeric(as.difftime(t2-t1))
		t1 = Sys.time()
		result[4] = evidence.obj$log.chib.jeliazkov.evidence(Its,BurnIn)
		t2 = Sys.time()
		times[4] = as.numeric(as.difftime(t2-t1))

		cat("\nDone Chib and Jeliazkov")

		
		t1 = Sys.time()
		result[5] = evidence.obj$log.annealed.importance.sampling.evidence(sample.size = 2000)
		t2 = Sys.time()
		times[5] = as.numeric(as.difftime(t2-t1))

		cat("\nDone Annealed importance sampling")

		
		t1 = Sys.time()
		result[6] = evidence.obj$log.power.posterior.evidence(Its = power.posterior.Its,BurnIn = power.posterior.BurnIn)
		t2 = Sys.time()
		times[6] = as.numeric(as.difftime(t2-t1))

		cat("\nDone power posterior")

		
		# t1 = Sys.time()
		# result[7] = evidence.obj$log.nested.sampling.evidence(N=2000,Nstep=25000)
		# t2 = Sys.time()
		# times[7] = as.numeric(as.difftime(t2-t1))
		# 
		# 
		# cat("\nDone Nested samlping")
		
		
		R = matrix(nrow = NumTests+1, ncol = 3)
		R[1,] = c("Method","Log evidence","Relative speed")
		colnames(R) = rep(" ",3)
		names = c("Laplace (Newton optimization)","Laplace-Metropolis","Harmonic mean","Chib & Jeliazkov method","Annealed Importance Sampling","Power posteriors","Nested sampling")
		R[2:nrow(R),1] = names
		R[2:nrow(R),2] = as.numeric(result)
		for(j in 1:nrow(R)){
			
			if(is.numeric(R[j,2])){
				R[j,2] = round(R[j,2],digits=2)
				}
			
			}
		
		if(result[4] == "-") R = R[-5,]
		rownames(R) = rep(" ",nrow(R))

		h = sort(times,index.return = TRUE)
		times = times/min(times)
		
		Q = R
		Q = Q[,-3]
		
		T = matrix(nrow = NumTests,ncol=2)
		T[,1] = names[h$ix]
		T[,2] = round(times[h$ix],digits = 2)
		
		R[2:nrow(R),3] = T[,2]
		R[2:nrow(R),1] = R[h$ix+1,1]
		R[2:nrow(R),2] = R[h$ix+1,2]
			
		#print a summary of the results to screen
		
		if(Its == 5500 && BurnIn == 500){
			
			str = "(default).\n\n"
			
			}else{
				
				str =  ". \n\n"
				
				}
				
	#	cat("\nSummary of different approaches for computing the log evidence of a Gaussian linear regression model.\n")
		
			R = as.table(R)	
			
			print.table(R,digits = 2,justify = "right",width = 1)
			
			Result = list()
			Result$method.name = Q[2:nrow(R),1]
			Result$method.log.evidence = Q[2:nrow(R),2]
			
			options(warn = 0,digits = 7)
			
			return(Result)
		
		}



#### HANDLES TO FUNCTIONS TO COMPUTE VARIOUS LOG EVIDENCES ####

evidence.obj = function(DataObs,DataCovariates,PriorMeanMu,PriorPrecisionMu){
	
	#args:	DataObs- observed (binary) data (dependent variable)
	#		DataCovariates- design matrix *with* a column of ones (independent variables) if constant to be included
	#		PriorMeanMu - the mean of Mu in the prior (in this case the regression parameters)
	#		PriorPrecisionMu - the prior precision matrix for Mu (the regression parameters) e.g. tau*I
	
	#Throughout this function theta is understood to be the vector beta (the regression parameters)
	
	#functions in the return:
	#	log.laplace.evidence: gives the Laplace apprxoimation to the marginal likelihood when the mode
	#		is found by using a Newton method
	#	log.laplace.map.evidence: gives the Laplace approximation to the marginal likelihood when the mode
	#		is taken as the MAP from a Monte Carlo run
	#	log.laplace.metropolis.evidence: evaluation of the Laplace approximation at a given pt
	#	log.laplace.kernel.density.mode.evidence: evaluation of the Laplace approximation by using a kernel 
	#		density estimate on the output of the Gibbs to get the mode of each full conditional
	#	log.harmonic.mean.estimator.evidence: performs a Monte Carlo run and computes the correspoding
	#		harmonic mean estimator of the log evidence
	#	log.chib.jeliazkov.evidence: computes the estimator of the log evidence due to Chib and Jeliazkov 		log.annealed.importance.sampling.evidence: computes the estimator of the evidence by annealed 
	#		importance sampling (Neal 2001)
	#	log.nested.sampling.evidence: computes the estimator of the log evidence by Nested Sampling (Skilling 2008)
	#	log.power.posterior.evidence: computes the estimator of the log evidence by Power posteriors (Friel and Pettit 2008)
	#	metropolis.hastings: does a vanilla Metropolis Hastings sampler of the logistic regression model for given number of
	#		iterations and burn-in
	#	reversible.jump: performs a reversible jump algorithm between two models differing by only one variable & subsequent three
	#		functions with suffix .rj compute the log posterior, likelihood and prior for this algorithms
	#	newton.method: do a Newton method to optimize the posterior
	#	log.posterior: compute log posterior at a pt
	#	log.likelihood: compute the log likelihood at a pt
	#	log.prior: compute the prior at a pt
	#	log.ais.posterior: (not used)
	#	log.posterior.gradient: copute the gradient of the log posterior at a pt
	#	hessian: compute the Hessian of the log posterior at a pt
	
		
	#redeifne variables for handiness
	
 	n = length(DataObs)
	y = DataObs
	z = rep(1,n) - DataObs
	X = as.matrix(DataCovariates)
	d = ncol(X)
	mu0 = PriorMeanMu
	tau0 = PriorPrecisionMu
	choltau0 = chol(tau0)
	log.dettau0 = 2*sum(log(diag(choltau0)))
	
	logPi = log(pi)
	log2Pi = log(2*pi)
	
	S = array(dim = c(d,d,n))
	for(i in 1:d){
		for(j in 1:d){
				S[i,j,] = X[,i]*X[,j]
			}
		}
	
	evidence.obj = list()
	
	
	#### LAPLACE APPROXIMATION USING A NEWTON ALGORITHM TO FIND MODE ####
	
	evidence.obj$log.laplace.evidence = function(initial = rep(10^(-8),d)){
				
		theta = evidence.obj$newton.method(initial =  initial)
		
		theta = evidence.obj$newton.method(initial = theta)
		
		#print(theta)
			
		G = -evidence.obj$hessian(theta)
		
		log.detG = log(det(G))
		
		log.laplace.evidence = 0.5*(d)*log2Pi - 0.5*log.detG + evidence.obj$log.posterior(theta)
		
		return(log.laplace.evidence)
		
		}

	#### LAPLACE AT MAP OF METROPOLIS RUN ####
	
	evidence.obj$log.laplace.map.evidence = function(Its = 5500,BurnIn = 500){

		fix = list()
		fix$vars = rep(FALSE,d)
		fix$values = numeric(d)
		
		T = evidence.obj$metropolis.hastings(Its = Its,BurnIn = BurnIn,fix=fix,return.log.posterior=TRUE)
		
		logLikelihood = T[,d+1]
		
		j.star = which(logLikelihood == max(logLikelihood))
		
		if(length(j.star)>1) j.star = j.star[1]
		
		theta.jstar = T[j.star,1:d]
		
		logEvidence = evidence.obj$log.laplace.metropolis.evidence(initial = theta.jstar)
		
		return(logEvidence)
		
		}

	evidence.obj$log.laplace.metropolis.evidence = function(initial = rep(10^(-8))){
		
		theta = initial
			
		G = -evidence.obj$hessian(theta)
		
		log.detG = log(det(G))
		
		log.laplace.evidence = 0.5*(d)*log2Pi - 0.5*log.detG + evidence.obj$log.posterior(theta)
		
		return(log.laplace.evidence)
		
		}
		
	#### LAPLACE AT KERNEL DENSITY MODE OF METROPOLIS RUN ####
	
	evidence.obj$log.laplace.kernel.density.mode.evidence = function(Its = 5500,BurnIn = 500){

		fix = list()
		fix$vars = rep(FALSE,d)
		fix$values = numeric(d)
		
		T = evidence.obj$metropolis.hastings(Its = Its,BurnIn = BurnIn,fix=fix)
				
		theta = numeric(d)
		
		for(j in 1:ncol(T)){
			
			J = density(T[,j])
			index.max.pt = which(J$y == max(J$y))
			theta[j] = J$x[index.max.pt]

			
			}

		logEvidence = evidence.obj$log.laplace.metropolis.evidence(initial = theta)
		
		return(logEvidence)
		
		}


	#### HARMONIC MEAN ESTIMATOR ####

	evidence.obj$log.harmonic.mean.estimator.evidence = function(Its = 5500,BurnIn = 500){
		
		fix = list()
		fix$vars = rep(FALSE,d)
		fix$values = numeric(d)
		
		T = evidence.obj$metropolis.hastings(Its = Its,BurnIn = BurnIn,fix=fix,return.log.likelihood=TRUE)
		
		logLikelihood = T[,d+1]
		
		#try to get stability
		
		a = max(logLikelihood)
		
		harmonic.mean.estimator = log(n) + a - log(sum(exp(-logLikelihood + a)))
		
		return(harmonic.mean.estimator)
		
		}
		
	#### FUNCTION TO PERFORM CHIB & JELIAZKOV METHOD ####
	
	evidence.obj$log.chib.jeliazkov.evidence = function(Its = 5500,BurnIn = 500,sd.prop=0.5){
		
		log.modal.densities = numeric(d)
		
		fix = list()
		fix$vars = rep(FALSE,d)
		fix$values = numeric(d)
		
		#do full Metropolis-Hastings run to begin
		T0 = evidence.obj$metropolis.hastings(Its = Its,BurnIn = BurnIn,fix = fix,sd.prop = sd.prop,quiet = TRUE,return.log.posterior = TRUE)
		
		for(j in 1:(d)){

		#use kernel density estimate to locate mode (this saves alot of time!)
		J = density(T0[,j])
		index.max.pt = which(J$y == max(J$y))
		max.pt = J$x[index.max.pt]
		
		init = c(T0[Its-BurnIn,])
		init = init[-(d+1)]
		#print(init)
		init[j] = max.pt
		
		fix$vars[j] = TRUE
		fix$values[j] = max.pt
		
		#cycle through other variables and sample keeping j^th fixed
		Q0 = evidence.obj$metropolis.hastings(Its = Its-BurnIn,BurnIn = 0,fix = fix,initial = init,quiet = TRUE,return.log.posterior = TRUE)

		#sample the new values from the proposal (rw in this case so ok to do separately)
		q = rnorm(Its-BurnIn,mean=fix$values[j],sd = sd.prop)
		
		#numerator 
		
		T1 = T0
		T1 = T1[,-(d+1)]
		T1[,j] = rep(fix$values[j],nrow(T1))
		
		
		tt = dnorm(T0[,j],mean=fix$values[j],sd = sd.prop,log=TRUE)
		
		r2 = apply(t(T1),2,function(x) evidence.obj$log.posterior(x))
		r1 = T0[,(d+1)] 
		
		logAcc = r2 - r1
				
		logAcc[logAcc > 0] = 0
		
		logAcc = logAcc + tt
		
		numerator = mean(exp(logAcc))
		
		#denominator 
		
		
		Q1 = Q0
		Q1 = Q1[,-(d+1)]
		Q1[,j] = q
	
		
		r1 = Q0[,(d+1)]
		r2 = apply(t(Q1),2,function(x) evidence.obj$log.posterior(x))
			
			
		logAcc = r2 - r1
		logAcc[logAcc > 0] = 0
		
		denominator = mean(exp(logAcc))
		
		log.modal.densities[j] = log(numerator) - log(denominator)
				
		T0 = Q0
		
		}
		
		logEvidence = evidence.obj$log.posterior(fix$values) - sum(log.modal.densities)
				
		return(logEvidence)
		
		}
	
	#### ANNEALED IMPORTANCE SAMPLING ####
	
	evidence.obj$log.annealed.importance.sampling.evidence = function(sample.size = 2000,sd.prop = 0.5,temps = seq(1,0,by=-0.01)^5){
		
		fix = list()
		fix$vars = rep(FALSE,d)
		fix$values = numeric(d)
		
		n.temps = length(temps)
		
		#here the index 1 corresponds to tempering by 1.
		
		N = n.temps
		
		a.tilde = log(tau0[1,1]/(1/sd.prop^2))/log(temps[N-1])
				
		weights = numeric(sample.size)
		
		sd.temps = sd.prop/sqrt((temps^a.tilde))
		
		weights = numeric(sample.size)
		
		for(r in 1:sample.size){
			
			log.weight = 0
			
			#sample from the prior
			
			ep = backsolve(choltau0,rnorm(d),k=d)
					
			theta = mu0 + ep
			
			if(N<4) stop("Terminated Annealed Importance Sampling. Reason: Not a large enough temperature scale")
				
			for(t in (N-1):1){
				
				currentLL = evidence.obj$log.likelihood(theta)
				
				currentAIS.LJL = temps[t]*currentLL + evidence.obj$log.prior(theta)
	
				log.weight = log.weight + (temps[t] - temps[t+1])*currentLL
				
				if(t>1){
				
				for(reps in 1:3){
					
					for(j in 1:d){
						
						theta.prop = theta
						
						theta.prop[j] = rnorm(1,theta[j],sd = sd.temps[t])
						
						proposedLL = evidence.obj$log.likelihood(theta.prop)
						
						proposedAIS.LJL = temps[t]*proposedLL + evidence.obj$log.prior(theta.prop)

						
						logRatio = proposedAIS.LJL - currentAIS.LJL
						
						if(logRatio > 0 || log(runif(1))<logRatio){
									
							theta = theta.prop
							
							currentAIS.LJL = proposedAIS.LJL
							
							currentLL = proposedLL
							
							}
						
						}
					
					}
					
					}
				
				
				}
				
				weights[r] = log.weight
			
			}
			
			max.weight = max(weights)
			
			weights = weights - max.weight
			
			evidence = exp(max.weight)*mean(exp(weights))
			
			return(log(evidence))
		}
	

	#### NESTED SAMPLING ####
	
	evidence.obj$log.nested.sampling.evidence = function(N = 1000,Nstep = 1000){
		
		#an R imitation of the pseudocode in the appendix of Skilling (2006).
				
		logw = log(1.0 - exp(-1.0/N))
		
		Theta = matrix(nrow = d,ncol = N)
		
		log.likelihoods = numeric(N)
		
		#do initial draws from the prior
		
		for(j in 1:N){
			
			Theta[,j] = mu0 + backsolve(choltau0,rnorm(d),k=d)
						
			}
		
		log.likelihoods = apply(Theta,2,function(x) evidence.obj$log.likelihood(x))
		
		#set logZ to the equivalent of log(0)... largest negative number possible.
		
		logZ = -.Machine$double.xmax
		
		log.contribution = 0
						
		nest = 1
		
		breakwhile = FALSE
		
		#same termination criterion as Chopin and Robert (2010)
		
		while((breakwhile==FALSE)){
						
			r = which(log.likelihoods == min(log.likelihoods))
			
			l.truncate = log.likelihoods[r]
			
			#print(l.truncate)
			
			log.contribution = logw + l.truncate
			
			if(log.contribution < -8*log(10) + logZ){
				
				breakwhile = TRUE
				
				}
			
			logZnew = evidence.obj$PLUS(logZ,log.contribution)
						
			#sample one of the qualifying likelihoods to start Metropolis steps
			
			s = sample(setdiff(1:N,r), size = 1)
			
			#starting value
			
			Theta[,r] = Theta[,s]
						
			Theta[,r] = evidence.obj$evolve(Theta[,r],l.truncate)
			
			log.likelihoods[r] = evidence.obj$log.likelihood(Theta[,r])
			
			#shrink interval
				
			logw = logw - 1/N
			
			nest = nest+1
			
			logZ = logZnew
			
			}
		
		#cat("\nBefore correction logZ = ",logZ)

		logw = -nest/N - log(N)
		
		for(i in 1:N){
			
			logZnew = evidence.obj$PLUS(logZ,logw+log.likelihoods[i])
			
			logZ = logZnew
			
			}
			
		return(logZ)	
		
		}
		

	evidence.obj$evolve = function(theta,l){
				
		currentLogPosterior = evidence.obj$log.prior(theta)
		
		acceptances = 0
		
		while(acceptances<100){
		
			for(j in 1:d){
								
				theta.prop = theta
				
				theta.prop[j] = rnorm(1,mean = theta[j],sd = 0.1)
				
				proposedLogPosterior = evidence.obj$log.prior(theta.prop)
				
				logRatio = proposedLogPosterior - currentLogPosterior
				
				lproposed = evidence.obj$log.likelihood(theta.prop)
				
				if((logRatio > 0 || log(runif(1))<logRatio) && lproposed > l){
										
					theta = theta.prop
					
					currentLogPosterior = proposedLogPosterior
					
					acceptances = acceptances+1
					
					}
								
				}		
		}
		
		return(theta)
		
		}


	#### POWER POSTERIORS ####
	
	evidence.obj$log.power.posterior.evidence = function(Its = 25000,BurnIn = 5000,temps = seq(0,1,by=0.01)^5,sd.prop = 0.5){
		# previously, temps = seq(0,1,by=0.1)^5
	  
		N = length(temps)
		
		a.tilde = log(tau0[1,1]/(1/sd.prop^2))/log(temps[2])
		
		#print(a.tilde)		
		sd.temps = sd.prop/sqrt((temps^a.tilde))
		
		if(N>1) sd.temps[1] = sd.temps[2]
		
		#cycle through temperatures estimating the expected deviance each time
		
		expected.deviance = numeric(N)
		
		expected.theta = numeric(d)
		
		#draw from prior to initialize first chain
		
		theta = mu0 + backsolve(choltau0,rnorm(d),k=d)
				
		acceptances = list()
		acceptances$proposed = numeric(N)
		acceptances$accepted = numeric(N)
		
		for(t in 1:N){
			
			expected.theta = numeric(d)

			currentLL = evidence.obj$log.likelihood(theta)

			currentLJL = temps[t]*currentLL + evidence.obj$log.prior(theta)
						
			for(It in 1:Its){
				
				for(j in 1:d){
						
						acceptances$proposed[t] = acceptances$proposed[t]+1
						
						theta.prop = theta
						
						theta.prop[j] = rnorm(1,mean = theta[j],sd = sd.temps[t])
												
						proposedLL = evidence.obj$log.likelihood(theta.prop)
						
						proposedLJL = temps[t]*proposedLL + evidence.obj$log.prior(theta.prop)
						
						logRatio = 	proposedLJL - currentLJL											
						if(logRatio > 0 || log(runif(1))<logRatio){
										
							acceptances$accepted[t] = acceptances$accepted[t]+1
							
							theta = theta.prop
							
							currentLJL = proposedLJL
							
							currentLL = proposedLL
							
							}
						}
			
				if(It > BurnIn){
					
					expected.deviance[t] = expected.deviance[t] + currentLL
					expected.theta = ((It-BurnIn-1)*expected.theta + theta)/(It-BurnIn)
					
					}
				
				
			
			}
		
		#initialize new chang with mean value of old
		
		theta = expected.theta

		}
		
		expected.deviance = expected.deviance/(Its - BurnIn)
				
		integration.weights = numeric(N-1)
		function.vals = numeric(N-1)
		
		for(i in 1:(N-1)){
			
			integration.weights[i] = temps[i+1] - temps[i]
			
			function.vals[i] = 0.5*(expected.deviance[i+1]+expected.deviance[i])			
			}
		
		log.evidence = t(integration.weights)%*%function.vals
		
		return(log.evidence)
		
		}
	
	
		
	#### FUNCTION TO DO A VANILLA METROPOLIS HASTINGS FOR LOGISTIC REGRESSION MODEL ####

	evidence.obj$metropolis.hastings = function(Its = 11000, BurnIn = 1000, initial = rep(10^(-6),d),fix,sd.prop = 0.5,quiet = TRUE,return.log.posterior=FALSE,return.log.likelihood=FALSE){
		
		if(length(which(fix$vars == FALSE)) == 0){
						
			tt = rep(initial,Its-BurnIn)
			
			T = matrix(tt,ncol=d,byrow = TRUE)
			
			h = evidence.obj$log.posterior(T[1,])
			
			log.posterior.values = rep(h,nrow(T))
			
			l = evidence.obj$log.likelihood(T[1,])
			
			log.likelihood.values = rep(l,nrow(T))
			
			if(return.log.posterior == TRUE){
			
			T = cbind(T,log.posterior.values)
			
			}
		
		if(return.log.likelihood == TRUE){
			
			T = cbind(T,log.likelihood.values)
			
			}		
			
						
			return(T)
			
			}
		
		acceptances = list()
		acceptances$proposed = numeric(d)
		acceptances$accepted = numeric(d)
		
		T = matrix(nrow = Its-BurnIn,ncol = d)
		
		log.posterior.values = numeric(Its-BurnIn)
		
		log.likelihood.values = numeric(Its-BurnIn)
		
		theta = initial
		
		currentLogPosterior = evidence.obj$log.posterior(theta)
		
		sample.vars = which(fix$vars == FALSE)
				
		fix.vars = which(fix$vars == TRUE)
		if(length(fix.vars)>0) theta[fix.vars] = fix$values[fix.vars]
		
		for(ItNum in 1:Its){
			
			#visit each parameter in turn
			
			for(j in sample.vars){
				
				acceptances$proposed[j] = acceptances$proposed[j] + 1
				
				theta.prop = theta
				
				theta.prop[j] = rnorm(1,mean = theta[j],sd = sd.prop)
				
				proposedLogPosterior = evidence.obj$log.posterior(theta.prop)
				
				logRatio = proposedLogPosterior - currentLogPosterior
				
				if(logRatio > 0 || log(runif(1))<logRatio){
					
					acceptances$accepted[j] = acceptances$accepted[j]+1
					
					theta = theta.prop
					
					currentLogPosterior = proposedLogPosterior
					
					}
				
				if(ItNum > BurnIn){
					
					T[ItNum-BurnIn,] = theta
					
					log.posterior.values[ItNum-BurnIn] = currentLogPosterior
					
					log.likelihood.values[ItNum-BurnIn] = currentLogPosterior+evidence.obj$log.prior(theta)
					
					}
				
				}
			
			}
		
		if(!quiet){
			cat("\nAcceptance rates for each of the regression parameters:\n")
			print(acceptances$accepted/acceptances$proposed)
		}
		
		
		if(return.log.posterior == TRUE){
			
			T = cbind(T,log.posterior.values)
			
			}
		
		if(return.log.likelihood == TRUE){
			
			T = cbind(T,log.likelihood.values)
			
			}
		
		return(T)
		
		}	
		
	#### FUNCTION TO DO A REVERSIBLE JUMP ALGORITHM TO EVALUATE BF BETWEEN TWO MODELS ####
		
	evidence.obj$reversible.jump = function(Its = 25000,BurnIn = 5000,model1,model2,prior1,prior2,sd.prop = 0.5){
		
		# model 1 & 2 will be a set of variable numbers corresponding to cols in X (including 1 if constant to be included) 
		# prior 1 & 2 will give the prior probability assigned to the models usually using a corrected RJ scheme to allow for frequent jumps between models
		#WARNING: model 1 should be the one with smaller dimension
		
		d1 = length(model1)
		d2 = length(model2)
		
		if(abs(d1-d2)>1) stop("\nReversible jump algorithm only currently able to handle models with a difference in dimension of 1 only\n")
		
		#start off in model 1 with draw from prior
		store.model = numeric(Its-BurnIn)
		
		acceptances = list()
		acceptances$j1to2 = 0
		acceptances$p1to2 = 0
		acceptances$j2to1 = 0
		acceptances$p2to1 = 0
		
		model = model1
		M = 1
		D = length(model)
		
		theta = mu0[1:D] + backsolve(choltau0[1:D,1:D],rnorm(D),k = D)
		
		print(theta)
		
		currentLJL = evidence.obj$log.posterior.rj(theta,model)
		
		for(It in 1:Its){
			
			for(j in 1:length(model)){
				
				theta.prop = theta
				
				theta.prop[j] = theta[j] + sd.prop*rnorm(1)
				
				proposedLJL = evidence.obj$log.posterior.rj(theta.prop,model)
				
				logRatio = proposedLJL - currentLJL
				
				if(logRatio > 0 || log(runif(1))< logRatio){
					
					theta = theta.prop
					
					currentLJL = proposedLJL
					
					}
				}
					
					
				#propose a swap to the other model
				
				if(M == 1){
					
					acceptances$p1to2 = acceptances$p1to2 + 1
					
					model.prop = model2
					
					r = rnorm(1,mean = 0,sd = sd.prop)
					
					theta.prop = c(theta,r)
					
					proposedLJL = evidence.obj$log.posterior.rj(theta.prop,model2)
					
					logRatio = proposedLJL - currentLJL - dnorm(r,0,sd.prop,log=TRUE) + log(prior2) - log(prior1)
					
					if(logRatio > 0 || log(runif(1))<logRatio){
						
						acceptances$j1to2 = acceptances$j1to2 + 1
						
						theta = theta.prop
						
						model = model2
						
						M = 2
						
						currentLJL = proposedLJL
						
						}
					
					}else{
						
						acceptances$p2to1 = acceptances$p2to1 + 1
						
						model.prop = model1
						
						r = theta[d2]
						
						theta.prop = theta[-d2]
						
						proposedLJL = evidence.obj$log.posterior.rj(theta.prop,model1)
						
						logRatio = proposedLJL - currentLJL + dnorm(r,0,sd.prop,log = TRUE) + log(prior1) - log(prior2)
						
						if(logRatio > 0 || log(runif(1))<logRatio){
							
							acceptances$j2to1 = acceptances$j2to1 + 1
							
							theta = theta.prop
							
							model = model1
							
							M = 1
							
							currentLJL = proposedLJL
							
							}
						
						
						}

				
				
				
				if(It > BurnIn){
					
					store.model[It-BurnIn] = M
					
					}
				
			

			}
		
		
		cat("\nSummary of Run:\n Acceptance rates from model 1 to 2 = ",acceptances$j1to2/acceptances$p1to2,"\n Acceptance rates from model 2 to 1 = ",acceptances$j2to1/acceptances$p2to1)
		
		bf = (prior2/prior1)*length(which(store.model == 1))/(length(which(store.model == 2)))
		
		return(bf)
		
		
		}
	
	evidence.obj$log.posterior.rj = function(theta,model){
		
		logPosterior = evidence.obj$log.likelihood.rj(theta,model) + evidence.obj$log.prior.rj(theta)
		
		}
	
	evidence.obj$log.likelihood.rj = function(theta,model){
			
			l = X[,model]%*%theta 
						
			g = exp(l)
			
			logLikelihood = t(y)%*%l - sum(log(1+g))
			
			return(logLikelihood)
			
		}
		
	evidence.obj$log.prior.rj = function(theta){
		
		d0 = length(theta)
		
		logPrior = -0.5*d0*log2Pi + 0.5*d0*log(tau0[1,1]) - 0.5*t(theta - mu0[1:d0])%*%tau0[1:d0,1:d0]%*%(theta-mu0[1:d0])
		
		}	
	
	
	#### FUNCTION TO CARRY OUT NEWTON METHOD ####
	
	evidence.obj$newton.method = function(tolerance = 0.01,maxsteps = 200,initial = rep(10^(-3),d)){
		
		#use the MLE as a starting value
		
		theta = initial
				
		numsteps = 0
		
		tolcriterion = 100
		
		step.size = 1
		
		while(step.size > 0){
			
		step.size = step.size/2
				
		while(tolcriterion>tolerance && numsteps < maxsteps){
			
			G = evidence.obj$hessian(theta)
						
			invG = solve(G)
			
			thetaNew = theta - step.size*invG%*%evidence.obj$log.posterior.gradient(theta)
						
			if(evidence.obj$log.posterior(thetaNew) < evidence.obj$log.posterior(theta)){
				break
				}
			
			tolcriterion = abs(evidence.obj$log.posterior(thetaNew)-evidence.obj$log.posterior(theta))
			
			theta = thetaNew
			
			numsteps = numsteps+1
			
			}
			
			}
			
			if(numsteps == maxsteps) warning('Maximum number of steps reached in Newton method\n\tFailed to reach tolerance level\n\tThis may affect accuracy of results')
			#print(theta)
			return(theta)
		
		}	
	

	evidence.obj$log.posterior = function(theta){
		
		logPosterior = evidence.obj$log.likelihood(theta) + evidence.obj$log.prior(theta)
		
		return(logPosterior)
		
		}
		
	evidence.obj$log.likelihood = function(theta){
			
			l = X%*%theta
			
			g = exp(l)
			
			logLikelihood = t(y)%*%l - sum(log(1+g))
			
			return(logLikelihood)
			
		}
		
	evidence.obj$log.prior = function(theta){
		
		logPrior = -0.5*d*log2Pi + 0.5*log.dettau0 - 0.5*t(theta - mu0)%*%tau0%*%(theta-mu0)
		
		}
	
		#### FUNCTION FOR EVALUATING TEMPERED LOG PRIOR POSTERIOR MIXTURE ####
	
	evidence.obj$log.ais.posterior = function(theta,t){
		
		logAISPosterior = t*evidence.obj$log.likelihood(theta) + evidence.obj$log.prior(theta)		
		return(logAISPosterior)
		
		}
	
	evidence.obj$log.posterior.gradient = function(theta){
		
			l = X%*%theta
			
			g = exp(l)
			
			mu = g/(1+g)
			
			logPosteriorGradient = t(X)%*%(y-mu) - tau0%*%(theta-mu0)
			
			return(logPosteriorGradient)	
		
		}
	
	evidence.obj$hessian = function(theta){
		
		H = matrix(nrow=d,ncol=d)
		
		l = X%*%theta
			
		g = exp(l)
			
		mu = g/(1+g)
		
		zmu = 1 - mu
		
		mu.zmu = mu*zmu
		
		for(i in 1:d){
			
			for(j in 1:d){
				
				H[i,j] = -t(S[i,j,])%*%mu.zmu
				
				}
			
			}
		
		#prior
		
		H = H - tau0
		
		return(H)
		
		}
		

	return(evidence.obj)
	
	}
	
	

		
