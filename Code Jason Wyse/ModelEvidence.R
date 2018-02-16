#written by Jason Wyse. Last modified 29/11/11.

#create handles to functions with parameters specified and the data to compute log of posterior etc.

	#### FUNCTION TO PRODUCE SUMMARY OF EVIDENCES FROM DIFFERENT APPROACHES ####

	compare.log.evidences = function(evidence.obj,Its = 5500, BurnIn = 500,power.posterior.Its = 25000,power.posterior.BurnIn = 5000){
		
		options(warn = -1,digits=2)
		
		NumTests = 9
		
		result = numeric(NumTests)
		times = numeric(NumTests)
		
		t1 = Sys.time()
		result[1] = evidence.obj$log.analytic.evidence()
		t2 = Sys.time()
		#cat("\nAnalytic evidence computation complete. Time taken was ",as.numeric(as.difftime(t2-t1))," seconds.\n")
		times[1] = as.numeric(as.difftime(t2-t1))
		t1 = Sys.time()
		result[2] = evidence.obj$log.laplace.evidence()
		t2 = Sys.time()
		#cat("\nLaplace approximation of evidence complete. Time taken was ",as.numeric(as.difftime(t2-t1))," seconds.\n")
		times[2] = as.numeric(as.difftime(t2-t1))
		t1 = Sys.time()
		result[3] = evidence.obj$log.laplace.evidence.map(Its,BurnIn)
		t2 = Sys.time()
		#cat("\nLaplace-Metropolis approximation of evidence complete. Time taken was ",as.numeric(as.difftime(t2-t1))," seconds.\n")
		times[3] = as.numeric(as.difftime(t2-t1))
		t1 = Sys.time()
		result[4] = evidence.obj$log.laplace.evidence.mle()
		t2 = Sys.time()
		#cat("\nLaplace MLE approximation of evidence complete. Time taken was ",as.numeric(as.difftime(t2-t1))," seconds.\n")
		times[4] = as.numeric(as.difftime(t2-t1))
		t1 = Sys.time()
		result[5] = evidence.obj$log.harmonic.mean.estimator.evidence(Its,BurnIn)
		t2 = Sys.time()
		#cat("\nHarmonic mean approximation of evidence complete. Time taken was ",as.numeric(as.difftime(t2-t1))," seconds.\n")
		times[5] = as.numeric(as.difftime(t2-t1))
		t1 = Sys.time()
		result[6] = evidence.obj$log.chib.evidence(Its,BurnIn)
		t2 = Sys.time()
		#cat("\nChib's approximation of evidence complete. Time taken was ",as.numeric(as.difftime(t2-t1))," seconds.\n")
		times[6] = as.numeric(as.difftime(t2-t1))
		t1 = Sys.time()
		result[7] = evidence.obj$log.annealed.importance.sampling.evidence(sample.size = 10000)
		t2 = Sys.time()
		times[7] = as.numeric(as.difftime(t2-t1))
		t1 = Sys.time()
		result[8] = evidence.obj$log.nested.sampling.evidence(N = Its,Nstep = 10000)
		t2 = Sys.time()
		times[8] = as.numeric(as.difftime(t2-t1))
		t1 = Sys.time()
		result[9] = evidence.obj$log.power.posterior.evidence(Its = power.posterior.Its,BurnIn = power.posterior.BurnIn)
		t2 = Sys.time()
		times[9] = as.numeric(as.difftime(t2-t1))
		
		
		R = matrix(nrow = NumTests+1, ncol = 3)
		R[1,] = c("Method","Log evidence","Relative speed")
		colnames(R) = rep(" ",3)
		names = c("Analytic", "Laplace (Newton optimization)","Laplace-Metropolis","Laplace at MLE","Harmonic mean","Chib's method","Annealed Importance Sampling","Nested Sampling","Power posteriors")
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
				
		
		
		cat("\nSummary of different approaches for computing the log evidence of a Gaussian linear regression model.\n")
		
#		cat("\nThe results dependent on a Monte Carlo experiment were carried out with a Gibbs sampler with ",BurnIn," burn-in iterations and ",Its - BurnIn," subsequent iterations ",str)
#		
			R = as.table(R)	
			
			print.table(R,digits = 2,justify = "right",width = 1)
			
			Result = list()
			Result$method.name = Q[2:nrow(R),1]
			Result$method.log.evidence = Q[2:nrow(R),2]
			
			options(warn = 0,digits = 7)
			
			return(Result)
		
		}

#### HANDLES TO FUNCTIONS TO COMPUTE VARIOUS LOG EVIDENCES ####

evidence.obj = function(DataObs,DataCovariates,PriorMeanMu,PriorSignalMu,PriorShapePrecision,PriorRatePrecision){
	
	#args:	DataObs- observed data (dependent variable)
	#		DataCovariates- design matrix *with* a column of ones (independent variables) if constant to be included
	#		PriorMeanMu - the mean of Mu in the prior
	#		PriorSignalMu - the ratio of the precision of the prior on Mu to the error precision
	#		PriorShapePrecision - the shape of the gamma prior for the error precision
	#		PriorRatePrecision - the rate (inverse scale) of the gamma prior for the error precision
	
	#Throughout this function theta is understood to be the vector (beta',possibly.transformed(tau))'
	#		where ' denotes transpose
	
	#functions in the return:
	#	log.analytic.evidence: computes the exact analytic value of the log evidence
	#	log.laplace.evidence: computes the Laplace approximation to the evidence
	#	log.laplace.evidence.map: computes the Laplace approximation at the MAP of a Monte Carlo run
	#	log.laplace.evidence.mle: computes the Laplace approximation a the MLE
	#	log.harmonic.mean.estimator.evidence: runs a Monte Carlo algorithm and computes the corresponding
	#		harmonic mean estimator of the log evidence
	#	log.chib.evidence: computes the estimator of the marginal likelihood due to Chib (1995)
	#	log.annealed.importance.sampling.evidence: computes the estimator of the marginal likelihood using
	#		annealed importance sampling (Neal 2001)
	#	log.power.posterior.evidence: computes the estimator of the log marginal likelihood using power posteriors
	#		(Friel and Pettit, 2008)
	#	log.nested.sampling.evidence: computes the estimator of the log marginal likelihood using Nested sampling (Skilling 2008)
	#	gibbs.sampler: runs a Gibbs sampler on the linear model for a given number of iterations and burn-in
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
	X = as.matrix(DataCovariates)
	tX = t(X)
	d = ncol(X)
	mu0 = PriorMeanMu
	tau0 = PriorSignalMu
	alpha = 2*PriorShapePrecision
	delta = 2*PriorRatePrecision
	
	#convenient constants
	logPi = log(pi)
	log2Pi = log(2*pi)
	XTX = tX%*%X
	XTy = tX%*%y
	M = XTX + tau0
	cholM = chol(M)
	log.detM = 2*sum(log(diag(cholM)))
	invM = chol2inv(cholM)
	choltau0 = chol(tau0)
	invtau0 = chol2inv(choltau0)
	log.dettau0 = 2*(sum(log(diag(choltau0))))
	P = diag(1,n) - X%*%invM%*%tX
	beta0 = invM%*%(tX%*%y + tau0%*%mu0)
	yTy = t(y)%*%y
	c0 = yTy + t(mu0)%*%(tau0%*%mu0) - t(beta0)%*%M%*%beta0
	c1 = t(mu0)%*%tau0%*%mu0
	
	evidence.obj = list()
	
	#make a call to all the functions and return an annotaed matrix with the log marginal values
	
	#### ANALYTIC CALCULATION FOR THE GAUSSIAN LINEAR MODEL ####
	
	evidence.obj$log.analytic.evidence = function(){
		
		log.evidence = -0.5*n*logPi + 0.5*log.dettau0 - 0.5*log.detM + 0.5*alpha*log(delta) + lgamma((n+alpha)/2) - lgamma(alpha/2) -0.5*(n+alpha)*log(c0+delta)
											
		return(log.evidence)
		
		}
	
	#### LAPLACE APPROXIMATION USING A NEWTON ALGORITHM TO FIND MODE ####
		
	evidence.obj$log.laplace.evidence = function(){
				
		theta = evidence.obj$newton.method()	
		
		G = -evidence.obj$hessian(theta)
		
		log.detG = log(det(G))
		
		log.laplace.evidence = 0.5*(d+1)*log2Pi - 0.5*log.detG + evidence.obj$log.posterior(theta)
		
		return(log.laplace.evidence)
		
		}

	#### LAPLACE APPROXIMATION USING MAP OF GIBBS RUN AS THE MODE ####

	evidence.obj$log.laplace.evidence.map = function(Its = 5500,BurnIn = 500){
		
		#using the map of a Gibbs run as the mode for Laplace
		
		fix = list()
		fix$vars = rep(FALSE,d+1)
		fix$values = numeric(d+1)
		
		T = evidence.obj$gibbs.sampler(Its,BurnIn,fix)
		
		log.posterior = numeric(Its - BurnIn)
		
		log.posterior = apply(t(T),2,function(x) evidence.obj$log.posterior(x))
		
		j.map = which(log.posterior == max(log.posterior))
		
		theta = t(T[j.map,])

		G = -evidence.obj$hessian(theta)
		
		log.detG = log(det(G))
		
		log.laplace.evidence.map = 0.5*(d+1)*log2Pi - 0.5*log.detG + evidence.obj$log.posterior(theta)
		
		return(log.laplace.evidence.map)
		
		}

	#### LAPLACE APPROXIMATION USING THE MLE AS THE MODE ####

	evidence.obj$log.laplace.evidence.mle = function(){
		
		#using the MLE as a mode for Laplace
		
		betaHat = solve(t(X)%*%X,method = "chol")%*%(t(X)%*%y)
		
		z = y - X%*%betaHat
		
		tauHat = n/(t(z)%*%z)
		
		theta = c(betaHat,tauHat)
		
		G = -evidence.obj$hessian(theta)
		
		log.detG = log(det(G))
		
		if(is.nan(log.detG)){
			cat("\n*The Laplace approximation at the MLE was not computed since the negative Hessian at the MLE was not invertible.*\n")
			return("-")
			}
		
		log.laplace.evidence.mle = 0.5*(d+1)*log2Pi - 0.5*log.detG + evidence.obj$log.posterior(theta)
		
		return(log.laplace.evidence.mle)		
		
		}
		
	#### HARMONIC MEAN ESTIMATOR ####	
		
	evidence.obj$log.harmonic.mean.estimator.evidence = function(Its = 5500,BurnIn = 500){
		
		fix = list()
		fix$vars = rep(FALSE,d+1)
		fix$values = numeric(d+1)
		
		T = evidence.obj$gibbs.sampler(Its,BurnIn,fix)
		
		logLikelihood = apply(t(T),2,function(x) evidence.obj$log.likelihood(x))

		#try to get stability
		
		a = min(logLikelihood)
		
		harmonic.mean.estimator = log(n) + a - log(sum(exp(-logLikelihood + a)))
		
		return(harmonic.mean.estimator)
		
		}
	
	#### CHIB'S METHOD ####
		
	evidence.obj$log.chib.evidence = function(Its = 5500,BurnIn = 500){
		
		#make a matrix to store all results
		T = matrix(nrow = Its - BurnIn,ncol = d+1)
		
		modal.densities = numeric(d+1)
		
		fix = list()
		fix$vars = rep(FALSE,d+1)
		fix$values = numeric(d+1)
		
		#do full gibbs to begin
		T[,1:(d+1)] = evidence.obj$gibbs.sampler(Its,BurnIn,fix)
		
		#find the mode of the precision
		
		sh = 0.5*(n+d+alpha)
		
		dens = list()
		dens$pts = T[,d+1]
		dens$density = numeric(Its - BurnIn)
		
		rt = numeric(Its - BurnIn)
		
		for(j in 1:(Its - BurnIn)){
			
			beta = T[j,1:d]
			
			rt[j] = 0.5*( t(beta-beta0)%*%M%*%(beta-beta0) + c0 + delta )
			
			}
		
		for(It in 1:(Its - BurnIn)){
			
			tau = T[It,d+1]
			
			sum = sum(dgamma(tau,shape=sh,rate=rt))
				
			dens$density[It] = sum/(Its - BurnIn)
			
			}
				
		max.tau = dens$pts[which(dens$density == max(dens$density))]
		
		modal.densities[d+1] = max(dens$density)
		
		#fix tau and sample other variables in turn
		
		fix$vars[d+1] = TRUE
		fix$values[d+1] = max.tau
		
		for(i in 1:d){
			
			#cycle through each of the variables
			
			T = evidence.obj$gibbs.sampler(Its-BurnIn,0,fix)
			
			mu = numeric(Its - BurnIn)
			
			sdev = rep(sqrt(1/(max.tau*M[i,i])),Its - BurnIn)
			
			for(j in 1:(Its - BurnIn)){
				
				beta = T[j,1:d]
				
				mu[j] = beta0[i] - 0.5*(M[i,]%*%(beta-beta0) - M[i,i]*(beta[i]-beta0[i]))/M[i,i]
				
				}
			
			dens$pts = T[,i]
			
			for(j in 1:(Its - BurnIn)){
				
				b = T[j,i]
				
				sum = sum(dnorm(b,mean = mu,sd = sdev))
				
				dens$density[j] = sum/(Its - BurnIn)
				
				}
						
			mode = dens$pts[which(dens$density == max(dens$density))]
			
			fix$vars[i] = TRUE
			fix$values[i] = mode
			
			modal.densities[i] = max(dens$density)
			
			}
		
		log.chib.evidence = evidence.obj$log.posterior(fix$values) - sum(log(modal.densities))
		
		return(log.chib.evidence)


		}
		
	#### ANNEALED IMPORTANCE SAMPLING ####	
		
	evidence.obj$log.annealed.importance.sampling.evidence = function(sample.size = 10000,temps = seq(1,0,by=-0.01)^5){
		
		n.temps = length(temps)
		
		#here the index 1 corresponds to tempering 1.
		
		N = n.temps
		
		S = array(dim = c(n.temps,d,d))
		
		b0 = array(dim = c(d,n.temps))
		
		for(i in 1:N){
			
			S[i,,] = temps[i]*XTX + tau0 
			
			}
		
		cholS = array(dim = c(n.temps,d,d))
		
		for(i in 1:N){
			
			cholS[i,,] = chol(S[i,,])
			
			invS = chol2inv(cholS[i,,])
			
			b0[,i] = invS%*%(temps[i]*XTy + tau0%*%mu0) 			
			}
			
		const = numeric(N)
		
		for(i in 1:N){
			
			const[i] = temps[i]*yTy + t(mu0)%*%tau0%*%mu0 + delta - t(b0[,i])%*%S[i,,]%*%b0[,i] 			
			}
		
		weights = numeric(sample.size)
		
		for(r in 1:sample.size){
			
			log.weight = 0
			
			tau = rgamma(1,shape = alpha/2,rate = delta/2)
		
			beta = rnorm(d,mean = mu0,sd = sqrt(1/(tau*diag(tau0))))
			
			theta = c(beta,tau)
						
			if(N < 4) stop("Terminated Annealed Importance Sampling. Reason: Not a large enough temperature scale")
		
			for(t in (N-1):1){

			log.weight = log.weight + (temps[t]- temps[t+1])*evidence.obj$log.likelihood(theta)
				
				if(t>1){
							
				for(reps in 1:5){

				ep = backsolve(sqrt(tau)*cholS[t,,],rnorm(d),k=d)
					
				beta = b0[,t] + ep
				
				sh = 0.5*(n*temps[t] + d + alpha)
								
				rt = 0.5*(t(ep)%*%S[t,,]%*%ep + const[t])
				
				tau = rgamma(1,shape = sh,rate = rt)
				
				theta = c(beta,tau)
				
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
	
	evidence.obj$log.nested.sampling.evidence = function(N = 1000,Nstep = 500,tolerance = 0.0001){
		
		#No. its within Metropolis algorithm to get new truncated pt
		
		logTol = log(tolerance)
		
		M = 300
				
		logw = log(1-exp(-1/N))
		
		Theta = matrix(nrow = d+1,ncol = N)
		
		log.likelihoods = numeric(N)
		
		#do initial draws from the prior
		
		for(j in 1:N){
			
			tau = rgamma(1,shape = alpha/2,rate = delta/2)
		
			beta = rnorm(d,mean = mu0,sd = sqrt(1/(tau*diag(tau0))))
			
			Theta[,j] = c(beta,tau)	
			
			log.likelihoods[j] = evidence.obj$log.likelihood(Theta[,j])
			
			}
				
		logZ = -.Machine$double.xmax
		
		log.contribution = 0
		
		#print(logZ)
		
		sh = d+alpha
		
		nest = 1
		
		breakwhile = FALSE
		
		#same termination criterion as Chopin and Robert (2010)
		
		while((breakwhile == FALSE )){
						
			r = which(log.likelihoods == min(log.likelihoods))
			
			l.truncate = log.likelihoods[r]
						
			log.contribution = logw + l.truncate
			
			logZnew = evidence.obj$PLUS(logZ,log.contribution)

			if(log.contribution < -8*log(10) + logZ){
				
				breakwhile = TRUE
								
				}
			
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

		
		logw = -nest/N - log(N)
		
		for(i in 1:N){
			
			logZnew = evidence.obj$PLUS(logZ,logw+log.likelihoods[i])
			
			logZ = logZnew
			
			}
			

		return(logZ)	
		
		}
		

	evidence.obj$PLUS = function(a,b){

		if(a > b){
			
			z = a + log(1+exp(b-a))
			
			}else{
				
				z = b + log(1+exp(a-b))
				
				}
				
		return(z)
		
		}
				

	evidence.obj$evolve = function(theta,l){
		
		currentLJL = evidence.obj$log.prior(theta)
		
		acceptances = 0
		
		tau.current = theta[d+1]
		
		beta.current = theta[1:d]
		
		k = 10
		
		while(acceptances<10){
		
			tau.prop = rgamma(1,shape = k^2,rate = k^2/tau.current)
			
			beta.prop = beta.current + backsolve(sqrt(10000)*choltau0,rnorm(d))
			
			theta.prop = c(beta.prop,tau.prop)
						
			proposedLJL = evidence.obj$log.prior(theta.prop)
			
			logRatio = proposedLJL - currentLJL + dgamma(tau.current,shape = k^2,rate=k^2/tau.prop,log=TRUE) - dgamma(tau.prop,shape=k^2,rate=k^2/tau.current,log=TRUE)
			
			lnew = evidence.obj$log.likelihood(theta.prop)
			
			if(lnew > l && (logRatio > 0 || logRatio > log(runif(1)))){
				
				acceptances = acceptances + 1
				
				theta = theta.prop
				
				currentLJL = proposedLJL
				
				}
				
		}
		
		return(theta)
		
		}
		

	#### POWER POSTERIORS ####

evidence.obj$log.power.posterior.evidence = function(Its = 25000,BurnIn = 5000,temps = seq(0,1,by=0.01)^5){
		
	    N = length(temps)
		
		#cycle through temperatures estimating the expected deviance each time
		
		expected.deviance = numeric(N)
		
		S = array(dim = c(N,d,d))
		
		b0 = array(dim = c(d,N))
		
		for(i in 1:N){
			
			S[i,,] = temps[i]*XTX + tau0
			
			}
		
		cholS = array(dim = c(N,d,d))
		
		for(i in 1:N){
			
			cholS[i,,] = chol(S[i,,])
			
			invS = chol2inv(cholS[i,,])
			
			b0[,i] = invS%*%(temps[i]*XTy + tau0%*%mu0)
			
			}
			
		const = numeric(N)
		
		for(i in 1:N){
			
			const[i] = temps[i]*yTy + t(mu0)%*%tau0%*%mu0 + delta - t(b0[,i])%*%S[i,,]%*%b0[,i]
			
			}				
		
		#draw from prior to initialize first chain
		
		tau = rgamma(1,shape = alpha/2,rate = delta/2)
		
		beta = rnorm(d,mean = mu0,sd = sqrt(1/(tau*diag(tau0))))		
		
		for(t in 1:N){
			
			expected.theta = numeric(d+1)
			
			for(It in 1:Its){
					
					ep = backsolve(sqrt(tau)*cholS[t,,],rnorm(d),k=d)
					
					beta = b0[,t] + ep
									
					sh = 0.5*(n*temps[t] + d + alpha)
				
					diff = beta - b0[,t]
				
					rt = 0.5*(t(diff)%*%S[t,,]%*%diff + const[t])
					
					tau = rgamma(1,shape = sh,rate = rt)
				
					theta = c(beta,tau)
				
				
				
			if(It > BurnIn){
				
				expected.deviance[t] = expected.deviance[t] + evidence.obj$log.likelihood(theta)
				expected.theta = ((It-BurnIn-1)*expected.theta + theta)/(It-BurnIn)
				
			}
			
			}
			
			#initialize new chain with mean value of old
			theta = expected.theta
			
			}
		
		expected.deviance = expected.deviance/(Its-BurnIn)
		
		integration.weights = numeric(N-1)
	    function.vals = numeric(N-1)
		
		for(i in 1:(N-1)){
			
			integration.weights[i] = temps[i+1]-temps[i]
			
			function.vals[i] = 0.5*(expected.deviance[i+1] + expected.deviance[i])
			
			}
			
		log.evidence = t(integration.weights)%*%function.vals	
		return(log.evidence)
		
		
		}
	
	
	#### GIBBS SAMPLER FOR REGRESSION PARAMETERS AND PRECISION ####
	
	evidence.obj$gibbs.sampler = function(Its,BurnIn,fix,initial = NULL,return.log.posterior=FALSE,return.log.likelihood = FALSE){
		
		#do site by site updates for fair comparison between methods
		
		T = matrix(nrow = Its - BurnIn,ncol = d+1)
		
		#inialize from prior
		
		if(is.null(initial)){
		
		tau = rgamma(1,shape = alpha/2,rate = delta/2)
		
		beta = rnorm(d,mean = mu0,sd = sqrt(1/(tau*diag(tau0))))
		
		}else{
			
			beta = initial[1:d]
			
			tau = intial[d+1]
			
			}
		
		sh = 0.5*(n+d+alpha)
		
		sample.vars = which(fix$vars[1:d] == FALSE)
		
		fix.vars = which(fix$vars[1:d] == TRUE)
		if(length(fix.vars)>0) beta[fix.vars] = fix$values[fix.vars]
		
		if(fix$vars[d+1] == TRUE) tau = fix$values[d+1]
		
		sample.tau = !fix$vars[d+1]
		
		for(ItNum in 1:Its){
			
			#visit each parameter in turn
			
			for(j in sample.vars){
				
				w = M[j,]%*%(beta-beta0) - M[j,j]*(beta[j]-beta0[j])
				
				mu = beta0[j] - w/M[j,j]
				
				sig = sqrt(1/(tau*M[j,j]))
				
				beta[j] = rnorm(1,mean=mu,sd=sig)
				
				
				}
			
			rt = 0.5*( t(beta-beta0)%*%M%*%(beta-beta0) + c0 + delta )
			
			if(sample.tau) tau = rgamma(1,shape=sh,rate = rt)
			
			if(ItNum > BurnIn){
				
				T[ItNum-BurnIn,] = c(beta,tau)
				
				}
			
			}		
		
		return(T)
		
		}
		
	#### FUNCTION FOR EVALUATING THE LOG LIKELIHOOD ####
		
	evidence.obj$log.likelihood = function(theta){
		
		beta = theta[1:d]
		
		tau = theta[d+1]
		
		z = y - X%*%beta
		
		logLikelihood = -0.5*n*log2Pi + 0.5*n*log(tau) - 0.5*tau*t(z)%*%z
		
		return(logLikelihood)
		
		}	
	
	#### FUNCTION FOR EVALUATING THE LOG POSTERIOR ####
		
	evidence.obj$log.posterior = function(theta){
		
		dist = theta[1:d] - beta0
		
		tau = theta[d+1]
		
		logPosterior = -0.5*(n+d)*log2Pi + 0.5*(n+d)*log(tau) + 0.5*log.dettau0 -0.5*tau*(t(dist)%*%M%*%dist) - 0.5*tau*c0 +0.5*alpha*log(0.5*delta) - lgamma(0.5*alpha) + (0.5*alpha-1)*log(tau) - 0.5*delta*tau
		
		return(logPosterior)
		
		}
		
	#### FUNCTION FOR EVALUATING THE LOG PRIOR ####	
		
	evidence.obj$log.prior = function(theta){
		
		dist = theta[1:d] - beta0
		
		tau = theta[d+1]
		
		logPrior = -0.5*(d)*log2Pi + 0.5*d*log(tau)+0.5*log.dettau0 - 0.5*tau*t(dist)%*%tau0%*%dist + dgamma(tau,shape = 0.5*alpha,rate = 0.5*delta,log=TRUE)
		
		return(logPrior)
		
		}
		
	#### FUNCTION FOR EVALUATING TEMPERED LOG PRIOR POSTERIOR MIXTURE ####
	
	evidence.obj$log.ais.posterior = function(theta,t){
		
		logAISPosterior = (1-t)*evidence.obj$log.prior(theta) + t*evidence.obj$log.posterior(theta)
		
		return(logAISPosterior)
		
		}

	#### FUNCTION FOR EVALUATING THE LOG POSTERIOR GRADIENT ####
		
	evidence.obj$log.posterior.gradient = function(theta){
		
		dist = theta[1:d] - beta0
		
		tau = theta[d+1]
		
		w = -tau*(M%*%dist)
		
		z = 0.5*(n+d)/tau - 0.5*t(dist)%*%M%*%dist - 0.5*c0 + (0.5*alpha-1)/tau - 0.5*delta
		
		return(c(w,z))
		
		}	
	
	#### FUNCTION TO COMPUTE THE HESSIAN ####

	evidence.obj$hessian = function(theta){
		
		dist = theta[1:d] - beta0
		
		tau = theta[d+1]	
		
		H = matrix(nrow=d+1,ncol=d+1)
		
		H[1:d,1:d] = -tau*M
		
		z = -M%*%dist
		
		H[1:d,d+1] = z
		
		H[d+1,1:d] = t(z)
		
		H[d+1,d+1] = -0.5*(n+d)/(tau^2) - (0.5*alpha-1)/(tau^2)
		
		#print(-H)
		
		#print(solve(-H))
		
		return(H)
		
		}
	
	#### NEWTON METHOD FOR OPTIMIZING THE LOG POSTERIOR ####
	
	evidence.obj$newton.method = function(tolerance = 0.01,maxsteps = 200){
		
		#use the MLE as a starting value
		
		betaHat = solve(t(X)%*%X,method = "chol")%*%(t(X)%*%y)
		
		z = y - X%*%betaHat
		
		tauHat = n/(t(z)%*%z)
		
		theta = c(betaHat,tauHat)
				
		numsteps = 0
		
		tolcriterion = 100
		
		step.size = 1
				
		while(tolcriterion>tolerance && numsteps < maxsteps){
						
			G = evidence.obj$hessian(theta)
						
			invG = solve(G)
			
			thetaNew = theta - step.size*invG%*%evidence.obj$log.posterior.gradient(theta)
						
			if(thetaNew[d+1] < 0 || evidence.obj$log.posterior(thetaNew) < evidence.obj$log.posterior(theta)){
				cat('\nBreaking out of Newton scheme: Tolerance reached on log scale =', tolcriterion)
				return(theta)
				}
			
			tolcriterion = abs(evidence.obj$log.posterior(thetaNew)-evidence.obj$log.posterior(theta))
			
			theta = thetaNew
			
			numsteps = numsteps+1
			
			}
			
			if(numsteps == maxsteps) warning('Maximum number of steps reached in Newton method\n\tFailed to reach tolerance level\n\tThis may affect accuracy of results')

			return(theta)
		
		}	
		
		
		evidence.obj$ldensmvnorm = function(x,mu,precision,chol.precision){
			

			D = length(x)
			
			deter = 2.0*sum(log(diag(chol.precision)))
			
			mahal = t(x-mu)%*%precision%*%(x-mu)
			
			dens = -0.5*D*log2Pi + deter - 0.5*mahal
			
			return(dens)
			
			}
		
		return(evidence.obj)
		
}


			
