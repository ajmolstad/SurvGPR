MM_Alg <- function(Z, Y, K, max.iter){
	
	# -- Alg 1 from Hua Zhou ----
	sigma2 <- rep(var(Y)/dim(K)[3], dim(K)[3])
	Omega <- matrix(0, nrow=dim(K)[1], ncol=dim(K)[2])
	for(j in 1:dim(K)[3]){
		Omega <- Omega + sigma2[j]*K[,,j]
	}
	OmegaInvZ <- solve(Omega, Z)
	beta <- solve(crossprod(OmegaInvZ, Z),crossprod(OmegaInvZ,Y))

	log.lik.orig <- crossprod(Y - Z%*%beta, solve(Omega, Y - Z%*%beta)) + determinant(Omega, logarithm=TRUE)$modulus[1]

	iterating <- TRUE
	k.iter <- 1
	log.lik.old <- log.lik.orig
	Omegat <- Omega

	while(iterating){

		OmegaInvZ <- solve(Omegat, Z)
		betat <- solve(crossprod(OmegaInvZ, Z),crossprod(OmegaInvZ,Y))
		temp <- Z%*%betat
		sigma2t <- rep(0, dim(K)[3])
		L <- solve(Omegat, Y - temp)
		for(j in 1:dim(K)[3]){
			sigma2t[j] <- sigma2[j]*sqrt(tcrossprod(crossprod(L, K[,,j]), t(L))/sum(diag(solve(Omegat, K[,,j]))))
		}
		Omegat <- matrix(0, nrow=dim(K)[1], ncol=dim(K)[2])
		for(j in 1:dim(K)[3]){
			Omegat <- Omegat + sigma2t[j]*K[,,j]
		}

		log.lik.new <- crossprod(Y - temp, solve(Omegat, Y - temp)) + determinant(Omegat, logarithm=TRUE)$modulus[1]

		# cat("k = ", k.iter, "; log-lik", log.lik.new, "\n")
		k.iter <- k.iter + 1
		if(k.iter > 2){
			if(log.lik.old - log.lik.new < 1e-8*abs(log.lik.orig)){
				iterating <- FALSE
			}
		}

		if(k.iter > max.iter){
			iterating <- FALSE
		}
		sigma2 <- sigma2t
		log.lik.old <- log.lik.new

	}

	return(list("Y" = Y, "sigma2" = sigma2t, "beta" = beta))


}



# keep <- MM_Alg(Z=Z, Y=Y, K=K, max.iter =200)


