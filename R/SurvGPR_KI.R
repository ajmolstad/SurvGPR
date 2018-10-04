SurvGPR_KI = function(time, status, Z, K, tol, max.iter, 
					max.iter.MM, 
                    quiet, max.samples){

	# ----------------------------------------
	# censored and uncensored observations 
	# ----------------------------------------
	train.obs <- which(status==1)
	train.cen <- which(status==0)
	M <- 2

	# ------------------------------------------
	# IPW mean-imputation for initial values  (Datta, 2005)
	# ------------------------------------------
	Y.train <- time
	time.train <- time
	status.train <- status
	time.sorted <- time.train[sort(time.train, index.return = TRUE)$ix]
	status.sorted <- status[sort(time.train, index.return = TRUE)$ix]
	status.sorted[length(status.sorted)] = 1

	failure.times <- unique(time.sorted[status.sorted==1])
	temp <- rep(0, length(failure.times))
	KM <- rep(0, length(failure.times))
	for(j in 1:length(failure.times)){
		temp[j] <- 1 - sum(time.sorted[status.sorted == 1] == failure.times[j])/sum(time.sorted >= failure.times[j])
		KM[j] <- prod(temp[1:j])
	}

	time.observed <- failure.times
	time.imputed <- time.sorted[status.sorted==0]
	time.updated <- time.imputed

	for(j in 1:length(time.imputed)){
		inds <- which(time.observed > time.imputed[j])
		intermed <- rep(0, length(inds))
		for(k in 1:length(inds)){
		  if(inds[k] == 1){
		    intermed[k] <- time.observed[inds[k]]*(1 - KM[inds[k]])
		  } else {
		    intermed[k] <- time.observed[inds[k]]*(KM[inds[k]-1] - KM[inds[k]])
		  }
		}

		if(any(time.observed == time.imputed[j])){
		  time.updated[j] <- sum(intermed)/KM[which(time.observed == time.imputed[j])]
		} else {
		  if(length(which(time.observed < time.imputed[j]))==0){
		    time.updated[j] <- sum(intermed)
		  } else {
		    time.updated[j] <- sum(intermed)/KM[max(which(time.observed < time.imputed[j]))]
		  } 
		}
	}

	time.sorted[status.sorted==0] <- time.updated
	Y.train[sort(time.train, index.return = TRUE)$ix] <- time.sorted
	Y.train <- pmax(Y.train, 1)
	Y.tot <- rep(0, length(time))
	Y.tot <- Y.train
	Yl <- matrix(log(Y.tot[train.cen]), nrow=1)
	# ----------------------------------------------
	# algorithm initialization 
	# ----------------------------------------------
	Y0 <- log(time[train.obs])
	Y1 <- log(pmax(time[train.cen], 1))
	Z0 <- Z[train.obs,, drop=FALSE]
	Z1 <- Z[train.cen,, drop=FALSE]
	Z.train <- rbind(Z0, Z1)
	rm.colind <- NULL
	Y.train <- Y.tot[c(train.obs, train.cen)]
	eo <- eigen(K[c(train.obs, train.cen),c(train.obs, train.cen),1])
	U <- eo$vec
	lambda <- eo$val
	Omega.temp <- matrix(0, nrow=dim(K)[1], ncol=dim(K)[1])
	alpha2.temp <- rep(1, M)
	Omega.iter <- alpha2.temp[1]*K[,,1] + alpha2.temp[2]*K[,,2]
	O.iter = solve(Omega.iter[c(train.obs,train.cen), c(train.obs,train.cen)])
	beta.iter <- solve(t(Z.train)%*%O.iter%*%Z.train, t(Z.train)%*%O.iter%*%log(Y.train))
	store <- Z.train%*%beta.iter
	alpha2.iter <- alpha2.temp
	Kbeta <- t(solve(Omega.iter[train.obs, train.obs], Omega.iter[train.obs, train.cen]))
	VtY <- Omega.iter[train.cen, train.cen] - t(solve(Omega.iter[train.obs, train.obs], Omega.iter[train.obs, train.cen]))%*%Omega.iter[train.obs, train.cen]
	converged <- FALSE
	broken <- FALSE
	nsamples <- 500
	loglik <- rep(0, max.iter)

  # -------------------------------------------------------
  # iterations 
  # -------------------------------------------------------
	for(qq in 1:max.iter){
		
		# --- E step --------
		EtY = Z1%*%beta.iter + Kbeta%*%(Y0 - Z0%*%beta.iter)
		Yl = rtmvnorm(nsamples, mean = c(EtY), sigma = VtY,  
		    lower=Y1, upper=rep(Inf, length(Y1)), start.value = Yl[dim(Yl)[1],], burn.in.samples = 1000, 
		    algorithm = "gibbs", thinning=10)
		barY = matrix(colSums(Yl)/dim(Yl)[1], ncol=1)
		Yup = matrix(c(Y0, barY), ncol=1)

		# --- M-step --------
		update <- TRUE
		beta.temp <- beta.iter
		alpha2.temp <- alpha2.iter
		O.temp <- O.iter
		inner.count <- 1

		while(update){

			update.MM <- TRUE
			iter.MM <- 1
			if(qq > 1){
				loglik.oldtemp <- loglik[qq-1]
			} else {
				loglik.oldtemp <- Inf
			}

			# --------------------------------------
			# Run blockwise coordinate descent
			# --------------------------------------
			while(update.MM){         

				# ----------------------------
				# update beta 
				# ----------------------------
				inner <- crossprod(Z.train, O.temp)
				beta.temp <- solve(inner%*%Z.train, inner%*%Yup)
				store <- Z.train%*%beta.temp

				# ----------------------------
				# update alpha 
				# ----------------------------
				out <- matrix(0, nrow=nsamples, ncol=length(Yup))
				out[, 1:length(Y0)] <- matrix(rep(Y0,nsamples), byrow=TRUE, nrow=nsamples)
				out[,(length(Y0)+1):length(Yup)] <- Yl
				inner <- out - tcrossprod(rep(1, nsamples), store)
				R <- rowSums(tcrossprod(t(U), inner)^2)/nsamples
				fn <- function(x){sum(R/(lambda*x[1] + x[2]) - log(1/(lambda*x[1] + x[2])))}
				out <- optimx(par=alpha2.temp, fn=fn, lower = c(0, 1e-6), upper=c(Inf, Inf), method="L-BFGS-B")
				alpha2.temp[1] <- out$p1
				alpha2.temp[2] <- out$p2 

				# ----------------------------
				# compute loglik 
				# ----------------------------
				Omega.temp <- alpha2.temp[1]*K[,,1] + alpha2.temp[2]*K[,,2]
				O.temp <- chol2inv(chol(Omega.temp[c(train.obs,train.cen), c(train.obs,train.cen)]))
				Omega.det <- determinant(Omega.temp[c(train.obs,train.cen), c(train.obs,train.cen)], logarithm=TRUE)$modulus[1]
				out.temp  <- rowSums(tcrossprod(inner, chol(O.temp))^2)
				loglik.temp <- - sum(out.temp) - nsamples*Omega.det
				cat(loglik.temp, "\n") # for development
				if(iter.MM == 1){
					loglik.orig <- loglik.temp
				}
				if(iter.MM > 1){
					if(abs(loglik.temp - loglik.oldtemp)/abs(loglik.orig) < tol){
						update.MM <- FALSE
						loglik.new <- loglik.temp
					} else {
						loglik.oldtemp <- loglik.temp
					}
				}
				
				loglik.oldtemp <- loglik.temp
				iter.MM <- iter.MM + 1

				if(iter.MM > max.iter.MM){
					update.MM <- FALSE
				}
			}

			# --- cat("Converged after", iter.MM, "iterations",  "\n")

			# ------------------------------
			out <- matrix(0, nrow=nsamples, ncol=length(Yup))
      		out[, 1:length(Y0)] <- matrix(rep(Y0,nsamples), byrow=TRUE, nrow=nsamples)
     		out[,(length(Y0)+1):length(Yup)] <- Yl
			loglik.new <- loglik.temp
     		loglik[qq] <- loglik.new
      		z.train.iter <- Z.train%*%beta.iter
      		inner <- out - tcrossprod(rep(1, nsamples), z.train.iter)
      		out.old <- rowSums(tcrossprod(inner, chol(O.iter))^2)
      		old.det <- determinant(Omega.iter[c(train.obs,train.cen), c(train.obs,train.cen)], logarithm=TRUE)$modulus[1]
      		loglik.old <- - sum(out.old) - nsamples*old.det
			ASE <- mcse(- .5*out.temp - .5*determinant(Omega.temp[c(train.obs,train.cen), c(train.obs,train.cen)], logarithm=TRUE)$modulus[1] + .5*out.old + .5*old.det, method="tukey")$se
			
			if(!quiet){
				cat(qq, ": ASE =", round(ASE,3), "; sigma2 =", round(alpha2.temp, 3), "; resid =", round(.5*loglik.new/nsamples - .5*loglik.old/nsamples - 1.96*ASE, 5),"; sk =", nsamples, "\n")
			}        
			
			if(c(.5*loglik.new/nsamples - .5*loglik.old/nsamples - 1.96*ASE) > 0 ){

				update <- FALSE
				beta.iter <- beta.temp
				alpha2.iter <- alpha2.temp
				Omega.iter <- Omega.temp
				O.iter <- O.temp
				Kbeta <- t(solve(Omega.iter[train.obs, train.obs], Omega.iter[train.obs, train.cen]))
				VtY <- Omega.iter[train.cen, train.cen] - t(solve(Omega.iter[train.obs, train.obs], Omega.iter[train.obs, train.cen]))%*%Omega.iter[train.obs, train.cen]

			} else {

				Yadd = rtmvnorm(nsamples, mean = c(EtY), sigma = VtY, 
				          lower=Y1, upper=rep(Inf, length(Y1)),  start.value = Yl[dim(Yl)[1],], algorithm = "gibbs", burn.in.samples = 1000, thinning=10)
				nsamples <- 2*nsamples
				Yl = rbind(Yl, Yadd)
				barY = matrix(colSums(Yl)/dim(Yl)[1], ncol=1)      
				Yup = matrix(c(Y0, barY), ncol=1)
				cat("Adding samples for ascent:" , nsamples, " samples", "\n")
				inner.count <-  inner.count + 1

				if(inner.count > max.iter){
					update <- FALSE
					broken <- TRUE
				} 
				if(nsamples > max.samples){
					update <- FALSE
					broken <- TRUE
				} 

			}
		}

		nsamples <- nsamples 
		if(nsamples > max.samples){
			broken <- TRUE
		}
		
		if(broken){
			break
		}

	}

	beta.new <- beta.temp



	# ----------------------------------------------
	# Compute AIC 
	# ----------------------------------------------
	Kbeta <- t(solve(Omega.temp[train.obs, train.obs], Omega.temp[train.obs, train.cen]))
	# --- get log-density for censored outcomes
	cond.mean <-  Z1%*%beta.temp + Kbeta%*%(Y0 - Z0%*%beta.temp)
	cond.Var <- Omega.temp[train.cen, train.cen] - t(solve(Omega.temp[train.obs, train.obs], Omega.temp[train.obs, train.cen]))%*%Omega.temp[train.obs, train.cen]
	d <- pmvnorm(lower=Y1, upper=rep(Inf, length(Y1)), mean=c(cond.mean),
	    sigma=cond.Var)$error
	cat(d, "\n")
	loglik <- log(d) -.5*determinant(Omega.temp[train.obs, train.obs], logarithm=TRUE)$modulus[1] - .5*t(Y0 - Z0%*%beta.temp)%*%solve(Omega.temp[train.obs,train.obs], Y0 - Z0%*%beta.temp) - .5*length(Y0)*log(2*pi)
	cat(loglik, "\n")
	AIC <- 2*(dim(Z)[2] + dim(K)[3]) - 2*loglik

	result <- list("beta" = beta.temp, "sigma2" = alpha2.temp, "Tout"  = Yup[match(1:length(time), c(train.obs, train.cen))], "AIC" = AIC)
	return(result)
	
}
