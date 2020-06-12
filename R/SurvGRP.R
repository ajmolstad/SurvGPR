# ------------------------------------------------
#
# MM algorithm for variance components models
#
# ------------------------------------------------

MM_Alg <- function(Z, Y, K, max.iter){
	
	# -- MM Algorithm for variance components ----
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


# --------------------------------------------
#
# Inner MCEM algorithm for kern.type = "K+I"
#
# --------------------------------------------


SurvGPR_KI = function(time, status, Z, K, tol, max.iter, 
					max.iter.MM, 
                    quiet, max.samples, initializer){

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

  
  result <- list("beta" = beta.temp, "sigma2" = alpha2.temp, "Tout" = Yup[match(1:length(time), c(train.obs, train.cen))], #"negloglik" = negloglikout, 
                 "Yimpute" = Y.tot)
  return(result)

}




# ---------------------------------------------
#
#  MCEM algorithm for kern.type = "multi-K"
#
# ---------------------------------------------

SurvGPR_MK = function(time, status, Z, K, tol, 
                      max.iter, 
                      max.iter.MM, 
                      quiet, max.samples, 
                      initializer){
  
  # ----------------------------------------
  # divide training and testing sets
  # ----------------------------------------
  train.obs <- which(status==1)
  train.cen <- which(status==0)
  train.inds <- 1:length(time)
  M <- dim(K)[3]
    
  # ------------------------------------------
  # IPW-imputation for initial values 
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
  Z0 <- Z[train.obs,,drop=FALSE]
  Z1 <- Z[train.cen,,drop=FALSE]
  Z.train <- rbind(Z0, Z1)
  Y.train <- Y.tot[c(train.obs, train.cen)]

  # -- error variance equal to unconditional variance, coefs nonzero 
  if(initializer == 0){
    
    alpha2.temp <- rep(var(log(Y.train))/M, M)
    Omega.temp <- matrix(0, nrow = dim(K)[1], ncol = dim(K)[1])
    K.chols <- array(0, dim=c(length(train.inds), length(train.inds), M))
    for(k in 1:M){
      Omega.temp <- Omega.temp + alpha2.temp[k]*K[,,k]
      K.chols[,,k] <- chol(K[c(train.obs, train.cen), c(train.obs, train.cen), k])
    }
    O.temp <- chol2inv(chol(Omega.temp[c(train.obs,train.cen), c(train.obs,train.cen)]))
    alpha2.iter <- alpha2.temp
    Omega.iter <- Omega.temp
    inner <- crossprod(Z.train, solve(Omega.iter[c(train.obs, train.cen),c(train.obs, train.cen)]))
    beta.iter <- ginv(tcrossprod(inner, t(Z.train)))%*%tcrossprod(inner, t(log(Y.train)))
    
  } 
  
  if(initializer == 1){
    
    out <- MM_Alg(Z.train, log(Y.train), K[c(train.obs, train.cen),c(train.obs, train.cen), ], max.iter = 1e3)
    Omega.temp <- matrix(0, nrow = dim(K)[1], ncol = dim(K)[1])
    beta.iter <- out$beta
    alpha2.temp <- out$sigma2
    K.chols <- array(0, dim=c(length(train.inds), length(train.inds), M))
    for(k in 1:M){
      Omega.temp <- Omega.temp + alpha2.temp[k]*K[,,k]
      K.chols[,,k] <- chol(K[c(train.obs, train.cen), c(train.obs, train.cen), k])
    }
    O.temp <- chol2inv(chol(Omega.temp[c(train.obs,train.cen), c(train.obs,train.cen)]))
    alpha2.iter <- alpha2.temp
    Omega.iter <- Omega.temp
    
  }


  loglik <- rep(0, max.iter)
  sigmas.iterates <- matrix(0, nrow=M, ncol=max.iter)
  nsamples.iterates <- rep(0, ncol=max.iter)
  Kbeta <- t(solve(Omega.iter[train.obs, train.obs], Omega.iter[train.obs, train.cen]))
  VtY <- Omega.iter[train.cen, train.cen] - t(solve(Omega.iter[train.obs, train.obs], Omega.iter[train.obs, train.cen]))%*%Omega.iter[train.obs, train.cen]
  O.iter <- solve(Omega.iter[c(train.obs,train.cen), c(train.obs,train.cen)])
  converged <- FALSE
  broken <- FALSE
  nsamples <- 500
  nsamples.inner.count <- 1

  # ---------------------------------------------------
  # iterations 
  # ---------------------------------------------------
  
  for(qq in 1:max.iter){
    
    # ---------------------------------------------------
    # E step 
    # ---------------------------------------------------
    EtY <- Z1%*%beta.iter + Kbeta%*%(Y0 - Z0%*%beta.iter)
    Yl <- rtmvnorm(nsamples, mean = c(EtY), sigma = VtY,
                  lower=Y1, upper=rep(Inf, length(Y1)), start.value = Yl[dim(Yl)[1],], 
                  algorithm = "gibbs", burn.in.samples = 1000, thinning = 10)
    barY <- matrix(colSums(Yl)/dim(Yl)[1], ncol=1)
    Yup <- matrix(c(Y0, barY), ncol=1)
    
    # ---------------------------------------------------
    # M-step
    # ---------------------------------------------------
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
      out <- matrix(0, nrow=nsamples, ncol=length(Yup))
      out[, 1:length(Y0)] <- matrix(rep(Y0,nsamples), byrow=TRUE, nrow=nsamples)
      out[,(length(Y0)+1):length(Yup)] <- Yl
      
      while(update.MM){         
        
        # -- solve beta -------
        inner <- crossprod(Z.train, O.temp)
        beta.temp <- solve(inner%*%Z.train, inner%*%Yup)
        store <- Z.train%*%beta.temp
        inner <- out - tcrossprod(rep(1, nsamples), store)
        
        # --- update sigmas -----
        alpha2.up <- rep(0, M)
        inner.temp <- tcrossprod(inner, O.temp)
        Omega.temp <- matrix(0, nrow=dim(Omega.iter)[1], ncol=dim(Omega.iter)[1])
        for(k in 1:M){
          if(alpha2.temp[k] > 1e-8){
            l.mat <- tcrossprod(inner.temp, K.chols[,,k])
            alpha2.up[k] <- alpha2.temp[k]*sqrt((nsamples^(-1)*sum(l.mat^2))/sum(O.temp*K[c(train.obs, train.cen),c(train.obs, train.cen),k]))
            Omega.temp <- Omega.temp + alpha2.up[k]*K[,,k]
          }
        }

        alpha2.temp <- alpha2.up
        O.temp <- chol2inv(chol(Omega.temp[c(train.obs,train.cen), c(train.obs,train.cen)]))
        Omega.det <- determinant(Omega.temp[c(train.obs,train.cen), c(train.obs,train.cen)], logarithm=TRUE)$modulus[1]
        out.temp  <- rowSums(tcrossprod(inner, chol(O.temp))^2)
        loglik.temp <- - sum(out.temp) - nsamples*Omega.det

        #if(M > 2){
         # --- extrapolation attempt ---
          if(iter.MM > 1){

            alpha2.extrapolation <- alpha2.temp + (1/sqrt(iter.MM-1))*(alpha2.temp - alpha2.temp.old)
            Omega.extrapolation <- matrix(0, nrow=dim(Omega.iter)[1], ncol=dim(Omega.iter)[1])
            for(k in 1:M){
              Omega.extrapolation <- Omega.extrapolation + alpha2.extrapolation[k]*K[,,k]
            }
            O.extrapolation <- chol2inv(chol(Omega.extrapolation[c(train.obs,train.cen), c(train.obs,train.cen)]))
            Omega.det.extrapolation <- determinant(Omega.extrapolation[c(train.obs,train.cen), c(train.obs,train.cen)], logarithm=TRUE)$modulus[1]
            out.extrapolation  <- rowSums(tcrossprod(inner, chol(O.extrapolation))^2)
            loglik.extrapolation <- - sum(out.extrapolation) - nsamples*Omega.det.extrapolation
            
            # ----- cat("Extrapolation step: extrap loglik = ",loglik.extrapolation, "prev loglik = ", loglik.temp, "\n")
            if(loglik.extrapolation > loglik.temp){
              
              alpha2.temp <- alpha2.extrapolation
              Omega.temp <- Omega.extrapolation
              O.temp <- O.extrapolation
              Omega.det <- Omega.det.extrapolation
              loglik.temp <- loglik.extrapolation

            }

          }
         #  }
        # cat(loglik.temp, "\n")

        if(iter.MM == 1){
          loglik.orig <- loglik.temp
        }

        if(iter.MM > 1){
          if(abs(loglik.temp - loglik.oldtemp)/abs(loglik.orig) < tol){
            update.MM <- FALSE
            loglik.new <- loglik.temp
          } 
        }
        
        loglik.oldtemp <- loglik.temp
        iter.MM <- iter.MM + 1
        alpha2.temp.old <- alpha2.temp
        
        if(iter.MM > max.iter.MM){
          update.MM <- FALSE
        }
      }


      loglik.new <- loglik.temp
      loglik[qq] <- loglik.new
      z.train.iter <- Z.train%*%beta.iter
      inner <- out - tcrossprod(rep(1, nsamples), z.train.iter)
      out.old <- rowSums(tcrossprod(inner, chol(O.iter))^2)
      old.det <- determinant(Omega.iter[c(train.obs,train.cen), c(train.obs,train.cen)], logarithm=TRUE)$modulus[1]
      loglik.old <- - sum(out.old) - nsamples*old.det
      sigmas.iterates[,qq] <-alpha2.iter
      nsamples.iterates[qq] <- nsamples

      ASE <- mcse(- .5*out.temp - .5*Omega.det + .5*out.old + .5*old.det, method="tukey")$se
      if(!quiet){
        cat(qq,": ASE =", round(ASE,3), "; sigma2 =", round(alpha2.temp[1], 3), "; resid =", round(.5*loglik.new/nsamples - .5*loglik.old/nsamples - 2.58*ASE, 5),"; sk =", nsamples, "\n")
      }        
      
      if(.5*loglik.new/nsamples - .5*loglik.old/nsamples - 2.58*ASE > 1e-4 & nsamples.inner.count <= 10){
        
        update <- FALSE
        beta.iter.old <- beta.iter
        alpha2.iter.old <- alpha2.iter
        beta.iter <- beta.temp
        alpha2.iter <- alpha2.temp
        Omega.iter <- Omega.temp
        O.iter <- O.temp
        Kbeta <- t(solve(Omega.iter[train.obs, train.obs], Omega.iter[train.obs, train.cen]))
        VtY <- Omega.iter[train.cen, train.cen] - t(solve(Omega.iter[train.obs, train.obs], Omega.iter[train.obs, train.cen]))%*%Omega.iter[train.obs, train.cen]
        nsamples.inner.count <- nsamples.inner.count + 1
        inner.count <-  inner.count + 1


      } else {
        
        Yadd = rtmvnorm(nsamples, mean = c(EtY), sigma = VtY, start.value = Yl[dim(Yl)[1],],
                        lower=Y1, upper=rep(Inf, length(Y1)), algorithm = "gibbs", burn.in.samples = 1000, thinning = 10)
        nsamples <- 2*nsamples
        Yl = rbind(Yl, Yadd)
        barY = matrix(colSums(Yl)/dim(Yl)[1], ncol=1)      
        Yup = matrix(c(Y0, barY), ncol=1)
        cat("Adding samples for ascent:" , nsamples, " samples", "\n")
        nsamples.inner.count <- 1
        
        
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
  
  Omega.new <- Omega.temp
  beta.new <- beta.temp
  
  result <- list("beta" = beta.temp, "sigma2" = alpha2.temp, "Tout" = Yup[match(1:length(time), c(train.obs, train.cen))], #"negloglik" = negloglikout, 
                 "Yimpute" = Y.tot)
  return(result)
  
}




# -------------------------------------------------
#
#   Prediction function for fitted model 
#
# --------------------------------------------------

SurvGPR_Predict <- function(results, barT, Z.full, K.full, train.inds, test.inds, times = NULL){
    
    # ------------------------------------------
    # Create covariance
    # ------------------------------------------
    K.temp <- array(0, dim=c(dim(K.full)[1], dim(K.full)[2], dim(K.full)[3]+1))
    for(j in 1:dim(K.full)[3]){
      K.temp[,,j] <- K.full[,,j]
    }
    K.temp[,,dim(K.full)[3]+1] <- diag(1, dim(K.full)[1])
    K.full <- K.temp
  
    # --------------------------------------------
    # Check dimensionality
    # --------------------------------------------
    if(dim(K.full)[3] != length(results$sigma2)){
      stop("Length of results$sigma2 != dim(K.full)[3]")
    }
  
    K.mat <- matrix(0, nrow=dim(K.full)[1], ncol=dim(K.full)[2])
    for(j in 1:length(results$sigma2)){
      K.mat <- K.mat + results$sigma2[j]*K.full[,,j]
    }
    
    # --------------------------------------------
    # Generate submatrices for prediction 
    # --------------------------------------------
    K11 <- K.mat[train.inds,train.inds, drop=FALSE]
    K12 <- K.mat[train.inds,test.inds, drop=FALSE]
    K22 <- K.mat[test.inds,test.inds, drop=FALSE]
    Kbeta = t(solve(K11, K12))
    
    # ---------------------------------------------
    # Compute conditional mean E(test|train)
    # ---------------------------------------------
    EtY = Z.full[test.inds,]%*%results$beta + Kbeta%*%(barT - Z.full[train.inds, ]%*%results$beta)
    
    # -------------------------------------------------------
    # If time!=NULL, compute individual survival curves
    # -------------------------------------------------------
    if(!is.null(times)){
   
      # -------------------------------------------------
      # get conditional variance
      # --------------------------------------------------
      VtY <- K22 - t(K12)%*%solve(K11, K12)
      
      out <- matrix(0, nrow=length(test.inds), ncol=length(times))
      for(j in 1:length(test.inds)){
        for(k in 1:length(times)){
          out[j,k] <- 1 - pnorm(log(times)[k], mean = EtY[j], sd = sqrt(VtY[j,j]))
        }
      }
    
    } else {
      out <- NULL
    }
    
    return(list("test.inds" = test.inds, "log.pred" = EtY, "survFunc" = out))
}




# -----------------------------------------------
# 
#   Complete function wrapper 
# 
# ------------------------------------------------

SurvGPR = function(time, status, Z, K, tol = 1e-7, 
                    max.iter = 100, max.iter.MM = 100, 
                    quiet = FALSE, max.samples = 1e5,
                    kern.type = c("multi-K"), initializer = 0){

  if(kern.type == "multi-K"){
      # - store
      K.input <- array(0, dim=c(dim(K)[1], dim(K)[2], dim(K)[3]+1))
      for(j in 1:dim(K)[3]){
        K.input[,,j] <- K[,,j]
      }
      K.input[,,dim(K)[3]+1] <- diag(1, dim(K[,,1])[1])
      K <- K.input
      
      if(!(initializer %in% 0:3)){
        stop("Initializer needs to one of {0,1,2,3}")
      }
      
      results <- SurvGPR_MK(time = time, status = status, Z = Z, K = K.input, 
                            tol = tol, max.iter = max.iter, max.iter.MM = max.iter.MM, 
                            quiet=quiet, max.samples = max.samples, initializer = initializer)
  }

  if(kern.type=="K+I"){
      K.input <- array(0, dim=c(dim(K)[1], dim(K)[2], 2))
      K.input[,,1] <- K[,,1]
      K.input[,,2] <- diag(1, dim(K[,,1])[2])
      results <- SurvGPR_KI(time = time, status = status, Z = Z, K = K.input, tol = tol, 
                            max.iter = max.iter, max.iter.MM = max.iter.MM, 
                            quiet=quiet, max.samples = max.samples, initializer = initializer)
  }
  
  return(results)
}
