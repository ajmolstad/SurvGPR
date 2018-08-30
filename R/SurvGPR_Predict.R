SurvGPR_Predict <- function(results, barT, Z.full, K.full, train.inds, test.inds, times = NULL, kern.type=c("K+I", "multi-K")){
    
    # ------------------------------------------
    #
    # ------------------------------------------
  
    if(kern.type=="K+I"){
      K.temp <- array(0, dim=c(dim(K.full)[1], dim(K.full)[2], 2))
      K.temp[,,1] <- K.full[,,1]
      K.temp[,,2] <- diag(1, dim(K.full)[1])
      K.full <- K.temp
    }
    if(kern.type=="multi-K"){
      K.temp <- array(0, dim=c(dim(K.full)[1], dim(K.full)[2], dim(K.full)[3]+1))
      for(j in 1:dim(K.full)[3]){
        K.temp[,,j] <- K.full[,,j]
      }
      K.temp[,,dim(K.full)[3]+1] <- diag(1, dim(K.full)[1])
      K.full <- K.temp
    }
  
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