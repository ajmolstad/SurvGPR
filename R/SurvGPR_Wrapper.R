# ---------------------------------------------------------
# Survival GPR function 
# ---------------------------------------------------------
# Inputs: 
#     time; n-dimensional failure time/censoring time
#		  status; n-dimensional indicator of censoring 
# 		Z: n by q matrix of covariates for mean function
# 		K: M-dimensional array of n by n similarity-kernels
#		  train.inds/test.inds; indicies for train/test
# 		tol: tolerance for change in objective function value
# 		max.iter: maximum number of outer iterations
#		  kern.type: type of kernel setup
# 		quiet: Print iteration progress?
# 		upper.bound: if/how to upper bound survival time
#		  upper.bound.vec: upper bound == "age", n-dimensional vec
#		  iter.MM: the number of MM iterations per update
# ----------------------------------------------------------
# Outputs:
#     EM.time: computing time of the EM algorithm in seconds
#     EM.Cind: the concordance measure from coxph
#     EM.Cuno: the Uno AUC measure with probability of censoring 
#              computed on training data
#     EM.Cuno_val: the Uno AUC measure with probability of censoring 
#                  computed on test data
#     Ypred: log-scale predicted survival times for the test set
#     Yimp: Imputed survival times for training data
# ----------------------------------------------------------------


SurvGPR = function(time, status, Z, K, tol = 1e-7, 
                    max.iter = 100, max.iter.MM = 100, 
                    quiet = FALSE, max.samples = 1e5,
                    kern.type = c("K+I", "multi-K"), initializer = 0){

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
