library(Rfast)
library(foreach)
library(doParallel)
library(doSNOW)


# Algorithms for Parameter Estimation in Normal Distibutions --------------

piEst = function(E){
  require(Rfast)
  # function to calculate pi values in EM model
  # E is n by k matrix (n: number of data points, k: number of cluster groups)
  colsums(E)/dim(E)[1]
}

NormParamMu = function(X, E){
  require(Rfast)
  # function to calculate mu in a normal model
  # X is matrix of n by m (n: number of data points, m: number of dimentions/replicates)
  # E is n by k matrix (n: number of data points, k: number of cluster groups)
  
  # this first if should never realy be true, this would mean that we fit only a single model, but for colsums compatibility...
  if(is.vector(E)){
    t(E) %*% X/sum(E)
  }else{
    # calculating mu as k by m matrix (k: number of cluster groups, m: number of dimentions/replicates)
    t(E) %*% X/colsums(E)
  }
}

NormParamSigSq1dim = function(X, E, mu){
  require(Rfast)
  # function to calculate mu in a normal model for 1 dimentional cases
  # X is vector of length n (n: number of data points)
  # E is n by k matrix (n: number of data points, k: number of cluster groups)
  # mu is k by 1 matrix as calculated by NormParamMu  (k: number of cluster groups)
  nModels = dim(mu)[1]  
  nSamples = length(X)
  # Here we return a vector of length k where each element is the sigma^2 for that model
  colsums(E*((matrix(X, nrow = nSamples, ncol = nModels) - matrix(mu, nrow = nSamples, ncol = nModels, byrow = T))^2))/colsums(E)
}

NormParamSigSqMultiDim = function(X, E, mu){
  require(Rfast)
  require(stats)
  # function to calculate mu in a normal model
  # X is matrix of n by m (n: number of data points, m: number of dimentions/replicates)
  # E is n by k matrix (n: number of data points, k: number of cluster groups)
  # mu is k by m matrix as calculated by NormParamMu  (k: number of cluster groups, m: number of dimentions/replicates)
  # declaring parameters before iteration => speed
  nModels = dim(mu)[1]  
  nDim = dim(X)[2]
  sigs = array(dim = c(nDim, nDim, nModels))
  # because of high dimentions, either a list or a tensor would be required to store the Covariance matrix, I opted for a list
  # looping through the k cluster models
  for(i in 1:nModels){
    # using the stats function cov.wt to calculate covariance matrices (was faster than self-build function with high numbers)
    s = tryCatch(cov.wt(X, wt = E[,i], center = mu[i,], method = "ML")$cov, error = function(w){F})
    if(!is.matrix(s)){
      s
    }else{
      sigs[,,i] = s
    }
  }
  sigs
}




# Algorithms for Parameter Estimation and Expectation in Binomial Distribution ------------

BinomLnTheta = function(N, S, E){log(crossprod(S,E)) - log(crossprod(N,E))}

BinomLnA = function(N, S, lnT, empi){
  if(is.vector(N)){
    nSamples = length(N)
  }else{
    nSamples = dim(N)[1]
  }
  S %*% lnT + (N-S) %*% log(1-exp(lnT)) + matrix(rep(log(empi), nSamples), nrow = nSamples, byrow = T)
}


# Algorithms to Evaluate Expectation in Normal distibution ----------------

NormLnAMatrix1dim = function(X, mu, sigsq, empi){
  # X is vector of length n (n: number of data points)
  # mu is k by 1 matrix as calculated by NormParamMu  (k: number of cluster groups)
  # sigsq is a vector of length k as calculated by NormParamSigSq1dim (k: number of cluster groups)
  # empi is a vector of length k as calculated by piEst
  nModels = dim(mu)[1]  
  nSamples = length(X)
  log(matrix(empi/sqrt(sigsq), nrow = nSamples, ncol = nModels, byrow = T)) - 
    (0.5*(((matrix(X, nrow = nSamples, ncol = nModels) - 
            matrix(mu, nrow = nSamples, ncol = nModels, byrow = T))^2)/matrix(sigsq, nrow = nSamples, ncol = nModels, byrow = T)))
}

NormLnAMatrixMultidim = function(X, mu, sigsq, empi){
  # X is matrix of n by m (n: number of data points, m: number of dimentions/replicates)
  # mu is k by 1 matrix as calculated by NormParamMu  (k: number of cluster groups)
  # sigsq is a array of m by m by k matrices of length k as calculated by NormParamSigSqMultiDim (k: number of cluster groups)
  # empi is a vector of length k as calculated by piEst
  nModels = dim(mu)[1]  
  nSamples = dim(X)[1]
  nDim = dim(X)[2]
  LnA = matrix(nrow = nSamples, ncol = nModels)
  for(i in 1:nModels){
    Xdash = X - matrix(mu[i,], nrow = nSamples, ncol = nDim, byrow = T)
    LnA[,i] = rep(log(empi[i]) - 0.5*(base::determinant(as.matrix(sigsq[,,i]))$modulus), nSamples) - 0.5*rowsums((Xdash %*% solve(sigsq[,,i])) * Xdash)
  }
  LnA
}

LnARowSums = function(lnA){
  require(Rfast)
  # lnA is n by k matrix as calculated by the lnA type functions (n: number of data points, k: number of cluster groups)
  lnARmax <- rowMaxs(lnA, value = T)
  lnARSum <- lnARmax + log(rowsums(exp(lnA - lnARmax)))
}


# Wrappers for one dimensional normal distribution ------------------------



Norm1dimEMcycle = function(X, k, maxiter = 1000, maxtol = 1E-20){
  # X is vector of length n (n: number of data points)
  # k is number of components/models to by fitted
  # maxiter and maxtol are thresholds for exiting the optimization
  
  nDim = 1
  nSamples = length(X)
  
  # start bby initilazing a potential E
  mu = matrix(runif(k*nDim), nrow = k)
  sigsq = runif(k, min = 0, max = 1)
  empi = runif(k)
  empi = empi/sum(empi)
  
  # setting likelyhood of initial (random) model
  postlnL = NA 
  # starting loop to find optimum
  niter = 0
  con = F
  
  while(!con){
    
    # writing over some variables for comparison
    niter = niter + 1
    priorlnL = postlnL
    
    # Ex step
    lnA = NormLnAMatrix1dim(X, mu, sigsq, empi)
    lnArs = LnARowSums(lnA) 
    E = exp(lnA - lnArs)
    
    # Max step
    mu = NormParamMu(X, E)
    sigsq = NormParamSigSq1dim(X, E, mu)
    empi = piEst(E)
    
    
    # checking for likelyhood convergence
    postlnL = sum(lnArs)
    if(niter > 5){
      if((postlnL - priorlnL)^2 < maxtol){
        converged = T
        
        con = T
      }
    }
    if(niter > maxiter){
      converged = f
      print("EM did not converge after maximum number of iterations")
      con = T
    }
  }
  list(Expectaion = E, GroupMeans = mu, GroupVariance = sigsq, GroupSizePi = empi, ModelLnLikelyhood = postlnL, converged = converged)
}

# Wrappers for multivariate normal distributions --------------------------

NormMultiDimEMcycle = function(X, k, maxiter = 1000, maxtol = 1E-20){
  require(Rfast)
  # X is vector of length n (n: number of data points)
  # k is number of components/models to by fitted
  # maxiter and maxtol are thresholds for exiting the optimization
  
  nDim = dim(X)[2]
  nSamples = dim(X)[1]
  
  # start bby initilazing a potential E
  Xcolmin = colMins(X, value = T)
  Xcolmax = colMaxs(X, value = T)
  Xcolvarmax = max((Xcolmax - Xcolmin)^2)
  
  mu = matrix(runif(k*nDim, min = Xcolmin, max = Xcolmax), nrow = k, byrow = T)
  sigsq = array(runif(k*nDim^2, min = 0, max = Xcolvarmax), dim = c(nDim, nDim, k))
  empi = runif(k)
  empi = empi/sum(empi)
  
  # setting likelyhood of initial (random) model
  postlnL = NA 
  # starting loop to find optimum
  niter = 0
  con = F
  
  while(!con){
    
    # writing over some variables for comparison
    niter = niter + 1
    priorlnL = postlnL
    
    # Ex step
    lnA = NormLnAMatrixMultidim(X, mu, sigsq, empi)
    lnArs = LnARowSums(lnA) 
    E = exp(lnA - lnArs)
    
    # Max step
    mu = NormParamMu(X, E)
    sigsq = NormParamSigSqMultiDim(X, E, mu)
    if(!is.array(sigsq)){
      simpleError("Covariance matrix was singular")
    }
    empi = piEst(E)
    
    
    # checking for likelyhood convergence
    postlnL = sum(lnArs)
    if(niter > 5){
      if((postlnL - priorlnL)^2 < maxtol){
        converged = T
        
        con = T
      }
    }
    if(niter > maxiter){
      converged = f
      print("EM did not converge after maximum number of iterations")
      con = T
    }
  }
  list(Expectaion = E, GroupMeans = mu, GroupVariance = sigsq, GroupSizePi = empi, ModelLnLikelyhood = postlnL, converged = converged)
}

Norm1dimEMwrapper = function(X, k, maxiter = 1000, maxtol = 0.0001, ntrys = 10){
  # X is vector of length n (n: number of data points)
  # k is number of components/models to by fitted
  # maxiter and maxtol are thresholds for exiting the optimization
  
  for(i in 1:ntrys){
    if(dim(X)[2] == 1){
      fit = tryCatch(Norm1dimEMcycle(X, k, maxiter = 1000, maxtol = 0.0001), function(w){list(ModelLnLikelyhood = -Inf)})
    }else{
      fit = NormMultiDimEMcycle(X, k, maxiter = 1000, maxtol = 0.0001)
    }
    
    if(i == 1){
      bestfit = fit
    }else{
      if(bestfit$ModelLnLikelyhood < fit$ModelLnLikelyhood){
        bestfit = fit
      }
    }
  }
  bestfit
}

ForeachCombineHandler = function(LL1, LL2){
  if(LL1$ModelLnLikelyhood < LL2$ModelLnLikelyhood){
    return(LL2)
  }else{
    return(LL1)
  }
}

NormEMwrapperParallel = function(X, k, maxiter = 1000, maxtol = 0.0001, ntrys = 10, ncores = parallel::detectCores()/2){
  require(foreach)
  # X is vector of length n (n: number of data points)
  # k is number of components/models to by fitted
  # maxiter and maxtol are thresholds for exiting the optimization
  cl = parallel::makeCluster(ncores)
  parallel::clusterExport(cl, c(lsf.str()))
  doParallel::registerDoParallel(cl)
    if(is.vector(X)){
      fit = foreach(i=1:ntrys, .packages = "Rfast", .combine = "ForeachCombineHandler") %dopar% {
        source("EMcalc.R")
        return(Norm1dimEMcycle(X, k, maxiter = 1000, maxtol = 0.0001))
      }
    }else{
      fit = foreach(i=1:ntrys, .packages = c("Rfast", "stats"), .combine = "ForeachCombineHandler") %dopar% {
        source("EMcalc.R")
        return(tryCatch(NormMultiDimEMcycle(X, k, maxiter = 1000, maxtol = 0.0001), error = function(w){list(ModelLnLikelyhood = -Inf)}))
      }
    }
  parallel::stopCluster(cl)
  fit$Group = rowMaxs(fit$Expectaion)
  return(fit)
}


# Wrappers for Binomials --------------------------------------------------

BinomMultiDimEMcycle = function(N, S, k, maxiter, maxtol){
  require(Rfast)
  # X is vector of length n (n: number of data points)
  # k is number of components/models to by fitted
  # maxiter and maxtol are thresholds for exiting the optimization
  
  nDim = dim(N)[2]
  nSamples = dim(N)[1]
  
  # start bby initilazing a potential E
  lnT <- runif(k, 0, 1)
  lnT <- lnT[order(lnT)]
  lnT <- matrix(rep(lnT, nDim), nrow = nDim, byrow = T)
  lnT = log(lnT)
  
  
  empi = runif(k)
  empi = empi/sum(empi)
  
  # setting likelyhood of initial (random) model
  postlnL = NA 
  # starting loop to find optimum
  niter = 0
  con = F
  
  while(!con){
    
    # writing over some variables for comparison
    niter = niter + 1
    priorlnL = postlnL
    
    # Ex step
    lnA = BinomLnA(N, S, lnT, empi)
    lnArs = LnARowSums(lnA) 
    E = exp(lnA - lnArs)
    
    # Max step
    lnT = BinomLnTheta(N, S, E)
    empi = piEst(E)
    
    
    # checking for likelyhood convergence
    postlnL = sum(lnArs)
    if(niter > 5){
      if((postlnL - priorlnL)^2 < maxtol){
        converged = T
        con = T
        return(list(Expectaion = E, BinomPropabilites = exp(lnT), GroupSizePi = empi, ModelLnLikelyhood = postlnL, converged = converged))
      }
    }
    if(niter > maxiter){
      converged = F
      con = T
      return(ModelLnLikelyhood = -Inf, converged = F)
    }
  }
}


BinomEMwrapperParallel = function(N, S, k, maxiter = 1000, maxtol = 0.0001, ntrys = 10, ncores = parallel::detectCores()/2){
  require(foreach)
  # k is number of components/models to by fitted
  # maxiter and maxtol are thresholds for exiting the optimization
  cl = parallel::makeCluster(ncores, setup_strategy = "sequential")
  parallel::clusterExport(cl, c(lsf.str()))
  doParallel::registerDoParallel(cl)
  fit = foreach(i=1:ntrys, .packages = "Rfast", .combine = "ForeachCombineHandler") %dopar% {
      source("/Users/streeck/Desktop/EM-Project/EMcalc.R")
      return(tryCatch(BinomMultiDimEMcycle(N, S, k, maxiter = maxiter, maxtol = maxtol), error = function(w){list(ModelLnLikelyhood = -Inf)}))
    }
  parallel::stopCluster(cl)
  fit$Group = as.character(rowMaxs(as.matrix(fit$Expectaion)))
  return(fit)
}