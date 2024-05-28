#library(compiler)


fun.obj <- function(vecL, vecL.ref, Delta, lambda){
  L <- matrix(vecL, nrow=nrow(Delta))
  beta <- diag(L %*% t(L))
  ones <- rep(1, length(beta))
  C <- Delta - (beta %*% t(ones) + ones %*% t(beta) -2*L %*% t(L))
  obj <- norm(C, 'f')^2 + lambda*sum(abs(vecL - vecL.ref))
  return(obj)
}

fun.grad.smooth <- function(vecL, Delta){
  L <- matrix(vecL, nrow=nrow(Delta))
  beta <- diag(L %*% t(L))
  ones <- rep(1, length(beta))
  C <- Delta - (beta %*% t(ones) + ones %*% t(beta) -2*L %*% t(L))
  grad <- 8*(C - diag(as.vector(C %*% ones))) %*% L
  grad <- c(as.vector(grad))
  return(grad)
}

fun.stepsize <- function(vecL, gradL, vecL.ref, lambda, eta, Delta){
  
  f <- function(x){
    vecL.new <- prox(vecL = vecL - x * gradL, vecL.ref = vecL.ref, lambda=lambda, eta=x)
    obj <- fun.obj(vecL = vecL.new, vecL.ref = vecL.ref, 
                   Delta = Delta, lambda = lambda)
    return(obj)
  }
  #size <- optimize(f, lower=0, upper=2*eta, tol=1e-10)$minimum
  size <- optimize(f, lower=0, upper=10*eta)$minimum
  return(size)
}


treelax <- function(Delta = NULL,
                    L.ref = NULL,
                    X = NULL,
                    popsize = NULL,
                    L.init = NULL,
                    lambda,
                    maxiter = 1e4,
                    rel.tol=1e-6,
                    eta.init = 100
                    ){
  
  if(is.null(Delta) & is.null(X)) stop("Either Delta or X must be provided")

  if(is.null(Delta)){
    P <- compute_P_from_X(X, popsize)
    Delta <- compute_Delta_from_P(P)
  }
  if(is.null(L.ref)){
    out_treemix <- estimate_L_from_X_treemix(X, popsize, m1=FALSE, m2=FALSE)
    L.ref <- out_treemix$L_m0[,-1] # remove the D0 (root)
    
    fun.L.scaling <- function(d){
      return (norm(Delta - compute_Delta_from_L(L.ref %*% diag(d)), 'f')^2)
    }
    out.L.scaling <- optim(par = rep(1, ncol(L.ref)),
                           fn = fun.L.scaling,
                           method = 'L-BFGS-B',
                           lower = 0,
                           control = list(factr = 1e5))
    for (j in 1:ncol(L.ref)){ L.ref[,j] <- L.ref[,j] * out.L.scaling$par[j] }
    out_include_L.ref <- TRUE
  }else{
    out_include_L.ref <- FALSE
  }

  vecL.ref <- as.vector(L.ref)
  if(is.null(L.init)){
    vecL <- vecL.ref
  }else{
    vecL <- as.vector(L.init)
  }
  
  eta <- eta.init
  
  obj <- fun.obj(vecL = vecL, vecL.ref = vecL.ref, Delta = Delta, lambda = lambda)
  objvals <- c(obj)
  
  for (iter in 1:maxiter){
    gradL <- fun.grad.smooth(vecL = vecL, Delta = Delta)
    if(sum(gradL^2)==0) break
    eta <- fun.stepsize(vecL = vecL, gradL = gradL, vecL.ref = vecL.ref, lambda = lambda, eta = eta, Delta = Delta)
    vecL <- prox(vecL = vecL - eta * gradL, vecL.ref = vecL.ref, lambda = lambda, eta = eta)
    
    obj <- fun.obj(vecL = vecL, vecL.ref = vecL.ref, Delta = Delta, lambda = lambda)
    objvals <- c(objvals, obj)
    if(iter>10 & (objvals[length(objvals)]/objvals[length(objvals)-1] > 1-rel.tol)) break
  }
  if(iter==maxiter){warning('maximum iteration reached')}
  
  L <- matrix(vecL, nrow=nrow(Delta))
  rownames(L) <- rownames(L.ref); colnames(L) <- colnames(L.ref)
  
  
  if(out_include_L.ref){
    out.list <- list(L = L, objvals = objvals, L.ref = L.ref)
  }else{
    out.list <- list(L = L, objvals = objvals)
  }
  return(out.list)
}


treelax.cv <- function(X,
                       popsize,
                       lambda.grid, 
                       nfold=5, maxiter=1e4, rel.tol=1e-6, eta.init=100){
  
  mat.mse <- matrix(0, nrow=length(lambda.grid), ncol=nfold)
  blocksize <- floor(nrow(X)/nfold)
  
  for (f in 1:nfold){
    index_from <- (f-1)*blocksize + 1
    index_to <- f*blocksize
    
    X.train <- X[-(index_from:index_to), ]
    X.test  <- X[index_from:index_to, ]
    
    P.train <- compute_P_from_X(X.train, popsize)
    P.test  <- compute_P_from_X(X.test, popsize)
    
    Delta.train <- compute_Delta_from_P(P.train)
    Delta.test  <- compute_Delta_from_P(P.test)
    
    # 'estimate_L_from_X_treemix' currently in 'drift_1kgp_Feb2024.Rmd' 
    # need to make it faster (only estimate tree)
    out_treemix <- estimate_L_from_X_treemix(X.train, popsize, m1=FALSE, m2=FALSE)
    L.ref <- out_treemix$L_m0[,-1] # remove the D0 (root)
    
    # adjust the scaling of drift sizes so that the approximation is tighter
    fun.L.scaling <- function(d){
      return (norm(Delta.train - compute_Delta_from_L(L.ref %*% diag(d)), 'f')^2)
    }
    out.L.scaling <- optim(par = rep(1, ncol(L.ref)),
                           fn = fun.L.scaling,
                           method = 'L-BFGS-B',
                           lower = 0)
    for (j in 1:ncol(L.ref)){ L.ref[,j] <- L.ref[,j] * out.L.scaling$par[j] }
    
    L.init <- NULL
    for (i in 1:length(lambda.grid)){
      lambda <- lambda.grid[i]
      out <- treelax(Delta = Delta.train, 
                     L.ref = L.ref,
                     L.init = L.init, # use warm start
                     lambda = lambda, 
                     maxiter = maxiter, rel.tol = rel.tol, eta.init = eta.init)
      Delta.hat <- compute_Delta_from_L(out$L)
      mat.mse[i,f] <- mean((Delta.hat - Delta.test)^2)
      L.init <- out$L
      
      #print(c(length(out$objvals), mat.mse[i,f]))
      
    }
  }
  return(mat.mse)
}



# util functions --------------------------

# soft thresholding operator
softthres <- function(z, lambda){
  return(sign(z) * pmax(abs(z) - lambda, 0))
}


# proximal operator
prox <- function(vecL, vecL.ref, lambda, eta){
  vecL <- softthres(vecL - vecL.ref, eta*lambda) + vecL.ref
  vecL <- pmax(vecL, 0)
  return(vecL)
}


# prepare P matrix (allele frequency) from X matrix (allele count)
compute_P_from_X <- function(X, popsize){
  
  P <- matrix(0, nrow=nrow(X), ncol=ncol(X))
  for (k in 1:ncol(X)){
    P[,k] <- X[,k]/(2*popsize[k])
  }
  colnames(P) <- colnames(X)
  return(P)
}


# prepare Delta matrix (distance matrix) from P matrix (allele frequency)
# P is a S-by-K matrix, where S is the number of SNPs and K is the number of populations
compute_Delta_from_P <- function(P){
  Delta <- dist(x = t(P), method = 'euclidean', diag = TRUE, upper = TRUE)^2 /nrow(P)
  Delta <- as.matrix(Delta)
  return(Delta)
}

# prepare Delta matrix (distance matrix) from L matrix (drift membership)
# L is a K-by-J matrix, where K is the number of populations and J is the number of drifts
compute_Delta_from_L <- function(L){
  Delta <- dist(x = L, method = 'euclidean', diag = TRUE, upper = TRUE)^2
  Delta <- as.matrix(Delta)
  return(Delta)
}

