library(kernlab)

#Weighted SVM solver
#Code adapted from DTRlearn2 package (source code: https://github.com/ychen178/DTRlearn2/blob/master/R/wsvm_solve.R)

wsvm_solve <-function(X, A, wR, kernel='linear', sigma=0.05, lambda=1, e=1e-7) {
  
  if (kernel=='linear') {
    K = X %*% t(X)
    if (is.vector(X)) K = t(X) %*% X
  }
  else if (kernel=='rbf'){
    rbf = rbfdot(sigma = sigma)
    K = kernelMatrix(rbf, X)
  }
  
  y = A * sign(wR)
  H = y %*% t(y) * K
  
  n = length(A)
  C = 1/(2*n*lambda)
  solution <- tryCatch(ipop(c = rep(-1, n), H = H, A = t(y), b = 0, l = numeric(n), u = C*abs(wR), r = 0), error=function(er) er)

  if ("error" %in% class(solution)) {
    return(list(beta0=NA, beta=NA, fit=NA, probability=NA, treatment=NA, sigma=NA, H=NA, alpha1=NA))
  }
  alpha = primal(solution)
  alpha1 = alpha * y 
  
  if (kernel=='linear'){
    w = t(X) %*% alpha1
    fitted = X %*% w
  } else if (kernel=='rbf'){
    fitted = K %*% alpha1
  }
  rm = y - fitted
  Imid = (alpha < C-e) & (alpha > e)
  rmid = rm[Imid==1]
  if (sum(Imid)>0){
    bias = mean(rmid)
  } else {
    Iup = ((alpha<e)&(A==-sign(wR)))|((alpha>C-e)&(A==sign(wR)))
    Ilow = ((alpha<e)&(A==sign(wR)))|((alpha>C-e)&(A==-sign(wR)))
    rup = rm[Iup]
    rlow = rm[Ilow]
    bias = (min(rup)+max(rlow))/2
  }
  
  
  if (kernel=='linear') {
    model = list(beta0=bias, beta=w, alpha1=alpha1) #, solution=solution) 
    class(model)<-'linear'
  } else if (kernel=='rbf') {
    model = list(beta0=bias, sigma=sigma, Z=X, alpha1=alpha1)
    class(model) = 'rbf'
  }
  return (model)
}

