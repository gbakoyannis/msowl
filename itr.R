# Estimation of an optimal individualized treatment rule

itr <- function(data, feat, w=c(0, 1, 0), tau=3, kernel='linear', 
                lambda=1, sigma=1, SE=TRUE){
  
  if(!(kernel %in% c('linear', 'rbf'))){
    stop("Only 'linear' and 'rbf' (radial basis function or Gaussian) kernels are supported")
  }
  ## state space
  S <- with(data, sort(unique(c(s1, s2))))
  
  ## absorbing state subspace
  T <- S[!(S %in% sort(unique(data$s1)))]
  
  if(length(S)!=length(w)){
    stop("length(w) not equal to the total number of states")
  }
  
  if(max(data$t1>data$t2)==1){
    stop("t1 > t2")
  }
  
  if(min(w)<=0 | max(w)>=1 | sum(w)>=length(w)){
    stop("min(w)<0 or max(w)>1 or sum(w)>=length(w)")
  }
  
  #Calculate reward
  #Estimate censoring distribution 
  c.delta <- 1*(data$s1==data$s2 & data$t2<tau)
  fit <- coxph(Surv(time=t1, time2=t2, event=c.delta,
                    type="counting") ~ 1, data=data)
  Haz <- basehaz(fit, centered=FALSE)
  tt <- sort(unique(c(Haz$time, data[data$s1!=data$s2,"t2"])))
  tt <- c(0,tt[tt<=tau])
  event <- matrix(sapply(tt, zeta_t, data=data, w=w, hazard=Haz, S=S, T=T), byrow=TRUE, nrow=length(tt))
  dm <- diff(c(tt, tau), lag=1)
  xi <- colSums(event*dm)
  xi <- matrix(tapply(X=xi, INDEX = data$id, FUN=sum))
  dat <- data[data$t1==0,]
  dat <- dat[order(dat$id),]
  
  pi_n <- mean(dat$A==1)
  Wt <- xi/(pi_n*(dat$A==1) + (1 - pi_n)*(dat$A==-1))
  X <- as.matrix(dat[,feat])
  A <- dat$A
  
  if(kernel=='linear'){
    fit <- wsvm_solve(X=X, A=A, wR=Wt, kernel='linear', lambda=lambda)
    
    if(!is.na(fit$beta0)){
      beta_opt <- c(fit$beta0, fit$beta)
      V_opt <- V_d(data=data, w=w, tau=tau, dec.fun=fit, feat=feat, SE=SE)
      res <- list(beta_opt=beta_opt, V_opt=V_opt$V_n, se_V_opt=V_opt$se, fit=fit)
    } else {
      res <- NA
      class(res) <- "error"
    }
  } else {
    fit <- wsvm_solve(X=X, A=A, wR=Wt, kernel='rbf', sigma=sigma, lambda=lambda)
    if(!is.na(fit$beta0)){
      V_opt <- V_d(data=data, w=w, tau=tau, dec.fun=fit, feat=feat, SE=FALSE)
      res <- list(V_opt=V_opt, fit=fit)
    } else {
      res <- NA
      class(res) <- "error"
    }
  }
  return(res)
}


