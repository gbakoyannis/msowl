

#Estimator of the value function of an ITR

V_d <- function(data, w=c(0, 1, 0), tau=3, dec.fun, feat=NULL, SE=TRUE){
  
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
  
  A <- dat$A
  g <- A*mean(A==1) + (1 - A)/2
  
  if(class(dec.fun)=="linear"){
    Z <- as.matrix(dat[,feat])
    f <- dec.fun$beta0 + Z%*%dec.fun$beta
  } else if (class(dec.fun)=="rbf"){
    rbf <- rbfdot(sigma = dec.fun$sigma)
    K <- kernelMatrix(rbf, dec.fun$Z, as.matrix(dat[,feat]))
    f <- dec.fun$beta0 + K%*%dec.fun$alpha1
  } else {
    f <- dec.fun
  }
  
  V_n <- mean((xi/g)*(A*f>=0))
  
  if(class(dec.fun)!="rbf" & isTRUE(SE)){
    
    ## SE estimation
    #Keep last record for each id
    myid.uni <- unique(data$id)
    a<-length(myid.uni)
    last <- c()
    for (i in 1:a) {
      temp<-subset(data, id==myid.uni[i])
      if (dim(temp)[1] > 1) {
        last.temp<-temp[dim(temp)[1],]
      }
      else {
        last.temp<-temp
      }
      last<-rbind(last, last.temp)
    }
    rm(temp, last.temp)
    last <- last[order(last$id),]
    
    last$d <- 1*(last$s1==last$s2 & last$t2<tau)
    
    ar_t<-function(t){
      mean(last$t2>=t)
    }
    dM_t<-function(t){
      last$d*(last$t2==t)-(last$t2>=t)*mean(last$d*(last$t2==t))/ar_t(t)
    }
    uft <- sort(unique(last[last$d==1, "t2"]))
    ar<-sapply(uft, ar_t, simplify=TRUE)
    dM<-sapply(uft, dM_t, simplify=TRUE) #Each column corresponds to a t
    
    #Infinite-dimensional parameter
    gamma_t <- function(t){
      if(sum(uft<=t)>1){
        gamma <- colSums(t(dM[,(uft<=t)])/ar[uft<=t])
      } else if(sum(uft<=t)==1){
        gamma <- dM[,(uft<=t)]/ar[uft<=t]
      } else {
        gamma <- rep(0, times=nrow(dM))
      }
      return(gamma)
    }
    if(length(tt)<=1000){
      z_t <- lapply(tt, zeta_h_t, data=data, w=w, hazard=Haz, S=S, T=T, h=gamma_t, dm=dm)
      xi2 <- t(Reduce('+', z_t))
      rm(z_t)
    } else {
      #In this case z_t is a huge object and RAM is not sufficient
      xi2 <-  0*zeta_h_t(t=1, data=data, w=w, hazard=Haz, S=S, T=T, h=gamma_t, dm=dm)
      times <- split(tt, cut(tt, quantile(tt, prob = 0:10 / 10, names = FALSE), include = TRUE))
      for(l in 1:10){
        z_t <- lapply(times[[l]], zeta_h_t, data=data, w=w, hazard=Haz, S=S, T=T, h=gamma_t, dm=dm)
        xi2 <- xi2 + Reduce('+', z_t)
        rm(z_t)
      }
      xi2 <- t(xi2)
    }
    xi2 <- aggregate(xi2,by=list(data$id), sum)[,(2:(ncol(xi2)+1))]
    eta <- colMeans((xi2/g)*(A*f>=0))
    
    W <- mean(A*(xi/g^2)*(A*f>=0))
    psi <- (xi/g)*(A*f>=0) - V_n + eta - W*((A==1) - mean(A==1))
    se <- sqrt(mean(psi^2)/length(psi))
  } else {
    if(isTRUE(SE)){
      warning("No standard error estimation with RBF kernel decision function")
    }
    se <- NA
  }
  res <- list(V_n=V_n, se=se)
  return(res)
}
