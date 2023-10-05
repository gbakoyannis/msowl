

# Multistate outcome weighted learning

library(survival)
library(kernlab)

#setwd("H:/msowl_hiv/programs")
#setwd("/Users/giorgosbakoyannis/Documents/MSOWL_HIV")

#dat <- read.csv("example_data.csv")

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
    model <- NA
    class(model) <- "ipop function error"
  } else {
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
      model = list(beta0=bias, beta=w, alpha1=alpha1)
      class(model)<-'linear'
    } else if (kernel=='rbf') {
      model = list(beta0=bias, sigma=sigma, Z=X, alpha1=alpha1)
      class(model) = 'rbf'
    }
  }
  return (model)
}



#Utility function for the optimal treatment rule estimation function
zeta_t <- function(t, data, w=c(0, 1, 0), fit, S, T){
  hazard <- basehaz(fit, centered=FALSE)
  #Get the value that corresponds to the maximum t that is <= "time"
  element <- function(t){
    max(1,max(1:length(hazard$time)*(hazard$time<=t)))
  }
  events <- list()
  for(j in 1:length(S)){
    if(!(S[j] %in% T)){
      events[[j]] <- (data$s1==S[j] & data$t1<=t & data$t2>t)
    } else {
      events[[j]] <- (data$s2==S[j] & data$t2<=t)
    }
  }
  event <- do.call(cbind, events)%*%w

  Tt <- (1*!(data$s2 %in% T))*t +
    (1*(data$s2 %in% T))*mapply(min, data$t2, t) #T^t

  elements <- sapply(Tt, element)
  G_t <- exp(-hazard$hazard[elements])

  res <- event/G_t
  return(res)
}




ms <- function(t1, t2, s1, s2, a){
  #if (missing(id))
  #  stop("Missing id argument")
  if (missing(t1))
    stop("Missing values in `t1'")
  if (missing(t2))
    stop("Missing values in `t2'")
  if (missing(s1))
    stop("Missing values in `s1'")
  if (missing(s2))
    stop("Missing values in `s2'")
  if (max(t1 >= t2)==1)
    stop("t1 >= t2")
  if (max(!(unique(a) %in% c(-1, 1)))==1)
    stop("Treatmetn var `a' is not equal to -1 or 1")
  msdat <- cbind(t1 = t1, t2 = t2, s1 = s1, s2 = s2, a = a)
  class(msdat) <- "msdat"
  msdat
}


#Estimator of the value function of an ITR
#function(formula, id, w = c(0, 1, 0), tau = 3, data, lambda = 1){
msowl.val <- function(formula, id, w = c(0, 1, 0), tau = 3, data,
                      rule, fixed.rules = FALSE, SE = TRUE, nboot = 100,
                      trim = NULL){
  if(class(rule)!="msowl")
    stop("Option `rule' should be the output of an msowl() run")
  if(!(fixed.rules %in% c(TRUE, FALSE)))
    stop("Option `fixed.rules' should be TRUE or FALSE")
  if(!(SE %in% c(TRUE, FALSE)))
    stop("Option `SE' should be TRUE or FALSE")
  if(!is.null(trim)){
    if(trim < 0 | trim > 0.5){
      stop("Option `trim' should be between >=0 and <= 0.5")
    }
  }
  V_d <- function(data){
    msout <- model.frame(formula, data)[[1]]
    covs <- model.matrix(formula, data)[, -1]
    #data <- as.data.frame(cbind(id = data[,id], msout, covs))
    data <- cbind.data.frame(id = data[,id],
                             as.data.frame(cbind(msout, covs)))
    feat <- colnames(data)[(ncol(msout) + 2):ncol(data)]
    rm(msout, covs)
    data <- data[order(data[,id], data[,"t1"]),]

    # form <- as.formula(paste("Surv(t1)", covs, sep = " ~ "))
    ## state space
    S <- with(data, sort(unique(c(s1, s2))))

    ## absorbing state subspace
    T <- S[!(S %in% sort(unique(data$s1)))]

    if(length(S)!=length(w)){
      stop("length(w) not equal to the total number of states")
    }

    if(min(w)<0 | max(w)>1 | sum(w)>=length(w)){
      stop("min(w)<0 or max(w)>1 or sum(w)>=length(w)")
    }

    #Calculate reward
    #Estimate censoring distribution
    c.delta <- 1*(data$s1==data$s2 & data$t2<tau)
    form <- as.formula("Surv(time=t1, time2=t2, event=c.delta,
                             type='counting') ~ 1")
    fit <- coxph(form, data=data)
    tt <- sort(unique(c(survfit(fit, se.fit = FALSE)$time,
                        data[data$s1!=data$s2, "t2"])))
    tt <- c(0,tt[tt<=tau])
    event <- matrix(sapply(tt, zeta_t, data = data, w = w, fit = fit, S = S, T = T),
                    byrow = TRUE, nrow = length(tt))
    dm <- diff(c(tt, tau), lag=1)
    xi <- colSums(event*dm)
    xi <- matrix(tapply(X=xi, INDEX = data$id, FUN=sum))
    dat <- data[data$t1==0,]
    dat <- dat[order(dat[, id]),]

    pi_n <- mean(dat$a==1)
    Wt <- xi/(pi_n*(dat$a==1) + (1 - pi_n)*(dat$a==-1))

    A <- dat$a

    form <- as.formula(paste("~", strsplit(as.character(rule$call), "~", fixed = TRUE))[[3]])
    Z <- model.matrix(form, dat)
    f_n <- Z%*%rule$beta_opt
    V_dn <- mean(Wt*(A*f_n>=0))

    if(fixed.rules){
      f <- 1
      V_1 <- mean(Wt*(A*f>=0))
      f <- -1
      V_m1 <- mean(Wt*(A*f>=0))
      res <- c(V_dn, V_1, V_m1)
      names(res) <- c("V(dn)", "V(1)", "V(-1)")
    } else {
      res <- V_dn
      names(res) <- "V(dn)"
    }
    return(res)
  }
  est <- V_d(data = data)
  if(SE){
    pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                         max = nboot,  # Maximum value of the progress bar
                         style = 3,    # Progress bar style (also available style = 1 and style = 2)
                         width = 50,   # Progress bar width. Defaults to getOption("width")
                         char = "=")   # Character used to create the bar

    if(!is.character(data$id)){
      data$id <- as.character(data$id)
    }
    clusters <- sort(unique(data[,id]))
    bres <- matrix(NA, nrow=nboot, ncol=(fixed.rules*3 + !fixed.rules*1))
    for(b in 1:nboot){
      index <- sample(1:length(clusters), length(clusters), replace=TRUE)
      cl <- clusters[index]
      bdat <- NULL
      for(j in 1:length(cl)){
        aa <- data[data[,id] %in% cl[j],]
        reps <- table(cl[1:j])[names(table(cl[1:j]))==cl[j]] - 1
        if(reps > 0){
          aa[,id] <- paste(aa[1,id], 1*reps, sep = ".")
        }
        bdat <- rbind(bdat, aa)
      }
      bVal <- try(V_d(data = bdat), silent = TRUE)
      if(class(bVal) != "try-error"){
        bres[b,] <- bVal
      }
      setTxtProgressBar(pb, b)
    }
    close(pb) # Close the connection

    if(fixed.rules){
      nfail <- sum(is.na(bres[,1]))
      SE <- apply(bres, 2, sd, na.rm = TRUE)

      res <- cbind(est, SE, est - qnorm(0.975)*SE, est + qnorm(0.975)*SE)
      colnames(res) <- c("estimate", "SE", "ll", "ul")
      res <- round(res, 3)

      C <- rbind(c(1, -1, 0),
                 c(1, 0, -1))
      est <- C%*%est
      SE <- sqrt(diag(C%*%var(bres, na.rm = TRUE)%*%t(C)))
      res2 <- cbind(est, est - qnorm(0.975)*SE, est + qnorm(0.975)*SE)
      colnames(res2) <- c("estimate", "ll", "ul")
      rownames(res2) <- c("V(dn) - V(1)", "V(dn) - V(-1)")
      res2 <- round(res2, 3)

      res <- list(res, res2, n.failed = nfail)
    } else {
      nfail <- sum(is.na(bres[,1]))
      SE <- sd(bres, na.rm = TRUE)
      res <- c(est, SE, est - qnorm(0.975)*SE, est + qnorm(0.975)*SE)
      names(res) <- c("estimate", "SE", "ll", "ul")
      res <- list(round(res, 3), n.failed = nfail)
    }
  } else {
    res <- est
  }

  return(res)
}


msowl <- function(formula, id, w = c(0, 1, 0), tau = 3, 
                  data, lambda = 1, jackknife = FALSE, trim = NULL){
  msout <- model.frame(formula, data)[[1]]
  covs <- model.matrix(formula, data)[, -1]
  data <- cbind.data.frame(id = data[,id],
                           as.data.frame(cbind(msout, covs)))
  feat <- colnames(data)[(ncol(msout) + 2):ncol(data)]
  rm(msout, covs)
  data <- data[order(data[,id], data[,"t1"]),]

  ## state space
  S <- with(data, sort(unique(c(s1, s2))))

  ## absorbing state subspace
  T <- S[!(S %in% sort(unique(data$s1)))]

  if(length(S)!=length(w)){
    stop("length(w) not equal to the total number of states")
  }

  if(min(w)<0 | max(w)>1 | sum(w)>=length(w)){
    stop("min(w)<0 or max(w)>1 or sum(w)>=length(w)")
  }

  #Calculate reward
  #Estimate censoring distribution
  c.delta <- 1*(data$s1==data$s2 & data$t2<tau)
  form <- as.formula("Surv(time=t1, time2=t2, event=c.delta,
                              type='counting') ~ 1")
  fit <- coxph(form, data=data)
  tt <- sort(unique(c(survfit(fit, se.fit = FALSE)$time,
                      data[data$s1!=data$s2, "t2"])))
  tt <- c(0,tt[tt<=tau])
  event <- matrix(sapply(tt, zeta_t, data = data, w = w, fit = fit, S = S, T = T),
                  byrow = TRUE, nrow = length(tt))
  dm <- diff(c(tt, tau), lag=1)
  xi <- colSums(event*dm)
  xi <- matrix(tapply(X=xi, INDEX = data$id, FUN=sum))
  dat <- data[data$t1==0,]
  dat <- dat[order(dat[, id]),]

  pi_n <- mean(dat$a==1)
  Wt <- xi/(pi_n*(dat$a==1) + (1 - pi_n)*(dat$a==-1))
 
  X <- as.matrix(dat[,feat])
  A <- dat$a

  fit <- wsvm_solve(X=X, A=A, wR=Wt, kernel='linear', lambda=lambda)

  if(unlist(gregexpr('error', class(fit)))[1]==-1){
    beta_opt <- c(fit$beta0, fit$beta)
    names(beta_opt) <- c("(constant)", feat)
    f_n <- fit$beta0 + X%*%fit$beta
    V_dn <- mean(Wt*(A*f_n>=0))

    if(isTRUE(jackknife)){
      ids <- sort(unique(dat$id))
      f <- rep(NA, times=nrow(dat))
      pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                           max = length(ids),  # Maximum value of the progress bar
                           style = 3,    # Progress bar style (also available style = 1 and style = 2)
                           width = 50,   # Progress bar width. Defaults to getOption("width")
                           char = "=")   # Character used to create the bar
      for(i in ids){
        fit1 <- wsvm_solve(X=X[dat$id!=i,], A=A[dat$id!=i], wR=Wt[dat$id!=i],
                           kernel='linear', lambda=lambda)
        if(unlist(gregexpr('error', class(fit1)))[1]==-1){
          f[dat$id==i] <- fit1$beta0 + X[dat$id==i,]%*%fit1$beta
        }
        setTxtProgressBar(pb, i)
      }
      close(pb) # Close the connection
      #V_jack <- mean(Wt*(A*f>=0))
      ###Sean code added. Issue with wsvm_solve giving NAs with some jackkniffed samples (add NA to output)
      vecVal=Wt*(A*f>=0)
      V_jack <- mean(vecVal,na.rm=TRUE)
      SD_jack=sd(vecVal,na.rm=TRUE)
      nMiss_jack=sum(is.na(vecVal))
      res <- list(beta_opt = beta_opt, Value = V_dn, Value_jack = V_jack,SD_jack=SD_jack,nMiss_jack=nMiss_jack, fit = fit, call = formula)
    } else {
      res <- list(beta_opt = beta_opt, Value = V_dn, fit = fit, call = formula)
    }
    class(res) <- "msowl"
  } else {
    res <- NA
    class(res) <- "error"
  }
  return(res)
}

#Select the optimal lambda through leave-one-out cross-validation
select.lambda <- function(formula, id, w, tau, data, 
                       lambdas = c(0.01, 0.1, 0.5, 1, 5, 10, 50, 100)){
  V_jack <- rep(NA,length(lambdas))
  SD_jack <- rep(NA,length(lambdas))
  nMiss_jack <- rep(NA,length(lambdas))
  for(i in 1:length(lambdas)){
    print(paste("lambda", lambdas[i], sep = " = "))
    fit <- msowl(formula, id = id, w = w, tau = tau, data = data,
                 lambda = lambdas[i], jackknife = TRUE)
    if(unlist(gregexpr('error', class(fit)))[1]==-1){
      V_jack[i] <- fit$Value_jack
      SD_jack[i] <- fit$SD_jack
      nMiss_jack <- fit$nMiss_jack
    }
  }
  #Should pick lambda better with standard deviation as well
  lambda <- min(lambdas[V_jack==max(V_jack, na.rm = TRUE)], na.rm = TRUE)
  res <- list(lambda = lambda, details = cbind(lambdas, V_jack,SD_jack,nMiss_jack))
  return(res)
}


#### Examples ####

#select.lambda(ms(t1 = t1, t2 = t2, s1 = s1, s2 = s2, a = A) ~ Z1 + Z2,
#                 id = "id", w = c(0, 1, 0), tau = 3, data = dat,
#                 lambdas = c(0.01, 0.1, 0.5, 1, 5, 10, 50, 100))

#fit <- msowl(ms(t1 = t1, t2 = t2, s1 = s1, s2 = s2, a = A) ~ Z1 + Z2,
#             id = "id", w = c(0, 1, 0), tau = 3, data = dat,
#             lambda = 1)

#Value <- msowl.val(ms(t1 = t1, t2 = t2, s1 = s1, s2 = s2, a = A) ~ Z1 + Z2,
#                   id = "id", w = c(0, 1, 0), tau = 3, rule = fit, fixed.rules = TRUE, data = dat)








