
#Utility function for SE estimation

zeta_h_t <- function(t, data, w=c(0, 1, 0), hazard, S, T, h, dm){
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
  event <- as.numeric(event)
  
  Tt <- (1*!(data$s2 %in% T))*t + 
    (1*(data$s2 %in% T))*mapply(min, data$t2, t) #T^t
  
  elements <- sapply(Tt, element)
  h_t <- sapply(Tt, h, simplify=TRUE)
  #h_t <- do.call(cbind, h_t)
  G_t <- exp(-hazard$hazard[elements])
  res <- (event/G_t)*h_t*dm[element(t)]
  return(res)
}
