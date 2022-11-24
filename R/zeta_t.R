
#Utility function for the optimal treatment rule estimation function 

zeta_t <- function(t, data, w=c(0, 1, 0), hazard, S, T){
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