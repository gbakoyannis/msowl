## Bakoyiannis 2023 Optimal ITR Example
# https://github.com/gbakoyannis/msowl

source("./R/itr.R")
source("./R/V_d.R")
source("./R/wsvm_solve.R")
source("./R/zeta_h_t.R")
source("./R/zeta_t.R")


options(error=traceback) 

# The artificial dataset `example_data.csv` (included in this repository) contains observations from an illness-death process without recovery. 


# The dataset can be obtained as follows
library(foreign)
data <- read.csv("./data/example_data.csv")
print("DATA")
head(data)

# id        t1         t2 s1 s2          Z1          Z2  A
# 1  1 0.0000000 0.01380424  1  3 -0.80030228 -0.02345331  1
# 2  2 0.0000000 0.25340078  1  3 -0.77119721 -0.34741893  1
# 3  3 0.0000000 0.61103624  1  2 -0.06412938  0.17254223 -1
# 4  3 0.6110362 0.75620895  2  3 -0.06412938  0.17254223 -1
# 5  4 0.0000000 0.54577648  1  2  0.14567147  0.66730508 -1
# 6  4 0.5457765 0.74995655  2  3  0.14567147  0.66730508 -1

# The variables `Z1` and `Z2` are covariates/features to be used for tailoring treatment to individuals. 
# Estimation of an optimal individual treatment rule for prolonging the time spent in State 2 based on a linear decision function and considering the process up to time τ = 3 can be performed as follows:

fit <- itr(
  data=data, 
  feat=c("Z1", "Z2"), 
  w = c(0, 1, 0), 
  tau = 3, 
  lambda=1, 
  kernel="linear", 
  SE=TRUE
)

print("The estimates of the coefficients of the optimal linear decision function:")
fit$beta_opt

print("The estimated value function of the estimated individualized treatment rule:")
fit$V_opt

print("The estimated standard error of the value function of the estimated rule:")
fit$se_V_opt

res = list()

# Estimation of the the value function of the latter estimated optimal treatment rule and its standard error using the function `V_d` can be performed as follows:
res$opt = V_d(data=data, w=c(0, 1, 0), tau=3, dec.fun=fit$fit, feat=c("Z1", "Z2"), SE = TRUE)
print(paste("Optimal Value of Indivualized Treatment:", res$opt[1], "+/-", res$opt[2]))

# Estimation of the the value function of the fixed rule that assigns treatment 1 to everyone, along with its standard error, can be performed as follows:
res$pos = V_d(data=data, w=c(0, 1, 0), tau=3, dec.fun=1, feat=c("Z1", "Z2"), SE = TRUE)

print(paste("Value of Treatment 1:", res$pos[1], "+/-", res$pos[2]))

# Estimation of the the value function of the fixed rule that assigns treatment -1 to everyone, along with its standard error, can be performed as follows:
res$neg = V_d(data=data, w=c(0, 1, 0), tau=3, dec.fun=-1, feat=c("Z1", "Z2"), SE = TRUE)

print(paste("Value of Treatment -1:", res$neg[1], "+/-", res$neg[2]))
