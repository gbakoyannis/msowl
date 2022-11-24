# msowl
Estimation of optimal individualized treatment rules with multistate processes via outcome weighted learning


## Description
This repository contains R functions for the estimation of optimal individualized treatment rules with multistate processes via outcome weighted learning. These functions also provide the estimated value function of an individualized treatment rule as well as its (estimated) standard error. The implemented methods are presented in <https://arxiv.org/abs/2204.09785>.

## Main Functions
`itr()` and `V_d()`. Both functions are beta version.

### Dependencies
Package `kernlab`

### Input data
The input data need to be a data frame in the long format containing the variables

* `id`: unique individual identifier.
* `t1`: starting time of the interval in the record.
* `t2`: ending time of the interval in record.
* `s1`: the state of the process at `t1`. The possible values are 1,...,k. 
* `s2`: the state of the process at `t2`. The possible values are 1,...,k.
* `A`: binary treatment. The possible values are -1 and 1.

### Function `itr()`

The function `itr()` estimates an optimal individualized treatment rule with multistate processes. The function has the following arguments:

* `data`: a data.frame in the long format described above.
* `feat`: a specification of the covariates/features to be used for tailoring treatment to individuals. `feat` is a vector containing the variable names for these variables which should be included in `data`.
* `w`: the weight of patient preferences that satisfies 0 <= min(`w`) < max(`w`) <= 1. length(`w`) should be equal to the number of states of the multistate process of interest.
* `tau`: Maximum time to be considered.
* `kernel`: a specification of the form of the decision function. kernel can be set to 'linear' or 'rbf' (radian basis function Gaussian kernel).
* `sigma`: the parameter σ of the Gaussian kernel if kernel = `rbf’.
* `lambda`: the penalty parameter λ.
* `SE`: logical value: if TRUE, the function returns the estimated standard error of the estimated value function of the estimated optimal individualized treatment rule.


### Function `V_d()`

The function `V_d()` computes the estimator the value function of a treatment rule and its standard error. The function has the following arguments:

* `data`: a data.frame in the long format described above.
* `feat`: a specification of the covariates/features to be used for tailoring treatment to individuals. `feat` is a vector containing the variable names for these variables which should be included in `data`.
* `w`: the weight of patient preferences that satisfies 0 <= min(`w`) < max(`w`) <= 1. length(`w`) should be equal to the number of states of the multistate process of interest.
* `tau`: Maximum time to be considered.
* `kernel`: a specification of the form of the decision function. kernel can be set to 'linear' or 'rbf' (radian basis function Gaussian kernel).
* `sigma`: the parameter σ of the Gaussian kernel if kernel = `rbf’.
* `lambda`: the penalty parameter λ.
* `SE`: logical value: if TRUE, the function returns the estimated standard error of the estimated value function of the estimated optimal individualized treatment rule.

## Example

The artificial dataset `example_data.csv` (included in this repository) contains observations from an illness-death process without recovery. The dataset can be obtained as follows
```
> library(foreign)
> data <- read.csv("example_data.csv")
> head(data)
  id        t1         t2 s1 s2          Z1          Z2  A
1  1 0.0000000 0.01380424  1  3 -0.80030228 -0.02345331  1
2  2 0.0000000 0.25340078  1  3 -0.77119721 -0.34741893  1
3  3 0.0000000 0.61103624  1  2 -0.06412938  0.17254223 -1
4  3 0.6110362 0.75620895  2  3 -0.06412938  0.17254223 -1
5  4 0.0000000 0.54577648  1  2  0.14567147  0.66730508 -1
6  4 0.5457765 0.74995655  2  3  0.14567147  0.66730508 -1
```
The variables `Z1` and `Z2` are covariates/features to be used for tailoring treatment to individuals. Estimation of an optimal individual treatment rule for prolonging the time spent in State 2 based on a linear decision function and considering the process up to time τ = 3 can be performed as follows:
```
> fit <- itr(data=data, feat=c("Z1", "Z2"), w = c(0, 1, 0), tau = 3, lambda=1, kernel=`linear’, SE=TRUE)
```
The estimates of the coefficients of the optimal linear decision function can be obtained as follows
```
> fit$beta_opt
```
The estimated value function of the estimated individualized treatment rule can be obtained as follows
```
> fit$V_opt
```
The estimated standard error of the value function of the estimated rule can be obtained as follows
```
> fit$se_V_opt
```
Estimation of the the value function of the latter estimated optimal treatment rule and its standard error using the function `V_d` can be performed as follows:
```
> V_d(data=data, w=c(0, 1, 0), tau=3, dec.fun=fit$fit, feat=c("Z1", "Z2"), SE = TRUE)
```
Estimation of the the value function of the fixed rule that assigns treatment 1 to everyone, along with its standard error, can be performed as follows:
```
> V_d(data=data, w=c(0, 1, 0), tau=3, dec.fun=1, feat=c("Z1", "Z2"), SE = TRUE)
```
Estimation of the the value function of the fixed rule that assigns treatment -1 to everyone, along with its standard error, can be performed as follows:
```
> V_d(data=data, w=c(0, 1, 0), tau=3, dec.fun=-1, feat=c("Z1", "Z2"), SE = TRUE)
```

