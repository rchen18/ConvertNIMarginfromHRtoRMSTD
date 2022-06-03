# Introduction

This example file demonstrate the utilities of our (Chen et al. 2022) proposed NI margin conversion methods when applied to both individual study and multiple studies. The example codes perform data generation and NI margin conversions from a HR measure to a RMSTD measure.

# Load Required Functions

First load the data generation functions. Two data generation models are included: (a) Weibull distribution (b) Piece-wise exponential distribution. In order to run the funtions, one need to install the R package 'cpsurvsim'.

## Weibull data generation function

``` r

######################################################################
# note that the lambda parameter (scale) is determined after         #
# fixing the shape, time tau, and overall failure rate               #  
######################################################################
getLambdaWeibull <- function(t, # year/months/days
                             k, # a given shape parameter
                             Ft # t-year event probability
) {
  t/((-log(1-Ft, base = exp(1)))^(1/k))
}
rngHCWeibull <- function(
  samp.size.hc, # sample size of the historical control group
  # follow.time.C.min, # minimum follow-up time for the historical control group 
  follow.t.max.hc, # maximum follow-up time for the historical control group
  shape.eve.hc, # shape parameter for the survival time distribution  
  lambda.eve.hc, # scale parameter for the survival time distribution
  lambda.cens.hc, # censoring rate parameter for the historical control group
  # recruit.period.hc, # recruit time for the historical control group
  details=F # whether to display
) {
  # simulate trial entry time
  # recruitment time follows a uniform distribution
  # recruit.t.hc <- runif(samp.size.hc, 0, recruit.period.hc)
  # simulate uncensored survival times
  eve.t.hc <- rweibull(samp.size.hc, shape=shape.eve.hc, scale=lambda.eve.hc)
  # true survival time equals trial entry time plus the uncensored survival time
  # eve.t.hc <- eve.t.hc + recruit.t.hc
  # simulate censoring times (loss to follow-up)
  # censoring time follows an exponential distribution
  cens.t.hc <- rexp(samp.size.hc, rate=lambda.cens.hc)
  # true censoring time equals trial entry time plus the censoring time
  # cens.t.hc <- cens.t.hc + recruit.t.hc 
  # true follow-up time equals entry time plus the maximum follow-up time
  follow.t <- follow.t.max.hc
  # + recruit.t.hc 
  # cut censoring at the maximum follow-up time
  cens.t.hc <- pmin(cens.t.hc, follow.t)
  # calculate observed times
  obs.data.hc <- pmin(eve.t.hc, cens.t.hc)
  # create the censoring status variable
  # status.C = 1 (event observed); status.C = 0 (event censored)
  status.hc <- numeric(samp.size.hc)
  # create right censoring
  for (i in 1:samp.size.hc) {
    if (eve.t.hc[i] < cens.t.hc[i])  {
      status.hc[i]=1
    }
    else {
      status.hc[i]=0
    }
  }
  # create the output data frame
  sim.data.hc <- data.frame(eve.t.hc, cens.t.hc, obs.data.hc, status.hc
                            # ,rep(0,samp.size.hc)
  )
  colnames(sim.data.hc) <-c("true.eve","true.cens","time","status"
                            # ,"trt"
  )
  if (details==T) {
    # realized proportion of events
    print(paste0("prop of events = ",length(which(sim.data.hc$status==1))/nrow(sim.data.hc)))
    # realized proportion of observed censorings
    print(paste0("prop of censored = ",length(which(sim.data.hc$status==0))/nrow(sim.data.hc)))
    # realized proportion of drop-outs
    print(paste0("lost to drop-out = ",length(intersect(which(sim.data.hc$status==0),which(sim.data.hc$time<120)))/nrow(sim.data.hc)))  
    # realized proportion of being censored due to the maximum follow-up time
    print(paste0("lost to follow-up = ",length(intersect(which(sim.data.hc$status==0),which(sim.data.hc$time==120)))/nrow(sim.data.hc)))   
    # detailed results
    res.list <- list(prop.eve = length(which(sim.data.hc$status==1))/nrow(sim.data.hc),
                     prop.cens = length(which(sim.data.hc$status==0))/nrow(sim.data.hc),
                     prop.drop.out = length(intersect(which(sim.data.hc$status==0),which(sim.data.hc$time<120)))/nrow(sim.data.hc),
                     prop.lost.follow = length(intersect(which(sim.data.hc$status==0),which(sim.data.hc$time==120)))/nrow(sim.data.hc),
                     data=sim.data.hc)
    return(res.list)
  }
  else{
    return(sim.data.hc)
  }
}

```

## Piece-wise Exponential (two-change-points) data generation function


``` r
require(cpsurvsim)
rngHCpwExp_cpsurvsim <- function(samp.size.hc, # sample size of the historical control group
                                 follow.t.max.hc, # maximum follow-up time for the historical control group
                                 rate.hc, # a vector of rate parameters of the survival time distribution (exponential)
                                 change.pts.hc, # a vector of change of time points values for the control group
                                 lambda.cens.hc # censoring rate parameter for the historical control group
) {
  require(cpsurvsim)
  samp.data <- exp_cdfsim(n = samp.size.hc, endtime = follow.t.max.hc,
                          theta = rate.hc, tau = change.pts.hc)
  random.drop.outs <- rexp(samp.size.hc, rate=lambda.cens.hc)
  random.drop.outs.cutoff <- pmin(rep(follow.t.max.hc, samp.size.hc), random.drop.outs)
  time.obs <- pmin(samp.data$time, random.drop.outs.cutoff)
  eve.ind <- ifelse(samp.data$time < random.drop.outs.cutoff, 1, 0)
  res.data <- data.frame(time=time.obs, status=eve.ind, true.eve=samp.data$time, true.cens=random.drop.outs.cutoff)
  return(res.data)
}
```

## NI margin conversion function: individual study level with Kaplan-Meier (KM) estimators

This function requires the `survival` R package

input arguments: 

* control.data: the historical active control group individual patient data
* the dataset should have four variables with exact these names:
  1. 'time': the observed event/censoring time
  2. 'status': the event status indicator; 1 means event, 0 means censoring
  3. 'true.eve': the true event time, which might not be observed
  4. 'true.cens': the true censoring time, which might not be observed
* hr.margin: the original NI margin measured in HR
* tau: up to which time point you want to convert the NI margin

``` r

KMcNIMarCov <- function(control.data, hr.margin, tau) {
  # load the necessary dependency library
  require(survival)
  # fit a KM curve (this step is for obtaining survival data info)
  KM.fit1 <- survfit(Surv(time, status) ~ 1, data=control.data)
  # survival info data frame
  KM.fit1.data <- data.frame(S.i = KM.fit1$surv,  # KM est
                             t.i = KM.fit1$time # extract event time info
  )
  ###############################
  # step 1: time length periods #
  ###############################
  # calculate time period length between each recorded time points
  time.len.cont <- numeric(nrow(KM.fit1.data))
  # first time period length is from t=0 to the first event time
  time.len.cont[1] <- KM.fit1.data$t.i[1]
  for (i in 2:length(KM.fit1.data$t.i)) {
    time.len.cont[i] <- KM.fit1.data$t.i[i]-KM.fit1.data$t.i[i-1]
  }
  KM.fit1.data <- cbind(KM.fit1.data, time.len.cont)
  #############################################
  # step 2: determine the converted NI Margin #
  #############################################
  # the following step is extremely important
  # we only need survival data up until the time tau (RMST is a function of time)
  # notice that the number of censoring at the end of the trial would be huge if
  # we don't exclude tau from our analysis
  KM.fit1.data <- KM.fit1.data[1:max(which(KM.fit1.data$t.i<=tau)),]
  # estimated RMST for the NIP-HC group
  mu.NIP.HC <- sum((KM.fit1.data$S.i)^hr.margin * KM.fit1.data$time.len.cont)
  # estimated RMST for the HC group
  mu.HC <- sum(KM.fit1.data$S.i* KM.fit1.data$time.len.cont)
  # the converted margin
  NI.margin.KM.c <- mu.NIP.HC - mu.HC
  return(NI.margin.KM.c)
}
```

## NI margin conversion function: individual study level assuming Weibull model

This function requires the `survival` R package

input arguments: 

* control.data: the historical active control group individual patient data
* the dataset should have four variables with exact these names:
  1. 'time': the observed event/censoring time
  2. 'status': the event status indicator; 1 means event, 0 means censoring
  3. 'true.eve': the true event time, which might not be observed
  4. 'true.cens': the true censoring time, which might not be observed
* hr.margin: the original NI margin measured in HR
* tau: up to which time point you want to convert the NI margin

``` r
weibullNIMarCov <- function(my.data, hr.margin, tau) {
  require(survival)
  # to exclude observations of time 0 (for simulated data)
  my.data <- my.data[which(my.data>0),]
  surv.fit.weibull <- survreg(Surv(my.data$time, my.data$status) ~ 1, dist='weibull')
  # survreg's scale = 1/(rweibull shape)
  est.shape.sim <- 1 / surv.fit.weibull$scale
  # survreg's intercept = log(rweibull scale)
  est.lambda.sim <- exp(as.numeric(surv.fit.weibull$coefficients))
  RMST.hc.sim <- integrate(function(x) exp(-(x/est.lambda.sim)^est.shape.sim), 0, tau)
  RMST.nip.hc.sim <- integrate(function(x) exp(-hr.margin*(x/est.lambda.sim)^est.shape.sim), 0, tau)
  RMST.d.sim <- RMST.nip.hc.sim$value-RMST.hc.sim$value
  return(RMST.d.sim)
}
```

## Variance estimator function of NI margin conversion function with KM estimators for multiple studies

This function calculates with-in study variations for NI margin conversions using multiple studies

``` r
KMVarEst <- function(control.data, hr.margin, tau) {
  #########################################
  # fit a KM curve (this step is for obtaining survival data info)
  KM.fit1 <- survfit(Surv(time, status) ~ 1, data=control.data)
  # survival info data frame
  KM.fit1.data <- data.frame(t.i = KM.fit1$time, # extract event time info
                             S.i = KM.fit1$surv, # survival probabilities
                             d.i = KM.fit1$n.event, # number of event at time t_i
                             n.i = KM.fit1$n.risk # risk size at time t_i
  )
  # we don't exclude tau from our analysis
  # leave out the cutoff point (tau)
  KM.fit1.data <- KM.fit1.data[1:min(max(which(KM.fit1.data$t.i<=tau)),max(which(KM.fit1.data$S.i>0))),]
  # the part inside the big brackets
  var.est.part1 <- (hr.margin^2)*KM.fit1.data$S.i^(2*hr.margin) + (KM.fit1.data$S.i)^2 + 2*hr.margin^2*(KM.fit1.data$S.i)^(hr.margin+1)
  # equation (11) the fraction
  var.est.part2 <- numeric(nrow(KM.fit1.data))
  for (i in 1:nrow(KM.fit1.data)) {
    var.est.part2[i] <- sum(KM.fit1.data$d.i[1:i]/(KM.fit1.data$n.i[1:i]*(KM.fit1.data$n.i[1:i]-KM.fit1.data$d.i[1:i])))                 
  }
  # calculate time period length between each recorded time points
  time.len.cont <- numeric(nrow(KM.fit1.data))
  # first time period length is from t=0 to the first event time
  time.len.cont[1] <- KM.fit1.data$t.i[1]
  for (i in 2:length(KM.fit1.data$t.i)) {
    time.len.cont[i] <- KM.fit1.data$t.i[i]-KM.fit1.data$t.i[i-1]
  }
  # check the time period length result
  # if (sum(time.len.cont) != max(KM.fit1.data$t.i)) {
  #   print("error! sum of all time periods' length should equal the maximum observed time")
  #   break
  # }
  # equation (11) 
  var.est.part3 <- time.len.cont^2
  # var est final (var est up to time tau)
  KMc.rmstd.var <- sum(var.est.part3*var.est.part2*var.est.part1)
  return(KMc.rmstd.var)
}

```

## NI margin conversion function with KM estimators for multiple studies where the variance is estimated by bootstrap

* K denotes the total number of historical active controlled studies used for NI margin conversion 
* n.K is a vector of sample sizes of the K historical studies
* data.K is the data set that contains all K historical studies where each observation has 
* boot.size: how many bootstrap samples are used to calculate the variance
* hr.margin: original NI margin measured in HR 
* tau: up to which time point you want to evaluate the NI margin conversion

``` r
KMcMultiNIMarCovBoot <- function(data.K, K, n.K, hr.margin, tau, boot.size) {
  # weight function container
  KM.C.weights.boot <- numeric(K)
  # est container
  KM.C.est.temp <- numeric(K)
  # estimate RMST using KM-C method for each study
  for (i in 1:K) {
    KM.C.est.temp[i] <- KMcNIMarCov(data.frame(data.K[which(data.K$study==i),]), hr.margin, tau)
  }
  # bootstrapped weighted functions
  # define containers for bootstrapping
  # dim 1 (row): boot est; dim 2 (column): study
  est.boots.mat <- matrix(rep(0,boot.size*K),nrow=boot.size,ncol=K)
  weights.boot <- numeric(K)
  for (i in 1:K) {
    data_K_i <- data.K[which(data.K$study==i),]
    for (j in 1:boot.size) {
      boot.index <- sample(1:nrow(data_K_i), prob=rep(1/nrow(data_K_i), nrow(data_K_i)), replace=T)
      data_K_i_j.boot <- data_K_i[boot.index,]
      est.boots.mat[j,i] <- KMcNIMarCov(data_K_i_j.boot, hr.margin, tau)
    }
    weights.boot[i] <- 1 / var(est.boots.mat[,i])
  }
  # weighted estimator of treatment effect
  KM.C.delta.w.bar <- sum(KM.C.est.temp * weights.boot) / sum(weights.boot)
  KM.C.pop.var.est <- max(0,(sum(weights.boot*(KM.C.est.temp-KM.C.delta.w.bar)^2)-(K-1))/(sum(weights.boot)-(sum(weights.boot^2)/sum(weights.boot))))                       
  # the weight combining within and between study variations
  KM.C.weights.boot.star <- 1/(KM.C.pop.var.est+1/weights.boot)
  #######################
  KM.C.est.RMSTD <- sum(KM.C.weights.boot.star * KM.C.est.temp)/sum(KM.C.weights.boot.star)
  # print(KM.C.est.temp)
  return(KM.C.est.RMSTD)
}

```

## NI margin conversion function with KM estimators for multiple studies where the variance is estimated by the derived equation (Chen et al. 2021)

* K denotes the total number of historical active controlled studies used for NI margin conversion 
* n.K is a vector of sample sizes of the K historical studies
* data.K is the data set that contains all K historical studies where each observation has 
* hr.margin: original NI margin measured in HR 
* tau: up to which time point you want to evaluate the NI margin conversion

``` r
KMcMultiNIMarCovNoBoot <- function(data.K, K, n.K, hr.margin, tau) {
  # weight function container
  KM.C.weights.boot <- numeric(K)
  # est container
  KM.C.est.temp <- numeric(K)
  # estimate RMST using KM-C method for each study
  for (i in 1:K) {
    KM.C.est.temp[i] <- KMcNIMarCov(data.frame(data.K[which(data.K$study==i),]), hr.margin, tau)
  }
  # estimate RMST est var for each study
  weights.noboot <- numeric(K)
  for (i in 1:K) {
    weights.noboot[i] <- 1 / KMVarEst(data.frame(data.K[which(data.K$study==i),]), hr.margin, tau)
  }
  # weighted estimator of treatment effect
  KM.C.delta.w.bar <- sum(KM.C.est.temp * weights.noboot) / sum(weights.noboot)
  KM.C.pop.var.est <- max(0,(sum(weights.noboot*(KM.C.est.temp-KM.C.delta.w.bar)^2)-(K-1))/(sum(weights.noboot)-(sum(weights.noboot^2)/sum(weights.noboot))))                       
  # the weight combining within and between study variations
  KM.C.weights.noboot.star <- 1/(KM.C.pop.var.est+1/weights.noboot)
  #######################
  KM.C.est.RMSTD <- sum(KM.C.weights.noboot.star * KM.C.est.temp)/sum(KM.C.weights.noboot.star)
  # print(KM.C.est.temp)
  return(KM.C.est.RMSTD)
}
```

# Conversion Examples

## Conversion at individual study level

First generate a sample of 200 observaitons that follow from a Weibull(shape=1.5, scale=log(110)) distribution. The generated data are subject to a Type-I right censoring mechanism with a maximum follow-up time of 60 months. A random drop-out mechanism is applied that follows an exponential distribution with lambda=0.001.

``` r
set.seed(123456)
samp.obs.weib <- rngHCWeibull(
    200, # sample size of the historical control group
    # follow.time.C.min, # minimum follow-up time for the historical control group 
    60, # maximum follow-up time for the historical control group
    1.5, # shape parameter for the survival time distribution  
    log(110), # scale parameter for the survival time distribution
    0.001, # censoring rate parameter for the historical control group
    # recruit.period.hc, # recruit time for the historical control group
    details=F # whether to display
)
head(samp.obs.weib)
#####################################
  true.eve true.cens     time status
1 1.743579        60 1.743579      1
2 2.025839        60 2.025839      1
3 4.505390        60 4.505390      1
4 4.930340        60 4.930340      1
5 4.756914        60 4.756914      1
6 6.477632        60 6.477632      1
```

Then convert a HR NI margin of 1.5 at tau=36 months using the KM method.

``` r
KMcNIMarCov(samp.obs.weib, 1.5, 36)
###################################
-1.002105
```
The converted NI margin measured in RMSTD is approximately -1 month.

Then generate a sample of 200 observaitons that follow from a two-change-points Piece-wise Exponential(theta_1=0.0025, r_1=1.3, r_2=0.7, t_1=18, t_2=36) distribution. The generated data are subject to a Type-I right censoring mechanism with a maximum follow-up time of 126 months. A random drop-out mechanism is applied that follows an exponential distribution with lambda=0.001.

``` r
set.seed(12345)
samp.obs.pwe <- rngHCpwExp_cpsurvsim(200, # sample size of the historical control group
                              126, # maximum follow-up time for the historical control group
                              c(0.0025,0.0025*1.3,0.0025*1.3*0.7), # a vector of rate parameters of the survival time distribution (exponential)
                              c(18,36), # a vector of change of time points values for the control group
                              0.001 # censoring rate parameter for the historical control group
)
head(samp.obs.pwe)
#############################
        time status   true.eve true.cens
1 126.000000      0 126.000000 126.00000
2 126.000000      0 126.000000 126.00000
3 126.000000      0 126.000000 126.00000
4   7.379469      1   7.379469 126.00000
5 126.000000      0 126.000000 126.00000
6   9.575256      1   9.575256  30.71664
```

Then convert a HR NI margin of 1.3 at tau=60 months using the KM method.

``` r
KMcNIMarCov(samp.obs.pwe, 1.3, 60)
#################################
-1.440775
```
The converted NI margin measured in RMSTD is approximately -1.4 month.

## Conversion using multiple studies

``` r
# create a sample data set

log.rate.eve.mean <- 0.04
re.var <- 0.5
samp.log.rate.vec <- exp(rnorm(3, log.rate.eve.mean, re.var))
hazard.delta.rates <- c(1.3,0.7)

n.K <- c(550,750,1000)
 set.seed(1234567)
 data_hc_ij_1 <- rngHCpwExp_cpsurvsim(n.K[1], # sample size of the historical control group
                                     follow.t.max.hc=126, # maximum follow-up time for the historical control group
                                     rate.hc=c(samp.log.rate.vec[1],samp.log.rate.vec[1]*hazard.delta.rates[1],samp.log.rate.vec[1]*hazard.delta.rates[1]*hazard.delta.rates[2]), # shape parameter for the survival time distribution
                                     change.pts.hc=c(30,54), # a vector of change of time points values for the control group
                                    lambda.cens.hc=0.01 # censoring rate parameter for the historical control group
)

data_hc_ij_1 <- cbind(data_hc_ij_1,rep(1,nrow(data_hc_ij_1)))
colnames(data_hc_ij_1)[5] <- "study"

data_hc_ij_2 <- rngHCpwExp_cpsurvsim(n.K[2], # sample size of the historical control group
                                     follow.t.max.hc=126, # maximum follow-up time for the historical control group
                                    rate.hc=c(samp.log.rate.vec[2],samp.log.rate.vec[2]*hazard.delta.rates[1],samp.log.rate.vec[2]*hazard.delta.rates[1]*hazard.delta.rates[2]), # shape parameter for the survival time distribution
                                      change.pts.hc=c(30,54), # a vector of change of time points values for the control group
                                      lambda.cens.hc=0.01 # censoring rate parameter for the historical control group
 )

data_hc_ij_2 <- cbind(data_hc_ij_2,rep(2,nrow(data_hc_ij_2)))
colnames(data_hc_ij_2)[5] <- "study"

data_hc_ij_3 <- rngHCpwExp_cpsurvsim(n.K[3], # sample size of the historical control group
                                   follow.t.max.hc=126, # maximum follow-up time for the historical control group
                                    rate.hc=c(samp.log.rate.vec[3],samp.log.rate.vec[3]*hazard.delta.rates[1],samp.log.rate.vec[3]*hazard.delta.rates[1]*hazard.delta.rates[2]), # shape parameter for the survival time distribution
                                      change.pts.hc=c(30,54), # a vector of change of time points values for the control group
                                      lambda.cens.hc=0.01 # censoring rate parameter for the historical control group
 )

data_hc_ij_3 <- cbind(data_hc_ij_3,rep(3,nrow(data_hc_ij_3)))
colnames(data_hc_ij_3)[5] <- "study"

# combine data
data.K <- rbind(data_hc_ij_1,data_hc_ij_2,data_hc_ij_3)
K = 3; tau = 60; hr.margin = 1.2

# run the conversions
hr.margin <- 1.1; tau <- 60; K <- 3; boot.size = 100
KMcMultiNIMarCovNoBoot(data.K, K, n.K, 1.4, 120)
#######################
[1] -0.2459887

# test function
hr.margin <- 1.1; tau <- 60; K <- 3; boot.size = 100
KMcMultiNIMarCovBoot(data.K, K, n.K, 1.4, 120, boot.size)
#######################
[1] -0.2457093
```




