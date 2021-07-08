# NI-margin-conversion-methods

# Introduction

This example file demonstrate the utilities of our (Chen et al. 2021) proposed NI margin conversion methods when applied to both individual study and multiple studies. The example codes perform data generation and NI margin conversions from a HR measure to a RMSTD measure.

# Load Required Functions

First load the data generation functions. Two data generation models are included: (a) Weibull distribution (b) Piece-wise exponential distribution. In order to run the funtions, one need to install the R package 'cpsurvsim'.

## Weibull data generation

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

## Piece-wise Exponential data generation


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



