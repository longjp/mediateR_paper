rm(list=ls())
set.seed(01042019)
library(parallel)
library(mediateR)
source("sim_indep_params.R")

## parameters input to simulate
sim_set <- cbind(rep(rep(nsi,each=N),length(nmmi)),
                 rep(nmmi,each=N*length(nsi)))
colnames(sim_set) <- c("n","nmm")

## compute indirect effects and CI
RunSim <- function(ii,mm_direct,indirect_t){
  print(paste0(ii,"/",nrow(sim_set)))
  sim_params <- QuickSimMultipleMediator(n=sim_set[ii,1],nmm=sim_set[ii,2],var_mm=var_mm,
                                         xx_direct=xx_direct,mm_direct=mm_direct,
                                         xx_prob=xx_prob,family)
  dat <- SimulateData(sim_params)
  fit <- ComputePath(dat,mmn=mmn)
  indirect <- ComputeEffectxx(dat,fit,"indirect",mmn=mmn,rmean=rmean)
  indirect_boot <- rep(NA_real_,B)
  for(jj in 1:B){
    ix <- sample(1:sim_params$n,replace=TRUE)
    dat_boot <- mediateR:::SubsetDat(dat,ix)
    fit_boot <- ComputePath(dat_boot,mmn=mmn)
    indirect_boot[jj] <- ComputeEffectxx(dat_boot,fit_boot,"indirect",mmn=mmn,rmean=rmean)
  }
  lb <- as.numeric(quantile(indirect_boot,(1-ci_prob)/2,na.rm=TRUE))
  ub <- as.numeric(quantile(indirect_boot,1-(1-ci_prob)/2,na.rm=TRUE))
  covP <- 1*(indirect_t > lb & indirect_t < ub)
  reject <- 1*(ub < 0 | lb > 0)
  return(c(indirect,covP,reject))
}


## SIMULATE WITH LARGE INDIRECT EFFECT
sim_params <- QuickSimMultipleMediator(n=n_true,nmm=5,var_mm=var_mm,
                                       xx_direct=xx_direct,mm_direct=mm_direct_large,
                                       xx_prob=xx_prob,family)
dat <- SimulateData(sim_params)
rmean <- max(as.matrix(dat$y)[,1])
fit <- ComputePath(dat,mmn=mmn)
indirect_t_large <- ComputeEffectxx(dat,fit,"indirect",mmn=mmn,rmean=rmean)
total_t_large <- ComputeEffectxx(dat,fit,"total",mmn=mmn,rmean=rmean)
out <- mclapply(1:nrow(sim_set),RunSim,
                mm_direct=mm_direct_large,indirect_t=indirect_t_large,
                mc.cores=mc.cores)
dat <- cbind(sim_set,matrix(unlist(out),ncol=3,byrow=TRUE))
colnames(dat) <- c("n","nmm","indirect","cov","reject")
dat <- as.data.frame(dat)
dat$n <- as.factor(dat$n)
dat$nmm <- as.factor(dat$nmm)
dat_large <- dat


## SIMULATE WITH SMALL INDIRECT EFFECT
sim_params <- QuickSimMultipleMediator(n=n_true,nmm=5,var_mm=var_mm,
                                       xx_direct=xx_direct,mm_direct=mm_direct_small,
                                       xx_prob=xx_prob,family)
dat <- SimulateData(sim_params)
rmean <- max(as.matrix(dat$y)[,1])
fit <- ComputePath(dat,mmn=mmn)
indirect_t_small <- ComputeEffectxx(dat,fit,"indirect",mmn=mmn,rmean=rmean)
total_t_small <- ComputeEffectxx(dat,fit,"total",mmn=mmn,rmean=rmean)
out <- mclapply(1:nrow(sim_set),RunSim,
                mm_direct=mm_direct_small,indirect_t=indirect_t_small,
                mc.cores=mc.cores)
dat <- cbind(sim_set,matrix(unlist(out),ncol=3,byrow=TRUE))
colnames(dat) <- c("n","nmm","indirect","cov","reject")
dat <- as.data.frame(dat)
dat$n <- as.factor(dat$n)
dat$nmm <- as.factor(dat$nmm)
dat_small <- dat

save(dat_large,dat_small,
     indirect_t_large,total_t_large,
     indirect_t_small,indirect_t_small,
     file="sim_indep.RData")