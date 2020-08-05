rm(list=ls())
set.seed(01042019)
library(parallel)
library(mediateR)
source("sim_null_params.R")

## parameters input to simulate
sim_set <- cbind(rep(rep(nsi,each=N),length(nmmi)),
                 rep(nmmi,each=N*length(nsi)))
colnames(sim_set) <- c("n","nmm")

## compute indirect effects and CI
RunSim <- function(ii,mm_direct,indirect_t){
  print(paste0(ii,"/",nrow(sim_set)))
  sim_params <- QuickSimMultipleMediator(n=sim_set[ii,1],nmm=sim_set[ii,2],var_mm=var_mm,
                                         xx_direct=xx_direct,mm_direct=mm_direct,
                                         xx_prob=xx_prob,family=family)
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
  covP <- 1*(0 > lb & 0 < ub)
  return(c(indirect,covP))
}

sim_params <- QuickSimMultipleMediator(n=n_true,nmm=5,var_mm=var_mm,
                                       xx_direct=xx_direct,mm_direct=mm_direct_0,
                                       xx_prob=xx_prob,family=family)
dat <- SimulateData(sim_params)
rmean <- max(as.matrix(dat$y)[,1])
fit <- ComputePath(dat,mmn=mmn)
indirect_t_0 <- ComputeEffectxx(dat,fit,"indirect",mmn=mmn,rmean=rmean)
total_t_0 <- ComputeEffectxx(dat,fit,"total",mmn=mmn,rmean=rmean)
out <- mclapply(1:nrow(sim_set),RunSim,
                mm_direct=mm_direct_0,indirect_t=indirect_t_0,
                mc.cores=mc.cores)

dat <- cbind(sim_set,matrix(unlist(out),ncol=2,byrow=TRUE))
colnames(dat) <- c("n","nmm","indirect","cov")
dat <- as.data.frame(dat)
dat$n <- as.factor(dat$n)
dat$nmm <- as.factor(dat$nmm)
dat_0 <- dat


save(dat_0,
     indirect_t_0,total_t_0,
     file="sim_null.RData")

