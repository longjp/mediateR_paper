## the following script roughly reproduces Figure 2 D
## in Gaynor 2018 "Mediation analysis for common binary outcomes"
## the results of this script suggest our code is working
rm(list=ls())
set.seed(05072019)
library(boot)
library(parallel)
library(mediateR)
source("rare_probit_funcs.R")
mc.cores <- 7

## simulation parameters
xx_direct <- 0.4
names(xx_direct) <- "xx1"
mm_direct <- 0.5
names(mm_direct) <- "mm1"
co_direct <- 0.25
names(co_direct) <- "co1"
path <- matrix(1,ncol=1,nrow=1)
rownames(path) <- names(xx_direct)
colnames(path) <- names(mm_direct)
path_model <- path
path_model[1,1] <- 0.5
co_mm <- matrix(0.4,ncol=1,nrow=1)
rownames(co_mm) <- names(co_direct)
colnames(co_mm) <- names(mm_direct)
const_mm <- 0.1
var_mm <- 0.75^2
family <- "binomial"


## loop across k
## loop across number of samples
ks <- seq(from=-3,to=1.5,length.out=100)
N <- 5000
n <- 500

### run simulation
RunSim <- function(ii){
  print(ii)
  res_de <- matrix(0,nrow=3,ncol=N)
  res_ind <- matrix(0,nrow=3,ncol=N)
  k <- ks[ii]
  dat <- SimulatePaper(1e6,k)
  dat2 <- ReformatData(dat)
  sim_params <- list(n=n,path=path,path_model=path_model,co_mm=co_mm,const_mm=const_mm,var_mm=var_mm,
                     xx_direct=xx_direct,mm_direct=mm_direct,co_direct=co_direct,const_direct=k,family=family)
  true_effect_de <- c(dat$prev,ComputeEffectxx(dat2,sim_params,"direct"))
  true_effect_ind <- c(dat$prev,ComputeEffectxx(dat2,sim_params,"indirect"))
  for(jj in 1:N){
    dat <- SimulatePaper(n,k)
    dat2 <- ReformatData(dat)
    fit <- ComputePath(dat2)
    prob_approx <- ComputeProbitApprox(dat)
    res_de[1,jj] <- prob_approx[1]
    res_de[2,jj] <- exp(fit$xx_direct[1])
    res_de[3,jj] <- ComputeEffectxx(dat2,fit,"direct")
    res_ind[1,jj] <- prob_approx[2]
    res_ind[2,jj] <- exp(fit$path_model[1,1]*fit$mm_direct[1])
    res_ind[3,jj] <- ComputeEffectxx(dat2,fit,"indirect")
  }
  return(list(res_ind=rowMeans(res_ind),
              res_de=rowMeans(res_de),
              true_effect_ind=true_effect_ind,
              true_effect_de=true_effect_de))
}

out <- mclapply(1:length(ks),RunSim,mc.cores=mc.cores)

prevs_true_de <- t(vapply(out,function(x){x$true_effect_de},c(0,0)))
colnames(prevs_true_de) <- c("prevalence","direct_effect")
prevs_true_ind <- t(vapply(out,function(x){x$true_effect_ind},c(0,0)))
colnames(prevs_true_ind) <- c("prevalence","indirect_effect")

res_ind <- t(vapply(out,function(x){x$res_ind},rep(0,3)))
colnames(res_ind) <- c("Probit","Rare","Numeric")
res_de <- t(vapply(out,function(x){x$res_de},rep(0,3)))
colnames(res_de) <- c("Probit","Rare","Numeric")

save(res_ind,res_de,prevs_true_ind,prevs_true_de,file="rare_probit.RData")
