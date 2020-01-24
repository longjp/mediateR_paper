### parameters for simulation
mc.cores <- 28
family <- "cox"
var_mm <- 1
xx_direct <- 0.5
xx_prob <- 0.5
R2 <- xx_prob*(1-xx_prob) / (xx_prob*(1-xx_prob)  + var_mm)

## should consider network mediators?
mmn <- TRUE

## constants, alter for speed
n_true <- 10000 ## sample size for approximating true ratio
N <- 500 ## number of simulation runs
nsi <- 2^(0:4)*50 ## sample sizes for simulation
nmmi <- c(5,10,20) ## number of mediators for simulation
B <- 1000 ## number of bootstrap samples
ci_prob <- 0.95

## indirect effect sizes
mm_direct_large <- 0.2
mm_direct_small <- 0.1