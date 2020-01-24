## compute
rm(list=ls())
library(xtable)
library(ggplot2)
library(survival)
library(mediateR)
load("sim_null.RData")
source("sim_null_params.R")

ylim <- quantile(dat_0$indirect,c(0.005,0.995))
p <- ggplot(dat_0, aes(x=n, y=indirect, fill=nmm)) +
  geom_violin() + labs(x="Sample Size", y = "Indirect Effect",fill="# Mediators") +
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", size=0.5,alpha=0.5) + 
  coord_cartesian(ylim=c(ylim[1], ylim[2]))
pdf("../ms/figs/sim_null_indirect_0.pdf",width=6,height=5)
print(p)
dev.off()

# data_0_sub <- dat_0[dat_0$n==50 & dat_0$nmm==5,"indirect",drop=FALSE]
# head(data_0_sub)
# summary(data_0_sub[,"indirect"])
# ggplot(data_0_sub, aes(indirect, fill = indirect)) +
#   geom_density(alpha = 0.3,bw="SJ") + 
#   xlab("Indirect Effect") + ylab("Density")

## output table
indirect_0 <- aggregate(cov~n+nmm+n,data=dat_0,FUN=mean)
indirect_0 <- indirect_0[order(indirect_0$n,indirect_0$nmm),]

colnames(indirect_0) <- c("n","No. Med.","CI Cov.")
dat <- indirect_0
##dat <- dat[(dat$n %in% c(50,200,800)),]
rownames(dat) <- NULL
caption <- paste0("Empirical coverage probabilities for ",
                  100*ci_prob,
                  "\\% confidence intervals in the null simulation (indirect effect=0).\\label{tab:sim_null}")

## print results table
out <- xtable(dat,caption=caption,align="cc|c|c")
print(out,include.rownames=FALSE,hline.after = 0:nrow(dat),file="../ms/figs/sim_null.tex")

