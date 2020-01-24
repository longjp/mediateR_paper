## compute
rm(list=ls())
set.seed(1234)
library(xtable)
library(ggplot2)
library(survival)
library(ggfortify)
library(graph)
library(mediateR)
load("sim_indep.RData")
source("sim_indep_params.R")
set.seed(1234)

sim_params <- mediateR:::QuickSimMultipleMediator(n=nsi[1],nmm=nmmi[1],var_mm=var_mm,
                                       xx_direct=xx_direct,mm_direct=mm_direct_small,
                                       xx_prob=xx_prob,family)
path_sub <- mediateR:::FindSubgraph(sim_params)
gR <- mediateR:::MakeGraphNELObject(mediateR:::FindSubgraph(sim_params),
                         sim_params$xx_direct[rownames(sim_params$path) %in% rownames(path_sub)],
                         sim_params$mm_direct[colnames(sim_params$path) %in% colnames(path_sub)])
attrs <- list()
attrs$edge <- list()
attrs$edge$fontsize <- 12
edgeAttrs <- mediateR:::MakeedgeAttrs(sim_params$path_model,sim_params$xx_direct,sim_params$mm_direct)
pdf("../ms/figs/sim_indep_graph.pdf",width=8,height=6)
par(mar=c(1,1,1,1))
plot(subGraph(c("x1","m1","m2","m3","m4","m5","y"),gR),edgeAttrs=edgeAttrs,attrs=attrs)
dev.off()


## kaplan-meier plot of 1 simulation
dat <- SimulateData(sim_params)
fit <- survfit(dat$y~as.factor(dat$xx[,1]))
p <- autoplot(fit,xlab="Survival Time (months)",ylab="Survival") +
  coord_cartesian(ylim=c(0,1))
p <- p + labs(colour = "x1") + guides(fill=FALSE)
pdf("../ms/figs/sim_indep_km.pdf",width=6,height=5)
print(p)
dev.off()


## results for large
ylim <- quantile(c(dat_large$indirect,dat_small$indirect),c(0.005,0.995))
p <- ggplot(dat_large, aes(x=n, y=indirect, fill=nmm)) +
  geom_violin() + labs(x="Sample Size", y = "Indirect Effect",fill="# Mediators") +
  geom_hline(yintercept=indirect_t_large, linetype="solid", 
             color = "black", size=0.5,alpha=0.5) + 
  coord_cartesian(ylim=c(ylim[1], ylim[2]))
pdf("../ms/figs/sim_indep_indirect_large.pdf",width=6,height=5)
print(p)
dev.off()
covD <- aggregate(cov~n+nmm+n,data=dat_large,FUN=mean)
rejectD <- aggregate(reject~n+nmm,data=dat_large,FUN=mean)
indirect_large <- data.frame(covD,power=rejectD[,3])
indirect_large <- indirect_large[order(indirect_large$n,indirect_large$nmm),]

## results for small
p <- ggplot(dat_small, aes(x=n, y=indirect, fill=nmm)) +
  geom_violin() + labs(x="Sample Size", y = "Indirect Effect",fill="# Mediators") +
  geom_hline(yintercept=indirect_t_small, linetype="solid", 
             color = "black", size=0.5,alpha=0.5) +
  coord_cartesian(ylim=c(ylim[1], ylim[2]))
pdf("../ms/figs/sim_indep_indirect_small.pdf",width=6,height=5)
print(p)
dev.off()
covD <- aggregate(cov~nmm+n,data=dat_small,FUN=mean)
rejectD <- aggregate(reject~nmm+n,data=dat_small,FUN=mean)
indirect_small <- data.frame(covD,power=rejectD[,3])
indirect_small <- indirect_small[order(indirect_small$n,indirect_small$nmm),]

## merge small and large and output
colnames(indirect_large) <- c("n","No. Med.","CI Cov.","Power")
colnames(indirect_small) <- c("n","No. Med.","CI Cov.","Power")
dat <- cbind(indirect_large,indirect_small[,3:4])
##dat <- dat[(dat$n %in% c(50,200,800)),]
rownames(dat) <- NULL
caption <- "Empirical coverage probabilities and power for simulation with strong mediators and weak mediators.\\label{tab:sim_indep}"

## print results table
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- c(" & & \\multicolumn{2}{c}{Strong Med.} & \\multicolumn{2}{c}{Weak Med.} \\\\\n")
out <- xtable(dat,caption=caption,align="cc|c|cc|cc")
print(out,
      add.to.row=addtorow,
      include.rownames=FALSE,
      hline.after=c(-1,0:nrow(dat)),
      file="../ms/figs/sim_indep.tex")

