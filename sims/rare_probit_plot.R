rm(list=ls())
load("rare_probit.RData")
library(ggplot2)

## Figure 2a
bias_ind <- res_ind - prevs_true_ind[,2]
bias_ind <- data.frame(rep(colnames(bias_ind),each=nrow(bias_ind)),
                      as.vector(bias_ind),
                      rep(prevs_true_ind[,1],3))
colnames(bias_ind) <- c("Method","Bias","Prevalence")
p1 <- ggplot(bias_ind, aes(x=Prevalence, y=Bias, color=Method,shape=Method)) + geom_point() + 
             theme(legend.position="none") + xlab("Prevalence of Disease (P(y=1))") + 
            ylab("Bias in Effect Estimate")
pdf("../ms/figs/logistic_approx_ind.pdf",width=5,height=4)
plot(p1)
dev.off()

## Figure 2b
bias_de <- res_de - prevs_true_de[,2]
bias_de <- data.frame(rep(colnames(bias_de),each=nrow(bias_de)),
                      as.vector(bias_de),
                      rep(prevs_true_de[,1],3))
colnames(bias_de) <- c("Method","Bias","Prevalence")
p2 <- ggplot(bias_de, aes(x=Prevalence, y=Bias, color=Method,shape=Method)) + geom_point() +
             theme(legend.position = c(0.8, 0.8)) + xlab("Prevalence of Disease (P(y=1))") + 
             ylab("Bias in Effect Estimate")
pdf("../ms/figs/logistic_approx_de.pdf",width=5,height=4)
plot(p2)
dev.off()

