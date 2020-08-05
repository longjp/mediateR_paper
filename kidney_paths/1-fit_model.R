## find protein mediators in small set of pathways
rm(list=ls())
set.seed(1234)
library(mediateR)
library(survival)
library(xtable)
library(ggplot2)
library(ggfortify)

##getwd()
##setwd("~/Desktop/mediateR_paper/kidney_meta")
load("0-get_data.RData")
load("StringV10_ppi_pathway.Rdata")

dat[dat$event_time==0,"event_time"] <- 0.5
y <- Surv(dat$event_time,dat$obs_death)

rna_n <- vapply(strsplit(colnames(rna),"|",fixed=TRUE),function(x){x[1]},c(""))

pten <- c("PTEN")
tca <- c("ACO2","DLAT","DLD","FH","OGDH","PDK1",
         "PDK2","PDK3","PDK4","PDP1","PDP2","SDHB",
         "SDHC","SDHD")
fatty <- c("FASN","ACACA","ACLY")
ampk <- c("PRKAA1","PRKAA2","PRKAB1","PRKAB2","PRKAG1")
pentose <- c("TALDO1","TKT","G6PD","PGLS")
mrna_paths <- list(pten=pten,tca=tca,fatty=fatty,ampk=ampk,pentose=pentose)
mrna_paths_names <- c(pten="PTEN",
                      tca="TCA cycle",
                      fatty="Fatty acid synthesis",
                      ampk="AMPK",
                      pentose="Pentose phosphate")

mrna_pc <- matrix(0,nrow=nrow(rna),ncol=length(mrna_paths))
colnames(mrna_pc) <- names(mrna_paths)
for(ii in 1:length(mrna_paths)){
  ix <- rna_n %in% mrna_paths[[ii]]
  out <- prcomp(rna[,ix,drop=FALSE],scale=TRUE)
  mrna_pc[,ii] <- out$x[,1]
  ## want larger expression to correspond with larger
  ## PC. so if column of rotation matrix sums to -
  ## flip sign of PC
  if(sum(out$rotation[,1]) < 0){
    mrna_pc[,ii] <- -mrna_pc[,ii]
  }
}

## normalize PC scores
for(ii in 1:ncol(mrna_pc)){
  mrna_pc[,ii] <- qqnorm(mrna_pc[,ii],plot.it=FALSE)$x
}



## EDA: plot survival curves, thresholded on cox risk
fit <- coxph(y~mrna_pc)
fit
preds <- predict(fit)
cl_which <- preds < median(preds)
cl <- rep("High",length(cl_which))
cl[cl_which] <- "Low"
cl <- as.factor(cl)
km.fit <- survfit(y~as.factor(cl))


#### RPPA data cleaning

## handle NAs in RPPA data 
## 1) remove frequently missing genes
## 2) impute with median
rppa <- rppa[,colSums(is.na(rppa)) < 26]
colmeds <- apply(rppa,2,function(x){median(x,na.rm=TRUE)})
for(ii in 1:ncol(rppa)){
  rppa[is.na(rppa[,ii]),ii] <- colmeds[ii]
}
sum(is.na(rppa))

## determine RPPA protein pathways
rppa_names <- colnames(rppa)
protein_paths <- unique(ppi$pathway[,1])
ixs <- matrix(FALSE,nrow=length(protein_paths),ncol=length(rppa_names))
rownames(ixs) <- protein_paths
for(ii in 1:nrow(ixs)){
  prot_names <- ppi$pathwaydat[ppi$pathwaydat[,1]==protein_paths[ii],3]
  ix <- rep(FALSE,length(rppa_names))
  for(jj in 1:length(prot_names)){
    ix <- ix | grepl(prot_names[jj],rppa_names)
  }
  ixs[ii,] <- ix
}



nrow(ixs)
range(rowSums(ixs))


for(jj in 1:nrow(ixs)){
  ## compute point estimates
  print(jj)
  mmn <- TRUE
  rmean <- 2000
  dat <- list()
  dat$y <- y
  dat$xx <- mrna_pc
  dat$mm <- scale(rppa[,ixs[jj,]])
  dat$path <- matrix(1,nrow=ncol(dat$xx),ncol=ncol(dat$mm))
  rownames(dat$path) <- colnames(dat$xx)
  colnames(dat$path) <- colnames(dat$mm)
  dat$family <- "cox"
  ## fit model
  fit <- ComputePath(dat,mmn=mmn)
  xp <- quantile(dat$xx[,1],.05)
  xpp <- quantile(dat$xx[,1],.95)
  ## compute indirect effects and confidence intervals
  indirect <- ComputeEffectxx(dat,fit,"indirect",
                              mmn=mmn,rmean=rmean,xp=xp,xpp=xpp)
  direct <- ComputeEffectxx(dat,fit,"direct",
                            mmn=mmn,rmean=rmean,xp=xp,xpp=xpp)
  total <- indirect + direct
  
  B <- 1000
  direct_boot <- matrix(0,ncol=ncol(dat$xx),nrow=B)
  indirect_boot <- matrix(0,ncol=ncol(dat$xx),nrow=B)
  for(ii in 1:B){
    ix <- sample(1:length(dat$y),replace=TRUE)
    dat_boot <- mediateR:::SubsetDat(dat,ix)
    fit_boot <- ComputePath(dat_boot,mmn=mmn)
    indirect_boot[ii,] <- ComputeEffectxx(dat_boot,fit_boot,"indirect",rmean=rmean,mmn=mmn,xp=xp,xpp=xpp)
    direct_boot[ii,] <- ComputeEffectxx(dat_boot,fit_boot,"direct",rmean=rmean,mmn=mmn,xp=xp,xpp=xpp)
  }
  total_boot <- indirect_boot + direct_boot
  
  quants <- apply(indirect_boot,2,function(x){quantile(x,c(0.025,0.975))})
  indirect_CI <- apply(quants,2,
                       function(x){paste0("[",round(x[1]),",",round(x[2]),"]")})
  
  quants <- apply(direct_boot,2,function(x){quantile(x,c(0.025,0.975))})
  direct_CI <- apply(quants,2,
                     function(x){paste0("[",round(x[1]),",",round(x[2]),"]")})
  
  quants <- apply(total_boot,2,function(x){quantile(x,c(0.025,0.975))})
  total_CI <- apply(quants,2,
                    function(x){paste0("[",round(x[1]),",",round(x[2]),"]")})
  
  out <- data.frame(rna=mrna_paths_names,
                    indirect=as.character(round(indirect)),
                    indirect_CI,
                    direct=as.character(round(direct)),
                    direct_CI,
                    total=as.character(round(indirect+direct)),
                    total_CI,
                    stringsAsFactors=FALSE)
  
  ## print results table
  addtorow <- list()
  addtorow$pos <- list(0)
  addtorow$command <- c("Pathway & \\multicolumn{2}{c}{Indirect} & \\multicolumn{2}{c}{Direct}  & \\multicolumn{2}{c}{Total} \\\\\n")
  caption <- paste0("Indirect, Direct, and Total effects (in days) of metabolomic mRNA expression as mediated by ",
                    rownames(ixs)[jj]," pathway protein expression.")
  out <- xtable(out,caption=caption)
  print(out,add.to.row=addtorow,include.colnames=FALSE,
        include.rownames=FALSE,file=paste0("../ms/supp/res",jj,".tex"))
}
