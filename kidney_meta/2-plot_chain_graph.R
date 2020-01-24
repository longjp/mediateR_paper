## plot DAG or fit
rm(list=ls())
library(Rgraphviz)

load("1-fit_model.RData")

## obtain variable names
mrna_names <- colnames(dat$xx)
prot_names <- colnames(dat$mm)
prot_names <- c("AMPK_a","AMPK_pT","ACC_p","ACC","PTEN")
mrna_names[5] <- "pent"

## graph setup
l <- c(mrna_names,prot_names,"y","h")
nmpy <- length(mrna_names) + length(prot_names) + 1
edL <- vector("list",length=length(l))
names(edL) <- l
mp_edge <- (length(mrna_names)+1):nmpy
for(ii in 1:length(mrna_names)){
  edL[[ii]] <- list(edges=mp_edge,rep(1,length(mp_edge)))
}
prot_edges <- (length(mrna_names)+1):(nmpy-1)
for(ii in (length(mrna_names)+1):(nmpy-1)){
  prot_edges_temp <- prot_edges[prot_edges!=ii]
  edL[[ii]] <- list(edges=c(prot_edges_temp,nmpy),rep(1,length(prot_edges_temp)))
}
edL[[nmpy+1]] <- list(edges=1:length(mrna_names),rep(1,length(mrna_names)))
gR <- graphNEL(nodes=l,edgeL=edL,edgemode="directed")


### create agopen object
nodes <- buildNodeList(gR)
edges <- buildEdgeList(gR)
## make protein edges undirected
for(ii in 1:(length(prot_names)-1)){
    for(jj in (ii+1):length(prot_names)){
      ename <- paste0(prot_names[ii],"~",prot_names[jj])
      edges[[ename]]@attrs$arrowsize <- "0"
    }
}
for(ii in prot_names){
  nodes[[ii]]@attrs$color <- "orange"
}
vv <- agopen(name="foo",nodes=nodes,edges=edges,edgeMode="directed")

ConvertLWD <- function(x){
  return(exp(3*x))
}

## use coefficients for line width
## direct
names(fit$xx_direct) <- mrna_names
names(fit$mm_direct) <- prot_names
direct_fits <- c(fit$xx_direct,fit$mm_direct)
direct_fits_color <-  as.character(2*(direct_fits > 0) + 2)
names(direct_fits_color) <- names(direct_fits)
direct_fits_color[direct_fits_color=="2"] <- "red"
direct_fits_color[direct_fits_color=="4"] <- "blue"
direct_fits_lwd <- ConvertLWD(direct_fits)
## path
path_model <- fit$path_model
colnames(path_model) <- prot_names
rownames(path_model) <- mrna_names
path_model_color <- 2*(path_model < 0) + 2
path_model_color
path_model_color[path_model_color==2] <- "red"
path_model_color[path_model_color=="4"] <- "blue"
path_model_color
path_model_lwd  <- ConvertLWD(path_model)

for(ii in 1:length(vv@AgEdge)){
  if(vv@AgEdge[[ii]]@head=="y"){
    vv@AgEdge[[ii]]@lwd <- direct_fits_lwd[vv@AgEdge[[ii]]@tail]
    vv@AgEdge[[ii]]@color <- direct_fits_color[vv@AgEdge[[ii]]@tail]
  }
  if(vv@AgEdge[[ii]]@head %in% prot_names & vv@AgEdge[[ii]]@tail %in% prot_names){
    vv@AgEdge[[ii]]@color <- "#00000030"
    vv@AgEdge[[ii]]@lwd <- 0.5
  }
  if(vv@AgEdge[[ii]]@head %in% prot_names & vv@AgEdge[[ii]]@tail %in% mrna_names){
    ##vv@AgEdge[[ii]]@color <- "blue"
    vv@AgEdge[[ii]]@lwd <- path_model_lwd[vv@AgEdge[[ii]]@tail,vv@AgEdge[[ii]]@head]
    vv@AgEdge[[ii]]@color <- path_model_color[vv@AgEdge[[ii]]@tail,vv@AgEdge[[ii]]@head]
  }
}


pdf("../ms/figs/TCGA_dag.pdf",width=8,height=6)
plot(vv)
dev.off()
