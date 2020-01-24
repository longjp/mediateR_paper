## acquire TCGA data
rm(list=ls())
source("TCGA-Assembler/Module_A.R")
source("TCGA-Assembler/Module_B.R")
library(KEGGREST)

sCancer <- "KIRC"


###### STEP 1: GET rna
sPath12 <- "tcga_kirc"
path_geneExp <- DownloadRNASeqData(cancerType = sCancer,
                                   assayPlatform = "gene.normalized_RNAseq",
                                   saveFolderName = sPath12)
rna <- read.table(path_geneExp,
                  header=TRUE,row.names=1)
dim(rna)
head(rna[,1])
rna <- t(as.matrix(rna))


#' Download RPPA protein expression data
path_protein_RPPA <-
  DownloadRPPAData(cancerType = sCancer,
                   assayPlatform = "protein_RPPA",
                   saveFolderName = sPath12)
path_protein_RPPA
rppa <- read.table(path_protein_RPPA,
                  header=TRUE,sep="\t",stringsAsFactors=FALSE,row.names=1)
dim(rppa)
colnames(rppa)
rppa <- t(as.matrix(rppa))
rppa[1:3,1:3]





###### STEP 2: GET CLINICAL
### downloaded from : https://portal.gdc.cancer.gov/exploration?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.primary_site%22%2C%22value%22%3A%5B%22Kidney%22%5D%7D%7D%5D%7D
### clinical -> .tsv
dat <- read.table("clinical.tsv",sep="\t",header=TRUE,
                  stringsAsFactors=FALSE,na.strings=c("NA","--"))
dat <- dat[dat$project_id=="TCGA-KIRC",]

## rule for survival
## if have death time, event time=death time AND not censored
## if no death time but time of last follow up, event time=time of last follow up AND censored
## individuals who do not have a death time or data of last follow up are removed
event_time <- dat$days_to_death
obs_death <- 1*(!is.na(event_time))
## everyone is either dead or has last follow up time
sum(is.na(dat$days_to_death))
sum(is.na(dat$days_to_death) & is.na(dat$days_to_last_follow_up))
event_time[is.na(event_time)] <- dat$days_to_last_follow_up[is.na(event_time)]

dat$event_time <- event_time
dat$obs_death <- obs_death
sum(is.na(dat$event_time))

## make kaplan meier survival curve for everyone
dim(dat)
library(survival)
fit <- survfit(Surv(dat$event_time/365,dat$obs_death)~1,
               conf.type="log-log")
par(mar=c(4,4,1,1))
plot(fit,mark.time=TRUE)
abline(h=.5,col='grey')
abline(h=.4,col='grey')


###### STEP 3: SELECT / MERGE DATA

## some samples are of normal tissue. for example
## in code "TCGA-CZ-5985-11A-01R-1672-07" the "11A" indicates normal.
## we remove these rna expressions
## see: https://www.biostars.org/p/313063/ for description of bar code meanings
rownames(rna) <- gsub(".","-",rownames(rna),fixed=TRUE)
head(rownames(rna))
## only keep "01A" samples for rna
temp <- strsplit(rownames(rna),"-",fixed="TRUE")
temp <- vapply(temp,function(x){x[[4]]},c(""))
table(temp)
ix <- temp=="01A"
rna <- rna[ix,]
dim(rna)




rownames(rppa) <- gsub(".","-",rownames(rppa),fixed=TRUE)
head(rownames(rppa))
## only keep "01A" samples for rppa
temp <- strsplit(rownames(rppa),"-",fixed="TRUE")
temp <- vapply(temp,function(x){x[[4]]},c(""))
table(temp)
ix <- temp=="01A"
rppa <- rppa[ix,]
dim(rppa)





## nearly all snp in rna and all rna has clinical
## all rna data has clinical, 98% of clinical has rna
pnames <- substr(rownames(rppa),1,12)
rnames <- substr(rownames(rna),1,12)
length(pnames)
length(unique(pnames))
length(rnames)
length(unique(rnames))
mean(pnames %in% rnames)
mean(rnames %in% pnames)
head(rnames)
head(dat$submitter_id)
length(rnames)
length(dat$submitter_id)
length(unique(dat$submitter_id))
length(unique(rnames))
mean(dat$submitter_id %in% rnames)
mean(rnames %in% dat$submitter_id)

## only keep subjects for whom we have rna,dna,clinical
common_samples <- intersect(rnames,intersect(pnames,dat$submitter_id))
length(common_samples)
dat <- dat[dat$submitter_id %in% common_samples,]
rppa <- rppa[pnames %in% common_samples,]
rna <- rna[rnames %in% common_samples,]
dim(dat)
dim(rppa)
dim(rna)

## order the clinical, rna, snp data
rna <- rna[order(rownames(rna)),]
rppa <- rppa[order(rownames(rppa)),]
dat <- dat[order(dat$submitter_id),]
identical(substr(rownames(rna),1,12),dat$submitter_id)
identical(substr(rownames(rppa),1,12),dat$submitter_id)


## supplement data with kegg pathway information
rna_gene_names <- vapply(strsplit(colnames(rna),"|",fixed=TRUE),function(x){x[2]},c(""))
kegg_name <- paste0("hsa:",rna_gene_names)
name_to_path <- keggLink("pathway", "hsa")

path_names <- names(table(name_to_path))
path_to_name <- vector("list",length=length(path_names))
names(path_to_name) <- path_names
for(ii in 1:length(path_names)){
  path_to_name[[ii]] <- names(name_to_path)[name_to_path==path_names[ii]]
}

save(rna,dat,rppa,name_to_path,path_to_name,file="0-get_data.RData")
