##############################################
#This R file applys each transferred GSs to DFCI RNAseq datasets
##############################################
library("geneClassifiers")
library(survival)
library(survC1)
#import DFCI expression data and survival information
load("~/downloads/DFCI_TPM_Survival.Rd")
row.names(DFCI_TPM_Gene) <- gsub("\\.\\d+","",row.names(DFCI_TPM_Gene))
#import survival information
export <- read.table("dfci_survival.csv",sep = ",",row.names = 1,header = TRUE)
#sample ids with survial information
DFCI_common <- intersect(export$clinicalIDs,colnames(DFCI_TPM_Gene))
#save to RData and load
load("~/Documents/Research/DFCI_express.RData")
#####################################
#create a calculation function for calculating risk scores
calculation <- function(GeneSignature){
  #load required gene signature
  EMC92symbol <- read.table(GeneSignature,sep = ",",header = TRUE,stringsAsFactors = FALSE)
  #Z-transformed expression data
  DFCI_TPM_Gene <- t(DFCI_TPM_Gene[EMC92symbol$id,DFCI_common])
  DFCI_TPM_Gene <- t(scale(DFCI_TPM_Gene,center = TRUE, scale = TRUE))
  #calculate risk scores
  EMC92_score <- rep(NA,327)
  for (i in 1:327){
    express <- DFCI_TPM_Gene[EMC92symbol$id,i]
    weight <- EMC92symbol$weight
    EMC92_score[i] <- sum(weight*express) 
  }
  EMC92_score <- as.data.frame(EMC92_score)
  rownames(EMC92_score) <- colnames(DFCI_TPM_Gene)
  #check common sample ids between expression data and survival data
  #DFCI_common <- intersect(export$clinicalIDs,rownames(EMC92_score))#327
  DFCI_info <- export[export$clinicalIDs %in% DFCI_common,]
  #assign risk scores
  DFCI_info$score <- EMC92_score[DFCI_common,]
  inputdata <- apply(DFCI_info,2,as.numeric)
  set.seed(123)
  ###############################survival
  #calculate C-statistics
  OS <- Est.Cval(inputdata[,c(11,9,12)],2*365.25,nofit = TRUE)[["Dhat"]] 
  PFS <- Est.Cval(inputdata[,c(10,8,12)],2*365.25,nofit = TRUE)[["Dhat"]]  
  ################################
  #base on top 21.7% risk score
  dicho <- sort(DFCI_info$score, decreasing = TRUE)[71]
  #assign top 21.7% risk scores to 1
  DFCI_info$percent2 <- as.numeric(DFCI_info$score > dicho)
  #create survival objects for calculating HR
  DFCI_info$SurvObj <- with(DFCI_info, Surv(OS_time, OS_indi))
  DFCI_info$SurvObj2 <- with(DFCI_info, Surv(PFS_time, PFS_indi))
  #Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
  OS_HR <- summary(coxph(Surv(DFCI_info$OS_time,DFCI_info$OS_indi)~DFCI_info$percent2)) 
  PFS_HR <- summary(coxph(Surv(DFCI_info$PFS_time,DFCI_info$PFS_indi)~DFCI_info$percent2)) 
  return(list(OS, PFS, OS_HR, PFS_HR))
}
calculation("EMC92_mmrf.csv")
calculation("UAMS70_prebs.csv")
calculation("UAMS80_mmrf.csv")
calculation("UAMS17_prebs.csv")
calculation("HM19_prebs.csv")
calculation("IFM15_prebs.csv")
calculation("M100_mmrf.csv")
#######################################
# Calculate risk score via MRCIX6 
#Z-transformed expression data
G_MRCIX6 <- read.table("MRCIX6_prebs.csv",sep = ",",header = TRUE,stringsAsFactors = FALSE)
DFCI_TPM_Gene <- t(DFCI_TPM_Gene[G_MRCIX6$ensembl_gene_id,DFCI_common])
DFCI_TPM_Gene <- t(scale(DFCI_TPM_Gene,center = TRUE, scale = TRUE))
EMC92_score <- rep(NA,327)
for (i in 1:327){
  if(DFCI_TPM_Gene["ENSG00000156970",i]/DFCI_TPM_Gene["ENSG00000171720",i] >1 | DFCI_TPM_Gene["ENSG00000170312",i]/DFCI_TPM_Gene["ENSG00000214253",i] >1 | DFCI_TPM_Gene["ENSG00000164754",i]/DFCI_TPM_Gene["ENSG00000136156",i] >1){
    EMC92_score[i] <- 1
  }
  else{
    EMC92_score[i] <- 0
  }
}
EMC92_score <- as.data.frame(EMC92_score)
rownames(EMC92_score) <- colnames(DFCI_TPM_Gene)
#check common sample ids between expression data and survival data
#DFCI_common <- intersect(export$clinicalIDs,rownames(EMC92_score))#327
DFCI_info <- export[export$clinicalIDs %in% DFCI_common,]
#assign risk scores
DFCI_info$score <- EMC92_score[DFCI_common,]
inputdata <- apply(DFCI_info,2,as.numeric)
set.seed(123)
###############################survival
#calculate C-statistics
Est.Cval(inputdata[,c(11,9,12)],2*365.25,nofit = TRUE)[["Dhat"]] 
Est.Cval(inputdata[,c(10,8,12)],2*365.25,nofit = TRUE)[["Dhat"]]  
################################
#create survival objects for calculating HR
DFCI_info$SurvObj <- with(DFCI_info, Surv(OS_time, OS_indi))
DFCI_info$SurvObj2 <- with(DFCI_info, Surv(PFS_time, PFS_indi))
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(DFCI_info$OS_time,DFCI_info$OS_indi)~DFCI_info$score)) 
summary(coxph(Surv(DFCI_info$PFS_time,DFCI_info$PFS_indi)~DFCI_info$score)) 
