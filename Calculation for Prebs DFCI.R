####################################################
#This R file applys each original gene signature and reduced GSs to Prebs (DFCI converted)
#####################################################
library("geneClassifiers")
library(survival)
library(survC1)
#import DFCI expression data 
load("~/downloads/DFCI_TPM_Survival.Rd")
row.names(DFCI_TPM_Gene) <- gsub("\\.\\d+","",row.names(DFCI_TPM_Gene))
#import Prebs dataset
prebs <- as.matrix(read.table("allexpress2.csv",sep = ",",row.names = 1,header = TRUE,stringsAsFactors = FALSE))
colnames(prebs) <- gsub("X","",colnames(prebs))
#check common sample ids between DFCI and Prebs
commonsampleids <- intersect(colnames(prebs),colnames(DFCI_TPM_Gene)) #309
#import DFCI survival information
export <- read.table("dfci_survival.csv",sep = ",",row.names = 1,header = TRUE,stringsAsFactors = FALSE)
#check common sample ids between expression data and survival data
DFCI_common <- intersect(export$clinicalIDs,colnames(prebs))
#save to RData and load
load("~/Documents/Research/DFCI_Prebs_express.RData")
########################################################
#Create a calculation function for calculating risk scores
#Use original gene signatures without removing any probesets
calculation <- function(GeneSignature){
  #load required gene signature
  EMC92 <- data.frame(getWeights(getClassifier(GeneSignature)))
  probes <- rownames(EMC92)
  #Z-transformed expression data
  prebs <- t(prebs[probes,DFCI_common])
  prebs <- t(scale(prebs,center = TRUE, scale = TRUE))
  #set NAs to 0
  prebs[is.na(prebs)] <- 0
  #calculate risk scores
  EMC92_score <- rep(NA,dim(prebs)[2])
  for (i in 1:dim(prebs)[2]){
    express <- prebs[probes,i]
    weight <- EMC92[,1]
    EMC92_score[i] <- sum(weight*express) 
  }
  EMC92_score <- as.data.frame(EMC92_score)
  rownames(EMC92_score) <- colnames(prebs)
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
  dicho <- sort(DFCI_info$score, decreasing = TRUE)[64]
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
calculation("EMC92")
calculation("UAMS70")
calculation("UAMS80")
calculation("UAMS17")
calculation("HM19")
calculation("IFM15")
calculation("MILLENNIUM100")
##########################################
# Calculate risk scores via MRCIX6 
#Z-transformed expression data
G_MRCIX6 <- read.table("MRCIX6_prebs.csv",sep = ",",header = TRUE,stringsAsFactors = FALSE)
prebs <- t(prebs[G_MRCIX6$affy_hg_u133_plus_2,DFCI_common])
prebs <- t(scale(prebs,center = TRUE, scale = TRUE))
#set NAs to 0
prebs[is.na(prebs)] <- 0
EMC92_score <- rep(NA,dim(prebs)[2])
for (i in 1:dim(prebs)[2]){
  if(prebs["203755_at",i]/prebs["216326_s_at",i] >1 | prebs["203213_at",i]/prebs["218034_at",i] >1 | prebs["200608_s_at",i]/prebs["217731_s_at",i] >1){
    EMC92_score[i] <- 1
  }
  else{
    EMC92_score[i] <- 0
  }
}
EMC92_score <- as.data.frame(EMC92_score)
rownames(EMC92_score) <- colnames(prebs)
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
summary(coxph(Surv(DFCI_info$PFS_time,DFCI_info$PFS_indi)~DFCI_info$score)) ###########################################
#Create a calculation function for calculating risk scores
#Use transformed gene signatures with reduced probesets
load("~/Documents/Research/DFCI_Prebs_express.RData")
calculation2 <- function(GeneSignature){
  #load required gene signature
  EMC92symbol <- read.table(GeneSignature,sep = ",",header = TRUE,stringsAsFactors = FALSE)
  probes <- EMC92symbol$probe
  #Z-transformed expression data
  prebs <- t(prebs[probes,DFCI_common])
  prebs <- t(scale(prebs,center = TRUE, scale = TRUE))
  #set NAs to 0
  prebs[is.na(prebs)] <- 0
  #calculate risk scores
  EMC92_score <- rep(NA,dim(prebs)[2])
  for (i in 1:dim(prebs)[2]){
    express <- prebs[probes,i]
    weight <- EMC92symbol$weight
    EMC92_score[i] <- sum(weight*express) 
  }
  EMC92_score <- as.data.frame(EMC92_score)
  rownames(EMC92_score) <- colnames(prebs)
  #check common sample ids between expression data and survival data
  DFCI_common <- intersect(export$clinicalIDs,rownames(EMC92_score))
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
  dicho <- sort(DFCI_info$score, decreasing = TRUE)[64]
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
calculation2("EMC92_mmrf.csv")
calculation2("UAMS70_prebs.csv")
calculation2("UAMS80_mmrf.csv")
calculation2("UAMS17_prebs.csv")
calculation2("HM19_prebs.csv")
calculation2("IFM15_prebs.csv")
calculation2("M100_mmrf.csv")
