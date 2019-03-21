##############################################
#This R file applys each transferred GSs to MMRF RNAseq datasets
##############################################
library("geneClassifiers")
library(survival)
library(survC1)
# #import MMRF expression data 
# MMRF <- as.matrix(read.table("MMRF_CoMMpass_IA10c_E74GTF_Salmon_Gene_TPM.txt",sep = "\t",row.names = 1,header = TRUE,stringsAsFactors = FALSE ))
# #change the sample names
# colnames(MMRF) <- gsub("_1_BM","",colnames(MMRF))
# #import survival information
# MMRF_infp <- read.table("~/desktop/MMRF.csv",sep = ",",header = TRUE,stringsAsFactors = FALSE)
# MMRF_infp <- na.omit(MMRF_infp)
# #sample ids with survial information
# MMRF_common <- intersect(MMRF_infp$public_id,colnames(MMRF)) #647
#save to RData and load
load("~/Documents/Research/MMRF_express.RData")
#####################################
#create a calculation function for calculating risk scores
calculation <- function(GeneSignature){
  #load required gene signature
  EMC92symbol <- read.table(GeneSignature,sep = ",",header = TRUE,stringsAsFactors = FALSE)
  #Z-transformed expression data
  MMRF <- t(MMRF[EMC92symbol$id,MMRF_common])
  MMRF <- t(scale(MMRF,center = TRUE, scale = TRUE))
  #calculate risk scores
  EMC92_score <- rep(NA,647)
  for (i in 1:647){
    express <- MMRF[EMC92symbol$id,i]
    weight <- EMC92symbol$weight
    EMC92_score[i] <- sum(weight*express) 
  }
  EMC92_score <- as.data.frame(EMC92_score)
  rownames(EMC92_score) <- colnames(MMRF)
  #check common sample ids between expression data and survival data
  #MMRF_common <- intersect(MMRF_infp$public_id,rownames(EMC92_score)) #647
  MMRF_infp2 <- MMRF_infp[MMRF_infp$public_id %in% MMRF_common,]
  #assign risk scores
  MMRF_infp2$score <- EMC92_score[MMRF_common,]
  inputdata <- apply(MMRF_infp2,2,as.numeric)
  set.seed(123)
  ###############################survival
  #calculate C-statistics
  OS <- Est.Cval(inputdata[,c(3,2,6)],2*365.25,nofit = TRUE)[["Dhat"]] 
  PFS <- Est.Cval(inputdata[,c(5,4,6)],2*365.25,nofit = TRUE)[["Dhat"]]  
  ################################
  #base on top 21.7% risk score
  dicho <- sort(MMRF_infp2$score, decreasing = TRUE)[141]
  #assign top 21.7% risk scores to 1
  MMRF_infp2$percent2 <- as.numeric(MMRF_infp2$score > dicho)
  #create survival objects for calculating HR
  MMRF_infp2$SurvObj <- with(MMRF_infp2, Surv(ttcos, censos))
  MMRF_infp2$SurvObj2 <- with(MMRF_infp2, Surv(ttcpfs, censpfs))
  #Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
  OS_HR <- summary(coxph(Surv(MMRF_infp2$ttcos,MMRF_infp2$censos)~MMRF_infp2$percent2)) 
  PFS_HR <- summary(coxph(Surv(MMRF_infp2$ttcpfs,MMRF_infp2$censpfs)~MMRF_infp2$percent2)) 
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
MMRF <- t(MMRF[G_MRCIX6$ensembl_gene_id,MMRF_common])
MMRF <- t(scale(MMRF,center = TRUE, scale = TRUE))
EMC92_score <- rep(NA,647)
for (i in 1:647){
  if(MMRF["ENSG00000156970",i]/MMRF["ENSG00000171720",i] >1 | MMRF["ENSG00000170312",i]/MMRF["ENSG00000214253",i] >1 | MMRF["ENSG00000164754",i]/MMRF["ENSG00000136156",i] >1){
    EMC92_score[i] <- 1
  }
  else{
    EMC92_score[i] <- 0
  }
}
EMC92_score <- as.data.frame(EMC92_score)
rownames(EMC92_score) <- colnames(MMRF)
#check common sample ids between expression data and survival data
#MMRF_common <- intersect(MMRF_infp$public_id,rownames(EMC92_score)) #647
MMRF_infp2 <- MMRF_infp[MMRF_infp$public_id %in% MMRF_common,]
#assign risk scores
MMRF_infp2$score <- EMC92_score[MMRF_common,]
inputdata <- apply(MMRF_infp2,2,as.numeric)
set.seed(123)
###############################survival
#calculate C-statistics
Est.Cval(inputdata[,c(3,2,6)],2*365.25,nofit = TRUE)[["Dhat"]] 
Est.Cval(inputdata[,c(5,4,6)],2*365.25,nofit = TRUE)[["Dhat"]]  
################################
#create survival objects for calculating HR
MMRF_infp2$SurvObj <- with(MMRF_infp2, Surv(ttcos, censos))
MMRF_infp2$SurvObj2 <- with(MMRF_infp2, Surv(ttcpfs, censpfs))
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(MMRF_infp2$ttcos,MMRF_infp2$censos)~MMRF_infp2$score)) 
summary(coxph(Surv(MMRF_infp2$ttcpfs,MMRF_infp2$censpfs)~MMRF_infp2$score)) 

