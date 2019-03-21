#This R file applys each original gene signature and reduced GSs to Prebs (MMRF converted)
# MMRF_info <- read.csv("MMRF sample information.csv",sep = ",",header = TRUE,stringsAsFactors = FALSE)
# MMRF_info$Sample_Name <- gsub("_1_BM.+$","",MMRF_info$Sample_Name)
# MMRF <- read.table("MMRF3.csv",sep = ",",row.names = 1,header = TRUE,stringsAsFactors = FALSE)
# colnames(MMRF) <- gsub(".bam","",colnames(MMRF))
# MMRF_infp <- read.table("~/desktop/MMRF.csv",sep = ",",header = TRUE,stringsAsFactors = FALSE)
# colnames(MMRF) <- MMRF_info[MMRF_info$Run %in% colnames(MMRF),2]
# rownames(MMRF) <- paste0(rownames(MMRF),"_at")
# write.csv(MMRF,file = "MMRF3.csv")
#####################################################
library("geneClassifiers")
library(survival)
library(survC1)
#import MMRF expression data
# MMRF <- as.matrix(read.table("MMRF3.csv",sep = ",",row.names = 1,header = TRUE,stringsAsFactors = FALSE))
# #import survival information
# MMRF_infp <- read.table("~/desktop/MMRF.csv",sep = ",",header = TRUE,stringsAsFactors = FALSE)
# MMRF_infp <- na.omit(MMRF_infp)
# #sample ids with survial information
# MMRF_common <- intersect(MMRF_infp$public_id,colnames(MMRF)) #607
#save to RData and load
load("~/Documents/Research/MMRF_Prebs_express.RData")

#Create a calculation function for calculating risk scores
#Use original gene signatures without removing any probesets
calculation <- function(GeneSignature){
  #load required gene signature
  EMC92 <- data.frame(getWeights(getClassifier(GeneSignature)))
  probes <- rownames(EMC92)
  #Z-transformed expression data
  MMRF <- t(MMRF[probes,MMRF_common])
  MMRF <- t(scale(MMRF,center = TRUE, scale = TRUE))
  #set NAs to 0
  MMRF[is.na(MMRF)] <- 0
  #calculate risk scores
  EMC92_score <- rep(NA,dim(MMRF)[2])
  for (i in 1:dim(MMRF)[2]){
    express <- MMRF[probes,i]
    weight <- EMC92[,1]
    EMC92_score[i] <- sum(weight*express) 
  }
  EMC92_score <- as.data.frame(EMC92_score)
  rownames(EMC92_score) <- colnames(MMRF)
  #check common sample ids between expression data and survival data
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
  dicho <- sort(MMRF_infp2$score, decreasing = TRUE)[round(length(MMRF_common)*0.217)]
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
calculation("EMC92")
calculation("UAMS70")
calculation("UAMS80")
calculation("UAMS17")
calculation("HM19")
calculation("IFM15")
calculation("MILLENNIUM100")
#######################################
# Calculate risk score via MRCIX6 
#Z-transformed expression data
G_MRCIX6 <- read.table("MRCIX6_prebs.csv",sep = ",",header = TRUE,stringsAsFactors = FALSE)
MMRF <- t(MMRF[G_MRCIX6$affy_hg_u133_plus_2,MMRF_common])
MMRF <- t(scale(MMRF,center = TRUE, scale = TRUE))
#set NAs to 0
MMRF[is.na(MMRF)] <- 0
EMC92_score <- rep(NA,dim(MMRF)[2])
for (i in 1:dim(MMRF)[2]){
  if(MMRF["203755_at",i]/MMRF["216326_s_at",i] >1 | MMRF["200608_s_at",i]/MMRF["217731_s_at",i] >1){
    EMC92_score[i] <- 1
  }
  else{
    EMC92_score[i] <- 0
  }
}
################################
#Create a calculation function for calculating risk scores
#Use transformed gene signatures with reduced probesets
load("~/Documents/Research/MMRF_Prebs_express.RData")
calculation2 <- function(GeneSignature){
  #load required gene signature
  EMC92symbol <- read.table(GeneSignature,sep = ",",header = TRUE,stringsAsFactors = FALSE)
  probes <- EMC92symbol$probe
  #Z-transformed expression data
  MMRF <- t(MMRF[probes,MMRF_common])
  MMRF <- t(scale(MMRF,center = TRUE, scale = TRUE))
  #set NAs to 0
  MMRF[is.na(MMRF)] <- 0
  #calculate risk scores
  EMC92_score <- rep(NA,dim(MMRF)[2])
  for (i in 1:dim(MMRF)[2]){
    express <- MMRF[probes,i]
    weight <- EMC92symbol$weight
    EMC92_score[i] <- sum(weight*express) 
  }
  EMC92_score <- as.data.frame(EMC92_score)
  rownames(EMC92_score) <- colnames(MMRF)
  #check common sample ids between expression data and survival data
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
  dicho <- sort(MMRF_infp2$score, decreasing = TRUE)[round(length(MMRF_common)*0.217)]
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
calculation2("EMC92_mmrf.csv")
calculation2("UAMS70_prebs.csv")
calculation2("UAMS80_mmrf.csv")
calculation2("UAMS17_prebs.csv")
calculation2("HM19_prebs.csv")
calculation2("IFM15_prebs.csv")
calculation2("M100_mmrf.csv")
