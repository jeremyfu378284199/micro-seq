###########################################
#This R file checks the performance of each gene signatures on 3 microarray datasets
#calculate c-statistics,HR for different GSs from 3 microarray datasets
##########################################
#load required packages
source("http://bioconductor.org/biocLite.R")
biocLite("geneClassifiers")
biocLite("GEOquery")
biocLite("EnrichmentBrowser")
biocLite("sva")
install.packages('survival')
install.packages('survC1')
library("GEOquery")
library("sva")
library("geneClassifiers")
library(survival)
library(survC1)
########################################
#load matrix files HO65, APEX, UAMS 
HO65 <- exprs(getGEO(filename="~/Downloads/GSE19784_series_matrix.txt.gz"))
UAMSfull <- getGEO(filename = "~/Downloads/GSE24080_series_matrix.txt.gz")
#All UAMS subjects information
UAMSinfo <- UAMSfull@phenoData@data
#seprate UAMS data to TT2 and TT3
TT2 <- rownames(UAMSinfo[UAMSinfo$`maqc_distribution_status:ch1` %in% "Training",])
TT3 <- rownames(UAMSinfo[UAMSinfo$`maqc_distribution_status:ch1` %in% "Validation",])
UAMS <- exprs(getGEO(filename = "~/Downloads/GSE24080_series_matrix.txt.gz"))

APEXA <- log2(exprs(getGEO(filename = "~/Downloads/GSE9782-GPL96_series_matrix.txt.gz")))
APEXB <- log2(exprs(getGEO(filename = "~/Downloads/GSE9782-GPL97_series_matrix.txt.gz")))
APEX <- rbind(APEXA,APEXB)

###############################################
#get different GSs weights
###############################################
#################HO65
#load HO65 survial information
HO65_sample_information <- read.table("~/desktop/HO65 sample information2.csv",sep = ",",header = TRUE,stringsAsFactors = FALSE,row.names = 1)
#######EMC92
#calculate the risk scores of each subject from HO65 
EMC92 <- data.frame(getWeights(getClassifier("EMC92")))
EMC92_score_HO65 <- rep(NA,328)
for (i in 1:328){
  express <- HO65[rownames(EMC92),i]
  weight <- EMC92$getWeights.getClassifier..EMC92...
  EMC92_score_HO65[i] <- sum(weight*express) 
}
EMC92_score_HO65 <- as.data.frame(EMC92_score_HO65)
rownames(EMC92_score_HO65) <- colnames(HO65)
HO65_sample_information$EMC92 <- EMC92_score_HO65[rownames(HO65_sample_information),]
#########################
#base on top 21.7% risk score
dicho <- sort(EMC92_score_HO65$EMC92_score_HO65,decreasing = TRUE)[round(328*0.217)]
HO65_sample_information$EMC92percent <- as.numeric(HO65_sample_information$EMC92 >dicho)
###########################
HO65_sample_information$SurvObj <- with(HO65_sample_information, Surv(survival_time, survival_indicatior))
HO65_sample_information$SurvObj2 <- with(HO65_sample_information, Surv(progression_time, progression_indicator))
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(HO65_sample_information$survival_time,HO65_sample_information$survival_indicatior)~HO65_sample_information$EMC92percent))
summary(coxph(Surv(HO65_sample_information$progression_time,HO65_sample_information$progression_indicator)~HO65_sample_information$EMC92percent))
#Hazard ratio of OS prediction using continuous risk scores
summary(coxph(Surv(HO65_sample_information$survival_time,HO65_sample_information$survival_indicatior)~HO65_sample_information$EMC92)) 
summary(coxph(Surv(HO65_sample_information$progression_time,HO65_sample_information$progression_indicator)~HO65_sample_information$EMC92))

###############################
#####UAMS70
UAMS70 <- data.frame(getWeights(getClassifier("UAMS70")))
UAMS70_score_HO65 <- rep(NA,328)
for (i in 1:328){
  express <- HO65[rownames(UAMS70),i]
  weight <- UAMS70$getWeights.getClassifier..UAMS70...
  UAMS70_score_HO65[i] <- sum(weight*express) 
}
UAMS70_score_HO65 <- as.data.frame(UAMS70_score_HO65)
rownames(UAMS70_score_HO65) <- colnames(HO65)
HO65_sample_information$UAMS70 <- UAMS70_score_HO65[rownames(HO65_sample_information),]
#########################
#base on top 12.4% risk score
dicho <- sort(UAMS70_score_HO65$UAMS70_score_HO65,decreasing = TRUE)[round(328*0.124)]
HO65_sample_information$UAMS70percent <- as.numeric(HO65_sample_information$UAMS70 >dicho)
###########################
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(HO65_sample_information$survival_time,HO65_sample_information$survival_indicatior)~HO65_sample_information$UAMS70percent)) 
summary(coxph(Surv(HO65_sample_information$progression_time,HO65_sample_information$progression_indicator)~HO65_sample_information$UAMS70percent))
#Hazard ratio of OS prediction using continuous risk scores
summary(coxph(Surv(HO65_sample_information$survival_time,HO65_sample_information$survival_indicatior)~HO65_sample_information$UAMS70)) 
summary(coxph(Surv(HO65_sample_information$progression_time,HO65_sample_information$progression_indicator)~HO65_sample_information$UAMS70))
####UAMS80
UAMS80 <- data.frame(getWeights(getClassifier("UAMS80")))
UAMS80_score_HO65 <- rep(NA,328)
for (i in 1:328){
  express <- HO65[rownames(UAMS80),i]
  weight <- UAMS80$getWeights.getClassifier..UAMS80...
  UAMS80_score_HO65[i] <- sum(weight*express) 
}
UAMS80_score_HO65 <- as.data.frame(UAMS80_score_HO65)
rownames(UAMS80_score_HO65) <- colnames(HO65)
HO65_sample_information$UAMS80 <- UAMS80_score_HO65[rownames(HO65_sample_information),]
#########################
#base on top 25.9% risk score
dicho <- sort(UAMS80_score_HO65$UAMS80_score_HO65,decreasing = TRUE)[round(328*0.259)]
HO65_sample_information$UAMS80percent <- as.numeric(HO65_sample_information$UAMS80 >dicho)
###########################
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(HO65_sample_information$survival_time,HO65_sample_information$survival_indicatior)~HO65_sample_information$UAMS80percent)) 
summary(coxph(Surv(HO65_sample_information$progression_time,HO65_sample_information$progression_indicator)~HO65_sample_information$UAMS80percent))
#Hazard ratio of OS prediction using continuous risk scores
summary(coxph(Surv(HO65_sample_information$survival_time,HO65_sample_information$survival_indicatior)~HO65_sample_information$UAMS80)) 
summary(coxph(Surv(HO65_sample_information$progression_time,HO65_sample_information$progression_indicator)~HO65_sample_information$UAMS80))
####UAMS17
UAMS17 <- data.frame(getWeights(getClassifier("UAMS17")))
UAMS17_score_HO65 <- rep(NA,328)
for (i in 1:328){
  express <- HO65[rownames(UAMS17),i]
  weight <- UAMS17$getWeights.getClassifier..UAMS17...
  UAMS17_score_HO65[i] <- sum(weight*express) 
}
UAMS17_score_HO65 <- as.data.frame(UAMS17_score_HO65)
rownames(UAMS17_score_HO65) <- colnames(HO65)
HO65_sample_information$UAMS17 <- UAMS17_score_HO65[rownames(HO65_sample_information),]
#########################
#base on top 12.1% risk score
dicho <- sort(UAMS17_score_HO65$UAMS17_score_HO65,decreasing = TRUE)[round(328*0.121)]
HO65_sample_information$UAMS17percent <- as.numeric(HO65_sample_information$UAMS17 >dicho)
###########################
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(HO65_sample_information$survival_time,HO65_sample_information$survival_indicatior)~HO65_sample_information$UAMS17percent)) 
summary(coxph(Surv(HO65_sample_information$progression_time,HO65_sample_information$progression_indicator)~HO65_sample_information$UAMS17percent))
#Hazard ratio of OS prediction using continuous risk scores
summary(coxph(Surv(HO65_sample_information$survival_time,HO65_sample_information$survival_indicatior)~HO65_sample_information$UAMS17)) 
summary(coxph(Surv(HO65_sample_information$progression_time,HO65_sample_information$progression_indicator)~HO65_sample_information$UAMS17))
####HM19
HM19 <- data.frame(getWeights(getClassifier("HM19")))
HM19_score_HO65 <- rep(NA,328)
for (i in 1:328){
  express <- HO65[rownames(HM19),i]
  weight <- HM19$getWeights.getClassifier..HM19...
  HM19_score_HO65[i] <- sum(weight*express) 
}
HM19_score_HO65 <- as.data.frame(HM19_score_HO65)
rownames(HM19_score_HO65) <- colnames(HO65)
HO65_sample_information$HM19 <- HM19_score_HO65[rownames(HO65_sample_information),]
#########################
#base on top 21.7% risk score
dicho <- sort(HM19_score_HO65$HM19_score_HO65,decreasing = TRUE)[round(328*0.217)]
HO65_sample_information$HM19percent <- as.numeric(HO65_sample_information$HM19 >dicho)
###########################
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(HO65_sample_information$survival_time,HO65_sample_information$survival_indicatior)~HO65_sample_information$HM19percent)) 
summary(coxph(Surv(HO65_sample_information$progression_time,HO65_sample_information$progression_indicator)~HO65_sample_information$HM19percent))
#Hazard ratio of OS prediction using continuous risk scores
summary(coxph(Surv(HO65_sample_information$survival_time,HO65_sample_information$survival_indicatior)~HO65_sample_information$HM19)) 
summary(coxph(Surv(HO65_sample_information$progression_time,HO65_sample_information$progression_indicator)~HO65_sample_information$HM19))
####IFM15
IFM15 <- data.frame(getWeights(getClassifier("IFM15")))
IFM15_score_HO65 <- rep(NA,328)
for (i in 1:328){
  express <- HO65[rownames(IFM15),i]
  weight <- IFM15$getWeights.getClassifier..IFM15...
  IFM15_score_HO65[i] <- sum(weight*express) 
}
IFM15_score_HO65 <- as.data.frame(IFM15_score_HO65)
rownames(IFM15_score_HO65) <- colnames(HO65)
HO65_sample_information$IFM15 <- IFM15_score_HO65[rownames(HO65_sample_information),]
#########################
#base on top 28.6% risk score
dicho <- sort(IFM15_score_HO65$IFM15_score_HO65,decreasing = TRUE)[round(328*0.286)]
HO65_sample_information$IFM15percent <- as.numeric(HO65_sample_information$IFM15 >dicho)
###########################
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(HO65_sample_information$survival_time,HO65_sample_information$survival_indicatior)~HO65_sample_information$IFM15percent)) 
summary(coxph(Surv(HO65_sample_information$progression_time,HO65_sample_information$progression_indicator)~HO65_sample_information$IFM15percent))
#Hazard ratio of OS prediction using continuous risk scores
summary(coxph(Surv(HO65_sample_information$survival_time,HO65_sample_information$survival_indicatior)~HO65_sample_information$IFM15)) 
summary(coxph(Surv(HO65_sample_information$progression_time,HO65_sample_information$progression_indicator)~HO65_sample_information$IFM15))
####MRCIX6
MRCIX6_score_HO65 <- rep(NA,328)
for (i in 1:328){
  if(HO65["203755_at",i]/HO65["216326_s_at",i] >1 | HO65["203213_at",i]/HO65["218034_at",i] >1 | HO65["200608_s_at",i]/HO65["217731_s_at",i] >1){
    MRCIX6_score_HO65[i] <- 1
  }
  else{
    MRCIX6_score_HO65[i] <- 0
  }
}
MRCIX6_score_HO65 <- as.data.frame(MRCIX6_score_HO65)
rownames(MRCIX6_score_HO65) <- colnames(HO65)
HO65_sample_information$MRCIX6 <- MRCIX6_score_HO65[rownames(HO65_sample_information),]
#########################
#Hazard ratio 
summary(coxph(Surv(HO65_sample_information$survival_time,HO65_sample_information$survival_indicatior)~HO65_sample_information$MRCIX6)) 
summary(coxph(Surv(HO65_sample_information$progression_time,HO65_sample_information$progression_indicator)~HO65_sample_information$MRCIX6))

####MILLENNIUM100
MILLENNIUM100<- data.frame(getWeights(getClassifier("MILLENNIUM100")))
MILLENNIUM100_score_HO65 <- rep(NA,328)
for (i in 1:328){
  express <- HO65[rownames(MILLENNIUM100),i]
  weight <- MILLENNIUM100$getWeights.getClassifier..MILLENNIUM100...
  MILLENNIUM100_score_HO65[i] <- sum(weight*express) 
}
MILLENNIUM100_score_HO65 <- as.data.frame(MILLENNIUM100_score_HO65)
rownames(MILLENNIUM100_score_HO65) <- colnames(HO65)
HO65_sample_information$MILLENNIUM100 <- MILLENNIUM100_score_HO65[rownames(HO65_sample_information),]
#########################
#base on top 55.9% risk score
dicho <- sort(MILLENNIUM100_score_HO65$MILLENNIUM100_score_HO65,decreasing = TRUE)[round(328*0.559)]
HO65_sample_information$MILLENNIUM100percent <- as.numeric(HO65_sample_information$MILLENNIUM100 >dicho)
###########################
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(HO65_sample_information$survival_time,HO65_sample_information$survival_indicatior)~HO65_sample_information$MILLENNIUM100percent)) 
summary(coxph(Surv(HO65_sample_information$progression_time,HO65_sample_information$progression_indicator)~HO65_sample_information$MILLENNIUM100percent))
#Hazard ratio of OS prediction using continuous risk scores
summary(coxph(Surv(HO65_sample_information$survival_time,HO65_sample_information$survival_indicatior)~HO65_sample_information$MILLENNIUM100)) 
summary(coxph(Surv(HO65_sample_information$progression_time,HO65_sample_information$progression_indicator)~HO65_sample_information$MILLENNIUM100))
#################################################
inputdata <- apply(HO65_sample_information,2,as.numeric)
set.seed(123)
#C-statistics of OS prediciton
#EMC92
Est.Cval(inputdata[,c(2,1,5)],2*365.25,nofit = TRUE)[["Dhat"]]
#UAMS70
Est.Cval(inputdata[,c(2,1,6)],2*365.25,nofit = TRUE)[["Dhat"]]
#UAMS17
Est.Cval(inputdata[,c(2,1,7)],2*365.25,nofit = TRUE)[["Dhat"]]
#UAMS80
Est.Cval(inputdata[,c(2,1,8)],2*365.25,nofit = TRUE)[["Dhat"]]
#HM19
Est.Cval(inputdata[,c(2,1,9)],2*365.25,nofit = TRUE)[["Dhat"]]
#IFM15
Est.Cval(inputdata[,c(2,1,10)],2*365.25,nofit = TRUE)[["Dhat"]]
#MRCIX6
Est.Cval(inputdata[,c(2,1,11)],2*365.25,nofit = TRUE)[["Dhat"]]
#MILLENNIUM100
Est.Cval(inputdata[,c(2,1,12)],2*365.25,nofit = TRUE)[["Dhat"]]
#######################################################
#C-statistics of PFS prediciton
#EMC92
Est.Cval(inputdata[,c(4,3,5)],5*365.25,nofit = TRUE)[["Dhat"]]
#UAMS70
Est.Cval(inputdata[,c(4,3,6)],5*365.25,nofit = TRUE)[["Dhat"]]
#UAMS17
Est.Cval(inputdata[,c(4,3,7)],5*365.25,nofit = TRUE)[["Dhat"]]
#UAMS80
Est.Cval(inputdata[,c(4,3,8)],5*365.25,nofit = TRUE)[["Dhat"]]
#HM19
Est.Cval(inputdata[,c(4,3,9)],5*365.25,nofit = TRUE)[["Dhat"]]
#IFM15
Est.Cval(inputdata[,c(4,3,10)],5*365.25,nofit = TRUE)[["Dhat"]]
#MRCIX6
Est.Cval(inputdata[,c(4,3,11)],5*365.25,nofit = TRUE)[["Dhat"]]
#MILLENNIUM100
Est.Cval(inputdata[,c(4,3,12)],5*365.25,nofit = TRUE)[["Dhat"]]

###########################################################
##########################################################
#########################################################
#load UAMS survial information
UAMS_sample_information <- read.table("~/desktop/UAMS sample information2.csv",sep = ",",row.names = 1,header = TRUE)
#######EMC92
EMC92 <- data.frame(getWeights(getClassifier("EMC92")))
EMC92_score_UAMS <- rep(NA,559)
for (i in 1:559){
  express <- UAMS[rownames(EMC92),i]
  weight <- EMC92$getWeights.getClassifier..EMC92...
  EMC92_score_UAMS[i] <- sum(weight*express) 
}
EMC92_score_UAMS <- as.data.frame(EMC92_score_UAMS)
rownames(EMC92_score_UAMS) <- colnames(UAMS)
UAMS_sample_information$EMC92 <- EMC92_score_UAMS[rownames(UAMS_sample_information),]
UAMS_TT2 <- UAMS_sample_information[TT2,]
UAMS_TT3 <- UAMS_sample_information[TT3,]
#########################
#base on top 19.4% risk score
dicho <- sort(UAMS_TT2$EMC92,decreasing = TRUE)[round(340*0.194)]
UAMS_TT2$EMC92percent <- as.numeric(UAMS_TT2$EMC92 >dicho)
#base on top 16.2% risk score
dicho <- sort(UAMS_TT3$EMC92,decreasing = TRUE)[round(214*0.162)]
UAMS_TT3$EMC92percent <- as.numeric(UAMS_TT3$EMC92 >dicho)
###########################
UAMS_sample_information$SurvObj <- with(UAMS_sample_information, Surv(survival_time, survival_indicatior))
UAMS_sample_information$SurvObj2 <- with(UAMS_sample_information, Surv(progression_time, progression_indicator))
#TT2 HR 
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(UAMS_TT2$survival_time,UAMS_TT2$survival_indicatior)~UAMS_TT2$EMC92percent)) 
summary(coxph(Surv(UAMS_TT2$progression_time,UAMS_TT2$progression_indicator)~UAMS_TT2$EMC92percent))
#Hazard ratio of OS prediction using continuous risk scores
summary(coxph(Surv(UAMS_TT2$survival_time,UAMS_TT2$survival_indicatior)~UAMS_TT2$EMC92)) 
summary(coxph(Surv(UAMS_TT2$progression_time,UAMS_TT2$progression_indicator)~UAMS_TT2$EMC92))
#TT3 HR
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(UAMS_TT3$survival_time,UAMS_TT3$survival_indicatior)~UAMS_TT3$EMC92percent)) 
summary(coxph(Surv(UAMS_TT3$progression_time,UAMS_TT3$progression_indicator)~UAMS_TT3$EMC92percent))
#Hazard ratio of OS prediction using continuous risk scores
summary(coxph(Surv(UAMS_TT3$survival_time,UAMS_TT3$survival_indicatior)~UAMS_TT3$EMC92)) 
summary(coxph(Surv(UAMS_TT3$progression_time,UAMS_TT3$progression_indicator)~UAMS_TT3$EMC92))
#####UAMS70
UAMS70 <- data.frame(getWeights(getClassifier("UAMS70")))
UAMS70_score_UAMS <- rep(NA,559)
for (i in 1:559){
  express <- UAMS[rownames(UAMS70),i]
  weight <- UAMS70$getWeights.getClassifier..UAMS70...
  UAMS70_score_UAMS[i] <- sum(weight*express) 
}
UAMS70_score_UAMS <- as.data.frame(UAMS70_score_UAMS)
rownames(UAMS70_score_UAMS) <- colnames(UAMS)
UAMS_sample_information$UAMS70 <- UAMS70_score_UAMS[rownames(UAMS_sample_information),]
UAMS_TT2 <- UAMS_sample_information[TT2,]
UAMS_TT3 <- UAMS_sample_information[TT3,]
#########################
#base on top 17.1% risk score
dicho <- sort(UAMS_TT2$UAMS70,decreasing = TRUE)[round(340*0.171)]
UAMS_TT2$UAMS70percent <- as.numeric(UAMS_TT2$UAMS70 >dicho)
#base on top 14.8% risk score
dicho <- sort(UAMS_TT3$UAMS70,decreasing = TRUE)[round(214*0.148)]
UAMS_TT3$UAMS70percent <- as.numeric(UAMS_TT3$UAMS70 >dicho)
###########################
#TT2 HR
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(UAMS_TT2$survival_time,UAMS_TT2$survival_indicatior)~UAMS_TT2$UAMS70percent)) 
summary(coxph(Surv(UAMS_TT2$progression_time,UAMS_TT2$progression_indicator)~UAMS_TT2$UAMS70percent))
#Hazard ratio of OS prediction using continuous risk scores
summary(coxph(Surv(UAMS_TT2$survival_time,UAMS_TT2$survival_indicatior)~UAMS_TT2$UAMS70)) 
summary(coxph(Surv(UAMS_TT2$progression_time,UAMS_TT2$progression_indicator)~UAMS_TT2$UAMS70))
#TT3 HR
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(UAMS_TT3$survival_time,UAMS_TT3$survival_indicatior)~UAMS_TT3$UAMS70percent)) 
summary(coxph(Surv(UAMS_TT3$progression_time,UAMS_TT3$progression_indicator)~UAMS_TT3$UAMS70percent))
#Hazard ratio of OS prediction using continuous risk scores
summary(coxph(Surv(UAMS_TT3$survival_time,UAMS_TT3$survival_indicatior)~UAMS_TT3$UAMS70)) 
summary(coxph(Surv(UAMS_TT3$progression_time,UAMS_TT3$progression_indicator)~UAMS_TT3$UAMS70))
####UAMS80
UAMS80 <- data.frame(getWeights(getClassifier("UAMS80")))
UAMS80_score_UAMS <- rep(NA,559)
for (i in 1:559){
  express <- UAMS[rownames(UAMS80),i]
  weight <- UAMS80$getWeights.getClassifier..UAMS80...
  UAMS80_score_UAMS[i] <- sum(weight*express) 
}
UAMS80_score_UAMS <- as.data.frame(UAMS80_score_UAMS)
rownames(UAMS80_score_UAMS) <- colnames(UAMS)
UAMS_sample_information$UAMS80 <- UAMS80_score_UAMS[rownames(UAMS_sample_information),]
UAMS_TT2 <- UAMS_sample_information[TT2,]
UAMS_TT3 <- UAMS_sample_information[TT3,]
#########################
#base on top 25.6% risk score
dicho <- sort(UAMS_TT2$UAMS80,decreasing = TRUE)[round(340*0.256)]
UAMS_TT2$UAMS80percent <- as.numeric(UAMS_TT2$UAMS80 >dicho)
#base on top 23.9% risk score
dicho <- sort(UAMS_TT3$UAMS80,decreasing = TRUE)[round(214*0.239)]
UAMS_TT3$UAMS80percent <- as.numeric(UAMS_TT3$UAMS80 >dicho)
###########################
#TT2 HR
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(UAMS_TT2$survival_time,UAMS_TT2$survival_indicatior)~UAMS_TT2$UAMS80percent)) 
summary(coxph(Surv(UAMS_TT2$progression_time,UAMS_TT2$progression_indicator)~UAMS_TT2$UAMS80percent))
#Hazard ratio of OS prediction using continuous risk scores
summary(coxph(Surv(UAMS_TT2$survival_time,UAMS_TT2$survival_indicatior)~UAMS_TT2$UAMS80)) 
summary(coxph(Surv(UAMS_TT2$progression_time,UAMS_TT2$progression_indicator)~UAMS_TT2$UAMS80))
#TT3 HR
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(UAMS_TT3$survival_time,UAMS_TT3$survival_indicatior)~UAMS_TT3$UAMS80percent)) 
summary(coxph(Surv(UAMS_TT3$progression_time,UAMS_TT3$progression_indicator)~UAMS_TT3$UAMS80percent))
#Hazard ratio of OS prediction using continuous risk scores
summary(coxph(Surv(UAMS_TT3$survival_time,UAMS_TT3$survival_indicatior)~UAMS_TT3$UAMS80)) 
summary(coxph(Surv(UAMS_TT3$progression_time,UAMS_TT3$progression_indicator)~UAMS_TT3$UAMS80))
####UAMS17
UAMS17 <- data.frame(getWeights(getClassifier("UAMS17")))
UAMS17_score_UAMS <- rep(NA,559)
for (i in 1:559){
  express <- UAMS[rownames(UAMS17),i]
  weight <- UAMS17$getWeights.getClassifier..UAMS17...
  UAMS17_score_UAMS[i] <- sum(weight*express) 
}
UAMS17_score_UAMS <- as.data.frame(UAMS17_score_UAMS)
rownames(UAMS17_score_UAMS) <- colnames(UAMS)
UAMS_sample_information$UAMS17 <- UAMS17_score_UAMS[rownames(UAMS_sample_information),]
UAMS_TT2 <- UAMS_sample_information[TT2,]
UAMS_TT3 <- UAMS_sample_information[TT3,]
#########################
#base on top 14% risk score
dicho <- sort(UAMS_TT2$UAMS17,decreasing = TRUE)[round(340*0.14)]
UAMS_TT2$UAMS17percent <- as.numeric(UAMS_TT2$UAMS17 >dicho)
#base on top 13.4% risk score
dicho <- sort(UAMS_TT3$UAMS17,decreasing = TRUE)[round(214*0.134)]
UAMS_TT3$UAMS17percent <- as.numeric(UAMS_TT3$UAMS17 >dicho)
###########################
#TT2 HR
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(UAMS_TT2$survival_time,UAMS_TT2$survival_indicatior)~UAMS_TT2$UAMS17percent)) 
summary(coxph(Surv(UAMS_TT2$progression_time,UAMS_TT2$progression_indicator)~UAMS_TT2$UAMS17percent))
#Hazard ratio of OS prediction using continuous risk scores
summary(coxph(Surv(UAMS_TT2$survival_time,UAMS_TT2$survival_indicatior)~UAMS_TT2$UAMS17)) 
summary(coxph(Surv(UAMS_TT2$progression_time,UAMS_TT2$progression_indicator)~UAMS_TT2$UAMS17))
#TT3 HR
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(UAMS_TT3$survival_time,UAMS_TT3$survival_indicatior)~UAMS_TT3$UAMS17percent)) 
summary(coxph(Surv(UAMS_TT3$progression_time,UAMS_TT3$progression_indicator)~UAMS_TT3$UAMS17percent))
#Hazard ratio of OS prediction using continuous risk scores
summary(coxph(Surv(UAMS_TT3$survival_time,UAMS_TT3$survival_indicatior)~UAMS_TT3$UAMS17)) 
summary(coxph(Surv(UAMS_TT3$progression_time,UAMS_TT3$progression_indicator)~UAMS_TT3$UAMS17))
####HM19
HM19 <- data.frame(getWeights(getClassifier("HM19")))
HM19_score_UAMS <- rep(NA,559)
for (i in 1:559){
  express <- UAMS[rownames(HM19),i]
  weight <- HM19$getWeights.getClassifier..HM19...
  HM19_score_UAMS[i] <- sum(weight*express) 
}
HM19_score_UAMS <- as.data.frame(HM19_score_UAMS)
rownames(HM19_score_UAMS) <- colnames(UAMS)
UAMS_sample_information$HM19 <- HM19_score_UAMS[rownames(UAMS_sample_information),]
UAMS_TT2 <- UAMS_sample_information[TT2,]
UAMS_TT3 <- UAMS_sample_information[TT3,]
#########################
#base on top 21.7% risk score
dicho <- sort(UAMS_TT2$HM19,decreasing = TRUE)[round(340*0.217)]
UAMS_TT2$HM19percent <- as.numeric(UAMS_TT2$HM19 >dicho)
#base on top 21.7% risk score
dicho <- sort(UAMS_TT3$HM19,decreasing = TRUE)[round(214*0.217)]
UAMS_TT3$HM19percent <- as.numeric(UAMS_TT3$HM19 >dicho)
###########################
#TT2 HR
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(UAMS_TT2$survival_time,UAMS_TT2$survival_indicatior)~UAMS_TT2$HM19percent)) 
summary(coxph(Surv(UAMS_TT2$progression_time,UAMS_TT2$progression_indicator)~UAMS_TT2$HM19percent))
#Hazard ratio of OS prediction using continuous risk scores
summary(coxph(Surv(UAMS_TT2$survival_time,UAMS_TT2$survival_indicatior)~UAMS_TT2$HM19)) 
summary(coxph(Surv(UAMS_TT2$progression_time,UAMS_TT2$progression_indicator)~UAMS_TT2$HM19))
#TT3 HR
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(UAMS_TT3$survival_time,UAMS_TT3$survival_indicatior)~UAMS_TT3$HM19percent)) 
summary(coxph(Surv(UAMS_TT3$progression_time,UAMS_TT3$progression_indicator)~UAMS_TT3$HM19percent))
#Hazard ratio of OS prediction using continuous risk scores
summary(coxph(Surv(UAMS_TT3$survival_time,UAMS_TT3$survival_indicatior)~UAMS_TT3$HM19)) 
summary(coxph(Surv(UAMS_TT3$progression_time,UAMS_TT3$progression_indicator)~UAMS_TT3$HM19))
####IFM15
IFM15 <- data.frame(getWeights(getClassifier("IFM15")))
IFM15_score_UAMS <- rep(NA,559)
for (i in 1:559){
  express <- UAMS[rownames(IFM15),i]
  weight <- IFM15$getWeights.getClassifier..IFM15...
  IFM15_score_UAMS[i] <- sum(weight*express) 
}
IFM15_score_UAMS <- as.data.frame(IFM15_score_UAMS)
rownames(IFM15_score_UAMS) <- colnames(UAMS)
UAMS_sample_information$IFM15 <- IFM15_score_UAMS[rownames(UAMS_sample_information),]
UAMS_TT2 <- UAMS_sample_information[TT2,]
UAMS_TT3 <- UAMS_sample_information[TT3,]
#########################
#base on top 24.8% risk score
dicho <- sort(UAMS_TT2$IFM15,decreasing = TRUE)[round(340*0.248)]
UAMS_TT2$IFM15percent <- as.numeric(UAMS_TT2$IFM15 >dicho)
#base on top 26.1% risk score
dicho <- sort(UAMS_TT3$IFM15,decreasing = TRUE)[round(214*0.261)]
UAMS_TT3$IFM15percent <- as.numeric(UAMS_TT3$IFM15 >dicho)

###########################
#TT2 HR
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(UAMS_TT2$survival_time,UAMS_TT2$survival_indicatior)~UAMS_TT2$IFM15percent)) 
summary(coxph(Surv(UAMS_TT2$progression_time,UAMS_TT2$progression_indicator)~UAMS_TT2$IFM15percent))
#Hazard ratio of OS prediction using continuous risk scores
summary(coxph(Surv(UAMS_TT2$survival_time,UAMS_TT2$survival_indicatior)~UAMS_TT2$IFM15)) 
summary(coxph(Surv(UAMS_TT2$progression_time,UAMS_TT2$progression_indicator)~UAMS_TT2$IFM15))
#TT3 HR
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(UAMS_TT3$survival_time,UAMS_TT3$survival_indicatior)~UAMS_TT3$IFM15percent)) 
summary(coxph(Surv(UAMS_TT3$progression_time,UAMS_TT3$progression_indicator)~UAMS_TT3$IFM15percent))
#Hazard ratio of OS prediction using continuous risk scores
summary(coxph(Surv(UAMS_TT3$survival_time,UAMS_TT3$survival_indicatior)~UAMS_TT3$IFM15)) 
summary(coxph(Surv(UAMS_TT3$progression_time,UAMS_TT3$progression_indicator)~UAMS_TT3$IFM15))
####MRCIX6
MRCIX6_score_UAMS <- rep(NA,559)
for (i in 1:559){
  if(UAMS["203755_at",i]/UAMS["216326_s_at",i] >1 | UAMS["203213_at",i]/UAMS["218034_at",i] >1 | UAMS["200608_s_at",i]/UAMS["217731_s_at",i] >1){
    MRCIX6_score_UAMS[i] <- 1
  }
  else{
    MRCIX6_score_UAMS[i] <- 0
  }
}
MRCIX6_score_UAMS <- as.data.frame(MRCIX6_score_UAMS)
rownames(MRCIX6_score_UAMS) <- colnames(UAMS)
UAMS_sample_information$MRCIX6 <- MRCIX6_score_UAMS[rownames(UAMS_sample_information),]
UAMS_TT2 <- UAMS_sample_information[TT2,]
UAMS_TT3 <- UAMS_sample_information[TT3,]
#########################

#TT2 HR
summary(coxph(Surv(UAMS_TT2$survival_time,UAMS_TT2$survival_indicatior)~UAMS_TT2$MRCIX6)) 
summary(coxph(Surv(UAMS_TT2$progression_time,UAMS_TT2$progression_indicator)~UAMS_TT2$MRCIX6))
#TT3 HR
summary(coxph(Surv(UAMS_TT3$survival_time,UAMS_TT3$survival_indicatior)~UAMS_TT3$MRCIX6)) 
summary(coxph(Surv(UAMS_TT3$progression_time,UAMS_TT3$progression_indicator)~UAMS_TT3$MRCIX6))

####MILLENNIUM100
MILLENNIUM100<- data.frame(getWeights(getClassifier("MILLENNIUM100")))
MILLENNIUM100_score_UAMS <- rep(NA,559)
for (i in 1:559){
  express <- UAMS[rownames(MILLENNIUM100),i]
  weight <- MILLENNIUM100$getWeights.getClassifier..MILLENNIUM100...
  MILLENNIUM100_score_UAMS[i] <- sum(weight*express) 
}
MILLENNIUM100_score_UAMS <- as.data.frame(MILLENNIUM100_score_UAMS)
rownames(MILLENNIUM100_score_UAMS) <- colnames(UAMS)
UAMS_sample_information$MILLENNIUM100 <- MILLENNIUM100_score_UAMS[rownames(UAMS_sample_information),]
UAMS_TT2 <- UAMS_sample_information[TT2,]
UAMS_TT3 <- UAMS_sample_information[TT3,]
#########################
#base on top 46.4% risk score
dicho <- sort(UAMS_TT2$MILLENNIUM100,decreasing = TRUE)[round(340*0.464)]
UAMS_TT2$MILLENNIUM100percent <- as.numeric(UAMS_TT2$MILLENNIUM100 >dicho)
#base on top 53.5% risk score
dicho <- sort(UAMS_TT3$MILLENNIUM100,decreasing = TRUE)[round(214*0.535)]
UAMS_TT3$MILLENNIUM100percent <- as.numeric(UAMS_TT3$MILLENNIUM100 >dicho)

###########################
#TT2 HR
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(UAMS_TT2$survival_time,UAMS_TT2$survival_indicatior)~UAMS_TT2$MILLENNIUM100percent)) 
summary(coxph(Surv(UAMS_TT2$progression_time,UAMS_TT2$progression_indicator)~UAMS_TT2$MILLENNIUM100percent))
#Hazard ratio of OS prediction using continuous risk scores
summary(coxph(Surv(UAMS_TT2$survival_time,UAMS_TT2$survival_indicatior)~UAMS_TT2$MILLENNIUM100)) 
summary(coxph(Surv(UAMS_TT2$progression_time,UAMS_TT2$progression_indicator)~UAMS_TT2$MILLENNIUM100))
#TT3 HR
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(UAMS_TT3$survival_time,UAMS_TT3$survival_indicatior)~UAMS_TT3$MILLENNIUM100percent)) 
summary(coxph(Surv(UAMS_TT3$progression_time,UAMS_TT3$progression_indicator)~UAMS_TT3$MILLENNIUM100percent))
#Hazard ratio of OS prediction using continuous risk scores
summary(coxph(Surv(UAMS_TT3$survival_time,UAMS_TT3$survival_indicatior)~UAMS_TT3$MILLENNIUM100)) 
summary(coxph(Surv(UAMS_TT3$progression_time,UAMS_TT3$progression_indicator)~UAMS_TT3$MILLENNIUM100))

######################################################
#inputdata <- apply(UAMS_TT2,2,as.numeric)
inputdata <- apply(UAMS_TT3,2,as.numeric)
#survival
#EMC92
Est.Cval(inputdata[,c(2,1,5)],2*12,nofit = TRUE)[["Dhat"]]
#UAMS70
Est.Cval(inputdata[,c(2,1,6)],2*12,nofit = TRUE)[["Dhat"]]
#UAMS17
Est.Cval(inputdata[,c(2,1,7)],2*12,nofit = TRUE)[["Dhat"]]
#UAMS80
Est.Cval(inputdata[,c(2,1,8)],2*12,nofit = TRUE)[["Dhat"]]
#HM19
Est.Cval(inputdata[,c(2,1,9)],2*12,nofit = TRUE)[["Dhat"]]
#IFM15
Est.Cval(inputdata[,c(2,1,10)],2*12,nofit = TRUE)[["Dhat"]]
#MRCIX6
Est.Cval(inputdata[,c(2,1,11)],2*12,nofit = TRUE)[["Dhat"]]
#MILLENNIUM100
Est.Cval(inputdata[,c(2,1,12)],2*12,nofit = TRUE)[["Dhat"]]
#####################################
#progression 
#EMC92
Est.Cval(inputdata[,c(4,3,5)],5*12,nofit = TRUE)[["Dhat"]]
#UAMS70
Est.Cval(inputdata[,c(4,3,6)],5*12,nofit = TRUE)[["Dhat"]]
#UAMS17
Est.Cval(inputdata[,c(4,3,7)],5*12,nofit = TRUE)[["Dhat"]] 
#UAMS80
Est.Cval(inputdata[,c(4,3,8)],5*12,nofit = TRUE)[["Dhat"]]
#HM19
Est.Cval(inputdata[,c(4,3,9)],5*12,nofit = TRUE)[["Dhat"]]
#IFM15
Est.Cval(inputdata[,c(4,3,10)],5*12,nofit = TRUE)[["Dhat"]] 
#MRCIX6
Est.Cval(inputdata[,c(4,3,11)],5*12,nofit = TRUE)[["Dhat"]]
#MILLENNIUM100
Est.Cval(inputdata[,c(4,3,12)],5*12,nofit = TRUE)[["Dhat"]] 
###########################################################
##########################################################
#########################################################
#load APEX survial information
APEX_sample_information <- read.table("~/desktop/Apex sample information.csv",sep = ",",row.names = 1,header = TRUE)
#######EMC92
EMC92 <- data.frame(getWeights(getClassifier("EMC92")))
EMC92_score_APEX <- rep(NA,264)
for (i in 1:264){
  express <- APEX[rownames(EMC92),i]
  weight <- EMC92$getWeights.getClassifier..EMC92...
  EMC92_score_APEX[i] <- sum(weight*express) 
}
EMC92_score_APEX <- as.data.frame(EMC92_score_APEX)
rownames(EMC92_score_APEX) <- colnames(APEX)
APEX_sample_information$EMC92 <- EMC92_score_APEX[rownames(APEX_sample_information),]
#########################
#base on top 16.3% risk score
dicho <- sort(EMC92_score_APEX$EMC92_score_APEX,decreasing = TRUE)[round(328*0.163)]
APEX_sample_information$EMC92percent <- as.numeric(APEX_sample_information$EMC92 >dicho)
###########################
APEX_sample_information$SurvObj <- with(APEX_sample_information, Surv(survival_time, survival_indicatior))
APEX_sample_information$SurvObj2 <- with(APEX_sample_information, Surv(progression_time, progression_indicator))
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(APEX_sample_information$survival_time,APEX_sample_information$survival_indicatior)~APEX_sample_information$EMC92percent)) 
summary(coxph(Surv(APEX_sample_information$progression_time,APEX_sample_information$progression_indicator)~APEX_sample_information$EMC92percent))
#Hazard ratio of OS prediction using continuous risk scores
summary(coxph(Surv(APEX_sample_information$survival_time,APEX_sample_information$survival_indicatior)~APEX_sample_information$EMC92)) 
summary(coxph(Surv(APEX_sample_information$progression_time,APEX_sample_information$progression_indicator)~APEX_sample_information$EMC92))
#####UAMS70
UAMS70 <- data.frame(getWeights(getClassifier("UAMS70")))
diff <- setdiff(rownames(UAMS70),rownames(APEX))
index <- rep(NA,5)
for (i in 1:5){
  index[i] <- which(rownames(UAMS70) == diff[i])
}
UAMS65 <- data.frame(UAMS70[!(rownames(UAMS70) %in% diff),])
rownames(UAMS65) <- rownames(UAMS70)[-index]
UAMS65$UAMS70...rownames.UAMS70...in..diff.... <- c(rep(1/48,48),rep(-1/17,17))
UAMS70_score_APEX <- rep(NA,264)
for (i in 1:264){
  express <- APEX[rownames(UAMS65),i]
  weight <- UAMS65$UAMS70...rownames.UAMS70...in..diff....
  UAMS70_score_APEX[i] <- sum(weight*express) 
}
UAMS70_score_APEX <- as.data.frame(UAMS70_score_APEX)
rownames(UAMS70_score_APEX) <- colnames(APEX)
APEX_sample_information$UAMS70 <- UAMS70_score_APEX[rownames(APEX_sample_information),]
#########################
#base on top 6.4% risk score
dicho <- sort(UAMS70_score_APEX$UAMS70_score_APEX,decreasing = TRUE)[round(264*0.064)]
APEX_sample_information$UAMS70percent <- as.numeric(APEX_sample_information$UAMS70 >dicho)
###########################
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(APEX_sample_information$survival_time,APEX_sample_information$survival_indicatior)~APEX_sample_information$UAMS70percent)) 
summary(coxph(Surv(APEX_sample_information$progression_time,APEX_sample_information$progression_indicator)~APEX_sample_information$UAMS70percent))
#Hazard ratio of OS prediction using continuous risk scores
summary(coxph(Surv(APEX_sample_information$survival_time,APEX_sample_information$survival_indicatior)~APEX_sample_information$UAMS70)) 
summary(coxph(Surv(APEX_sample_information$progression_time,APEX_sample_information$progression_indicator)~APEX_sample_information$UAMS70))
####UAMS80
UAMS80 <- data.frame(getWeights(getClassifier("UAMS80")))
diff <- setdiff(rownames(UAMS80),rownames(APEX))
index <- rep(NA,5)
for (i in 1:5){
  index[i] <- which(rownames(UAMS80) == diff[i])
}
UAMS75 <- data.frame(UAMS80[!(rownames(UAMS80) %in% diff),])
rownames(UAMS75) <- rownames(UAMS80)[-index]
UAMS75$UAMS80...rownames.UAMS80...in..diff.... <- c(rep(1/42,42),rep(-1/33,33))
UAMS80_score_APEX <- rep(NA,264)
for (i in 1:264){
  express <- APEX[rownames(UAMS75),i]
  weight <- UAMS75$UAMS80...rownames.UAMS80...in..diff....
  UAMS80_score_APEX[i] <- sum(weight*express) 
}
UAMS80_score_APEX <- as.data.frame(UAMS80_score_APEX)
rownames(UAMS80_score_APEX) <- colnames(APEX)
APEX_sample_information$UAMS80 <- UAMS80_score_APEX[rownames(APEX_sample_information),]
#########################
#base on top 12.1% risk score
dicho <- sort(UAMS80_score_APEX$UAMS80_score_APEX,decreasing = TRUE)[round(264*0.121)]
APEX_sample_information$UAMS80percent <- as.numeric(APEX_sample_information$UAMS80 >dicho)
###########################
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(APEX_sample_information$survival_time,APEX_sample_information$survival_indicatior)~APEX_sample_information$UAMS80percent)) 
summary(coxph(Surv(APEX_sample_information$progression_time,APEX_sample_information$progression_indicator)~APEX_sample_information$UAMS80percent))
#Hazard ratio of OS prediction using continuous risk scores
summary(coxph(Surv(APEX_sample_information$survival_time,APEX_sample_information$survival_indicatior)~APEX_sample_information$UAMS80)) 
summary(coxph(Surv(APEX_sample_information$progression_time,APEX_sample_information$progression_indicator)~APEX_sample_information$UAMS80))

####UAMS17
UAMS17 <- data.frame(getWeights(getClassifier("UAMS17")))
diff <- setdiff(rownames(UAMS17),rownames(APEX))
UAMS16 <- data.frame(UAMS17[!(rownames(UAMS17) %in% diff),])
rownames(UAMS16) <- rownames(UAMS17)[-2]
UAMS17_score_APEX <- rep(NA,264)
for (i in 1:264){
  express <- APEX[rownames(UAMS16),i]
  weight <- UAMS16$UAMS17...rownames.UAMS17...in..diff....
  UAMS17_score_APEX[i] <- sum(weight*express) 
}
UAMS17_score_APEX <- as.data.frame(UAMS17_score_APEX)
rownames(UAMS17_score_APEX) <- colnames(APEX)
APEX_sample_information$UAMS17 <- UAMS17_score_APEX[rownames(APEX_sample_information),]
#########################
#base on top 10.6% risk score
dicho <- sort(UAMS17_score_APEX$UAMS17_score_APEX,decreasing = TRUE)[round(264*0.106)]
APEX_sample_information$UAMS17percent <- as.numeric(APEX_sample_information$UAMS17 >dicho)
###########################
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(APEX_sample_information$survival_time,APEX_sample_information$survival_indicatior)~APEX_sample_information$UAMS17percent)) 
summary(coxph(Surv(APEX_sample_information$progression_time,APEX_sample_information$progression_indicator)~APEX_sample_information$UAMS17percent))
#Hazard ratio of OS prediction using continuous risk scores
summary(coxph(Surv(APEX_sample_information$survival_time,APEX_sample_information$survival_indicatior)~APEX_sample_information$UAMS17)) 
summary(coxph(Surv(APEX_sample_information$progression_time,APEX_sample_information$progression_indicator)~APEX_sample_information$UAMS17))
####HM19
HM19 <- data.frame(getWeights(getClassifier("HM19")))
HM19_score_APEX <- rep(NA,264)
for (i in 1:264){
  express <- APEX[rownames(HM19),i]
  weight <- HM19$getWeights.getClassifier..HM19...
  HM19_score_APEX[i] <- sum(weight*express) 
}
HM19_score_APEX <- as.data.frame(HM19_score_APEX)
rownames(HM19_score_APEX) <- colnames(APEX)
APEX_sample_information$HM19 <- HM19_score_APEX[rownames(APEX_sample_information),]
#########################
#base on top 21.7% risk score
dicho <- sort(HM19_score_APEX$HM19_score_APEX,decreasing = TRUE)[round(264*0.217)]
APEX_sample_information$HM19percent <- as.numeric(APEX_sample_information$HM19 >dicho)
###########################
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(APEX_sample_information$survival_time,APEX_sample_information$survival_indicatior)~APEX_sample_information$HM19percent)) 
summary(coxph(Surv(APEX_sample_information$progression_time,APEX_sample_information$progression_indicator)~APEX_sample_information$HM19percent))
#Hazard ratio of OS prediction using continuous risk scores
summary(coxph(Surv(APEX_sample_information$survival_time,APEX_sample_information$survival_indicatior)~APEX_sample_information$HM19)) 
summary(coxph(Surv(APEX_sample_information$progression_time,APEX_sample_information$progression_indicator)~APEX_sample_information$HM19))
####IFM15
IFM15 <- data.frame(getWeights(getClassifier("IFM15")))
diff <- setdiff(rownames(IFM15),rownames(APEX))
which(rownames(IFM15) == diff)#15
IFM14 <- data.frame(IFM15[!(rownames(IFM15) %in% diff),])
rownames(IFM14) <- rownames(IFM15)[-15]
IFM15_score_APEX <- rep(NA,264)
for (i in 1:264){
  express <- APEX[rownames(IFM14),i]
  weight <- IFM14$IFM15...rownames.IFM15...in..diff....
  IFM15_score_APEX[i] <- sum(weight*express) 
}
IFM15_score_APEX <- as.data.frame(IFM15_score_APEX)
rownames(IFM15_score_APEX) <- colnames(APEX)
APEX_sample_information$IFM15 <- IFM15_score_APEX[rownames(APEX_sample_information),]
#########################
#base on top 24.6% risk score
dicho <- sort(IFM15_score_APEX$IFM15_score_APEX,decreasing = TRUE)[round(264*0.246)]
APEX_sample_information$IFM15percent <- as.numeric(APEX_sample_information$IFM15 >dicho)
###########################
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(APEX_sample_information$survival_time,APEX_sample_information$survival_indicatior)~APEX_sample_information$IFM15percent)) 
summary(coxph(Surv(APEX_sample_information$progression_time,APEX_sample_information$progression_indicator)~APEX_sample_information$IFM15percent))
#Hazard ratio of OS prediction using continuous risk scores
summary(coxph(Surv(APEX_sample_information$survival_time,APEX_sample_information$survival_indicatior)~APEX_sample_information$IFM15)) 
summary(coxph(Surv(APEX_sample_information$progression_time,APEX_sample_information$progression_indicator)~APEX_sample_information$IFM15))
####MRCIX6
MRCIX6_score_APEX <- rep(NA,264)
for (i in 1:264){
  if(APEX["203755_at",i]/APEX["216326_s_at",i] >1 | APEX["203213_at",i]/APEX["218034_at",i] >1 | APEX["200608_s_at",i]/APEX["217731_s_at",i] >1){
    MRCIX6_score_APEX[i] <- 1
  }
  else{
    MRCIX6_score_APEX[i] <- 0
  }
}
MRCIX6_score_APEX <- as.data.frame(MRCIX6_score_APEX)
rownames(MRCIX6_score_APEX) <- colnames(APEX)
APEX_sample_information$MRCIX6 <- MRCIX6_score_APEX[rownames(APEX_sample_information),]
#########################
#HR
summary(coxph(Surv(APEX_sample_information$survival_time,APEX_sample_information$survival_indicatior)~APEX_sample_information$MRCIX6)) 
summary(coxph(Surv(APEX_sample_information$progression_time,APEX_sample_information$progression_indicator)~APEX_sample_information$MRCIX6))
####MILLENNIUM100
MILLENNIUM100<- data.frame(getWeights(getClassifier("MILLENNIUM100")))
MILLENNIUM100_score_APEX <- rep(NA,264)
for (i in 1:264){
  express <- APEX[rownames(MILLENNIUM100),i]
  weight <- MILLENNIUM100$getWeights.getClassifier..MILLENNIUM100...
  MILLENNIUM100_score_APEX[i] <- sum(weight*express) 
}
MILLENNIUM100_score_APEX <- as.data.frame(MILLENNIUM100_score_APEX)
rownames(MILLENNIUM100_score_APEX) <- colnames(APEX)
APEX_sample_information$MILLENNIUM100 <- MILLENNIUM100_score_APEX[rownames(APEX_sample_information),]
#########################
#base on top 53.4% risk score
dicho <- sort(MILLENNIUM100_score_APEX$MILLENNIUM100_score_APEX,decreasing = TRUE)[round(264*0.534)]
APEX_sample_information$MILLENNIUM100percent <- as.numeric(APEX_sample_information$MILLENNIUM100 >dicho)
###########################
#Hazard ratio of OS predicition by dichotomizing risk scores to 0 and 1
summary(coxph(Surv(APEX_sample_information$survival_time,APEX_sample_information$survival_indicatior)~APEX_sample_information$MILLENNIUM100percent)) 
summary(coxph(Surv(APEX_sample_information$progression_time,APEX_sample_information$progression_indicator)~APEX_sample_information$MILLENNIUM100percent))
#Hazard ratio of OS prediction using continuous risk scores
summary(coxph(Surv(APEX_sample_information$survival_time,APEX_sample_information$survival_indicatior)~APEX_sample_information$MILLENNIUM100)) 
summary(coxph(Surv(APEX_sample_information$progression_time,APEX_sample_information$progression_indicator)~APEX_sample_information$MILLENNIUM100))
#################################################
inputdata <- apply(APEX_sample_information,2,as.numeric)
#survival
#EMC92
Est.Cval(inputdata[,c(2,1,5)],2*365.25,nofit = TRUE)[["Dhat"]]
#UAMS70
Est.Cval(inputdata[,c(2,1,6)],2*365.25,nofit = TRUE)[["Dhat"]]
#UAMS17
Est.Cval(inputdata[,c(2,1,7)],2*365.25,nofit = TRUE)[["Dhat"]]
#UAMS80
Est.Cval(inputdata[,c(2,1,8)],2*365.25,nofit = TRUE)[["Dhat"]]
#HM19
Est.Cval(inputdata[,c(2,1,9)],2*365.25,nofit = TRUE)[["Dhat"]]
#IFM15
Est.Cval(inputdata[,c(2,1,10)],2*365.25,nofit = TRUE)[["Dhat"]]
#MRCIX6
Est.Cval(inputdata[,c(2,1,11)],2*365.25,nofit = TRUE)[["Dhat"]]
#MILLENNIUM100
Est.Cval(inputdata[,c(2,1,12)],2*365.25,nofit = TRUE)[["Dhat"]]

#progression
#EMC92
Est.Cval(inputdata[,c(4,3,5)],5*365.25,nofit = TRUE)[["Dhat"]]
#UAMS70
Est.Cval(inputdata[,c(4,3,6)],5*365.25,nofit = TRUE)[["Dhat"]]
#UAMS17
Est.Cval(inputdata[,c(4,3,7)],5*365.25,nofit = TRUE)[["Dhat"]]
#UAMS80
Est.Cval(inputdata[,c(4,3,8)],5*365.25,nofit = TRUE)[["Dhat"]]
#HM19
Est.Cval(inputdata[,c(4,3,9)],5*365.25,nofit = TRUE)[["Dhat"]]
#IFM15
Est.Cval(inputdata[,c(4,3,10)],5*365.25,nofit = TRUE)[["Dhat"]]
#MRCIX6
Est.Cval(inputdata[,c(4,3,11)],5*365.25,nofit = TRUE)[["Dhat"]]
#MILLENNIUM100
Est.Cval(inputdata[,c(4,3,12)],5*365.25,nofit = TRUE)[["Dhat"]]




