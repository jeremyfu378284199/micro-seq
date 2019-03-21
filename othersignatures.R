load("~/downloads/DFCI_TPM_Survival.Rd")
row.names(DFCI_TPM_Gene) <- gsub("\\.\\d+","",row.names(DFCI_TPM_Gene))
source("http://bioconductor.org/biocLite.R")
biocLite("geneClassifiers")
biocLite("biomaRt")
library("geneClassifiers")
library("biomaRt")
library(survival)
library(survC1)
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",  host = 'www.ensembl.org',ensemblRedirect = FALSE)
##################################################
#EMC92 gene signature  (prebs &mmrf)
EMC92 <- data.frame(getWeights(getClassifier("EMC92")))
E92 <- rownames(EMC92)
#get the EMC92 matching information
G_EMC92 <- getBM(attributes = c("affy_hg_u133_plus_2", "hgnc_symbol","ensembl_gene_id"),filters = "affy_hg_u133_plus_2",values = E92,mart = mart)
common <- intersect(G_EMC92$ensembl_gene_id,rownames(DFCI_TPM_Gene))
G_EMC92 <- G_EMC92[G_EMC92$ensembl_gene_id %in% common,]
write.csv(G_EMC92,file = "EMC92_prebs.csv")
#################################################
#UAMS70 gene signature
UAMS70 <- data.frame(getWeights(getClassifier("UAMS70")))
U70 <- rownames(UAMS70)
#get the UAMS70 matching information
G_UAMS70 <- getBM(attributes = c("affy_hg_u133_plus_2", "hgnc_symbol","ensembl_gene_id"),filters = "affy_hg_u133_plus_2",values = U70,mart = mart)#use unique #68
common <- intersect(G_UAMS70$ensembl_gene_id,rownames(DFCI_TPM_Gene))#81
G_UAMS70 <- G_UAMS70[G_UAMS70$ensembl_gene_id %in% common,]
#write.csv(G_UAMS70,file = "UAMS70_prebs.csv")
prebs <- read.table("allexpress2.csv",sep = ",",row.names = 1,header = TRUE,stringsAsFactors = FALSE)
colnames(prebs) <- gsub("X","",colnames(prebs))
commonsampleids <- intersect(colnames(prebs),colnames(DFCI_TPM_Gene)) #309
########################################
#caculate the correlation between one probeset and multiple ensembl ids
cor(as.numeric(prebs["200750_s_at",commonsampleids]),as.numeric(DFCI_TPM_Gene["ENSG00000132341",commonsampleids]))
############################
#assign weights
diff <- setdiff(rownames(UAMS70),G_UAMS70$affy_hg_u133_plus_2)
diff <- c(diff,"1557277_a_at")
diff <- c("237964_at","1557277_a_at")
#"244686_at" "237964_at"
which(rownames(UAMS70) =="1557277_a_at")#sd=0
UAMS68 <- data.frame(UAMS70[!(rownames(UAMS70) %in% diff),])
rownames(UAMS68) <- rownames(UAMS70)[-c(49,61,68)]
rownames(UAMS68) <- rownames(UAMS70)[-c(61,68)]
UAMS68$UAMS70...rownames.UAMS70...in..diff.... <- c(rep(1/51,51),rep(-1/17,17))
U70 <- read.table("UAMS70_prebs.csv",sep = ",",row.names = 1,header = TRUE,stringsAsFactors = FALSE)
which(U70$probe =="1557277_a_at")
U70 <- U70[-23,]
probe <- U70$probe
U70$weight <- UAMS68[probe,]                            
U70 <- U70[,-4]        
colnames(U70) <- c("x","probe","gene","id","weight")
write.csv(U70,file = "UAMS70_prebs2.csv")
U70 <- read.table("UAMS70_prebs.csv",sep = ",",header = TRUE,stringsAsFactors = FALSE)

#############################################
#####################
#UAMS80    (prebs &mmrf)
UAMS80 <- data.frame(getWeights(getClassifier("UAMS80")))
U80 <- rownames(UAMS80)
#get the UAMS70 matching information
G_UAMS80 <- getBM(attributes = c("affy_hg_u133_plus_2", "hgnc_symbol","ensembl_gene_id"),filters = "affy_hg_u133_plus_2",values = U80,mart = mart)#use unique #71
common1 <- intersect(G_UAMS80$ensembl_gene_id,rownames(DFCI_TPM_Gene))#68
common <- intersect(G_UAMS80$ensembl_gene_id,rownames(MMRF))#66
G_UAMS80 <- G_UAMS80[G_UAMS80$ensembl_gene_id %in% common,]
write.csv(G_UAMS80,file = "UAMS80_mmrf.csv")
prebs <- read.table("allexpress2.csv",sep = ",",row.names = 1,header = TRUE,stringsAsFactors = FALSE)
colnames(prebs) <- gsub("X","",colnames(prebs))
commonsampleids <- intersect(colnames(prebs),colnames(DFCI_TPM_Gene)) #309
########################################
#caculate the correlation between one probeset and multiple ensembl ids
cor(as.numeric(prebs["200882_s_at",commonsampleids]),as.numeric(DFCI_TPM_Gene["ENSG00000159352",commonsampleids]))
############################
#assign weights
diff <- setdiff(rownames(UAMS80),G_UAMS80$affy_hg_u133_plus_2)
#37,38,40,/44,48,49,56,64,69 for prebs
#14,37,38,40,42,/44,48,49,56,64,69 for mmrf
index <- rep(NA,11)
for (i in 1:11){
  index[i] <- which(rownames(UAMS80) == diff[i])
}
UAMS71 <- data.frame(UAMS80[!(rownames(UAMS80) %in% diff),])
rownames(UAMS71) <- rownames(UAMS80)[-index]
UAMS71$UAMS80...rownames.UAMS80...in..diff.... <- c(rep(1/37,37),rep(-1/32,32))
U80 <- read.table("UAMS80_mmrf.csv",sep = ",",row.names = 1,header = TRUE,stringsAsFactors = FALSE)
#which(U70$probe =="1557277_a_at")#sd = 0
#U70 <- U70[-23,]
probe <- U80$affy_hg_u133_plus_2
U80$weight <- UAMS71[probe,]                            
colnames(U80) <- c("probe","gene","id","weight")
write.csv(U80,file = "UAMS80_mmrf.csv")
UAMS79 <- data.frame(UAMS80[-38,])
rownames(UAMS79) <- rownames(UAMS80)[-38]
UAMS79$UAMS80..38... <- c(rep(1/41,41),rep(-1/38,38))
##########################################
###############################
#####################
#UAMS17 
UAMS17 <- data.frame(getWeights(getClassifier("UAMS17")))
U17 <- rownames(UAMS17)
#get the UAMS70 matching information
G_UAMS17 <- getBM(attributes = c("affy_hg_u133_plus_2", "hgnc_symbol","ensembl_gene_id"),filters = "affy_hg_u133_plus_2",values = U17,mart = mart)#use unique #17
common <- intersect(G_UAMS17$ensembl_gene_id,rownames(DFCI_TPM_Gene))#19
G_UAMS17 <- G_UAMS17[G_UAMS17$ensembl_gene_id %in% common,]
write.csv(G_UAMS17,file = "UAMS17_prebs.csv")
prebs <- read.table("allexpress2.csv",sep = ",",row.names = 1,header = TRUE,stringsAsFactors = FALSE)
colnames(prebs) <- gsub("X","",colnames(prebs))
commonsampleids <- intersect(colnames(prebs),colnames(DFCI_TPM_Gene)) #309
########################################
#caculate the correlation between one probeset and multiple ensembl ids
cor(as.numeric(prebs["201897_s_at",commonsampleids]),as.numeric(DFCI_TPM_Gene["ENSG00000173207",commonsampleids]))
############################
#assign weights
diff <- setdiff(rownames(UAMS80),G_UAMS80$affy_hg_u133_plus_2)
#37,38,40,/44,48,49,56,64,69
index <- rep(NA,9)
for (i in 1:9){
  index[i] <- which(rownames(UAMS80) == diff[i])
}
UAMS71 <- data.frame(UAMS80[!(rownames(UAMS80) %in% diff),])
rownames(UAMS71) <- rownames(UAMS80)[-index]
UAMS71$UAMS80...rownames.UAMS80...in..diff.... <- c(rep(1/39,39),rep(-1/32,32))
U17 <- read.table("UAMS17_prebs.csv",sep = ",",row.names = 1,header = TRUE,stringsAsFactors = FALSE)
#which(U70$probe =="1557277_a_at")#sd = 0
#U70 <- U70[-23,]
probe <- U17$affy_hg_u133_plus_2
U17$weight <- UAMS17[probe,]                            
colnames(U17) <- c("probe","gene","id","weight")
write.csv(U17,file = "UAMS17_prebs.csv")
##########################################
###############################
#####################
#HM19
HM19 <- data.frame(getWeights(getClassifier("HM19")))
H19 <- rownames(HM19)
#get the UAMS70 matching information
G_HM19 <- getBM(attributes = c("affy_hg_u133_plus_2", "hgnc_symbol","ensembl_gene_id"),filters = "affy_hg_u133_plus_2",values = H19,mart = mart)#use unique #18
common <- intersect(G_HM19$ensembl_gene_id,rownames(DFCI_TPM_Gene))#19
write.csv(G_HM19,file = "HM19_prebs.csv")
prebs <- read.table("allexpress2.csv",sep = ",",row.names = 1,header = TRUE,stringsAsFactors = FALSE)
colnames(prebs) <- gsub("X","",colnames(prebs))
commonsampleids <- intersect(colnames(prebs),colnames(DFCI_TPM_Gene)) #309
############################
#assign weights
diff <- setdiff(rownames(HM19),G_HM19$affy_hg_u133_plus_2)
which(rownames(HM19) == diff)#5

HM18 <- data.frame(HM19[!(rownames(HM19) %in% diff),])
rownames(HM18) <- rownames(HM19)[-5]
H19 <- read.table("HM19_prebs.csv",sep = ",",row.names = 1,header = TRUE,stringsAsFactors = FALSE)
#which(U70$probe =="1557277_a_at")#sd = 0
#U70 <- U70[-23,]
probe <- H19$affy_hg_u133_plus_2
H19$weight <- HM18[probe,]                            
colnames(H19) <- c("probe","gene","id","weight")
write.csv(H19,file = "HM19_prebs.csv")
##########################################
###############################
#####################
#IFM15
IFM15 <- data.frame(getWeights(getClassifier("IFM15")))
I15 <- rownames(IFM15)
#get the UAMS70 matching information
G_IFM15 <- getBM(attributes = c("affy_hg_u133_plus_2", "hgnc_symbol","ensembl_gene_id"),filters = "affy_hg_u133_plus_2",values = I15,mart = mart)#use unique #16
common <- intersect(G_IFM15$ensembl_gene_id,rownames(DFCI_TPM_Gene))#16
G_IFM15 <- G_IFM15[G_IFM15$ensembl_gene_id %in% common,]
write.csv(G_IFM15,file = "IFM15_prebs.csv")
prebs <- read.table("allexpress2.csv",sep = ",",row.names = 1,header = TRUE,stringsAsFactors = FALSE)
colnames(prebs) <- gsub("X","",colnames(prebs))
commonsampleids <- intersect(colnames(prebs),colnames(DFCI_TPM_Gene)) #309
########################################
#caculate the correlation between one probeset and multiple ensembl ids
cor(as.numeric(prebs["200779_at",commonsampleids]),as.numeric(DFCI_TPM_Gene["ENSG00000128272",commonsampleids]))
############################
#assign weights
I15 <- read.table("IFM15_prebs.csv",sep = ",",row.names = 1,header = TRUE,stringsAsFactors = FALSE)
#which(U70$probe =="1557277_a_at")#sd = 0
#U70 <- U70[-23,]
probe <- I15$affy_hg_u133_plus_2
I15$weight <- IFM15[probe,]                            
colnames(I15) <- c("probe","gene","id","weight")
write.csv(I15,file = "IFM15_prebs.csv")
##########################################
#####################
#MRCIX6
MRCIX6 <- data.frame(getWeights(getClassifier("MRCIX6")))
M6 <- rownames(MRCIX6)
G_MRCIX6 <- getBM(attributes = c("affy_hg_u133_plus_2", "hgnc_symbol","ensembl_gene_id"),filters = "affy_hg_u133_plus_2",values = M6,mart = mart)#use unique #6
common <- intersect(G_MRCIX6$ensembl_gene_id,rownames(DFCI_TPM_Gene))
write.csv(G_MRCIX6,file = "MRCIX6_prebs.csv")
##########################################
#MILLWNNIUM100 gene signature  (prebs &mmrf)
MILLENNIUM100 <- data.frame(getWeights(getClassifier("MILLENNIUM100")))
M100 <- rownames(MILLENNIUM100)
#get the MILLWNNIUM100 matching information
G_M100 <- getBM(attributes = c("affy_hg_u133_plus_2", "hgnc_symbol","ensembl_gene_id"),filters = "affy_hg_u133_plus_2",values = M100,mart = mart)#use unique #99
common <- intersect(G_M100$ensembl_gene_id,rownames(DFCI_TPM_Gene))#101 ensembl ids
#check common ensembl id numbers for mmrf
common <- intersect(G_M100$ensembl_gene_id,rownames(MMRF))#95 ensembl ids
G_M100 <- G_M100[G_M100$ensembl_gene_id %in% common,]
write.csv(G_M100,file = "M100_mmrf.csv")
prebs <- read.table("allexpress2.csv",sep = ",",row.names = 1,header = TRUE,stringsAsFactors = FALSE)
colnames(prebs) <- gsub("X","",colnames(prebs))
commonsampleids <- intersect(colnames(prebs),colnames(DFCI_TPM_Gene)) #309
########################################
#caculate the correlation between one probeset and multiple ensembl ids
cor(as.numeric(prebs["211674_x_at",commonsampleids]),as.numeric(DFCI_TPM_Gene["ENSG00000126890",commonsampleids]))
############################
#assign weights
diff <- setdiff(rownames(MILLENNIUM100),G_M100$affy_hg_u133_plus_2) #215811_at
#diff <- c("237964_at","1557277_a_at")
#"244686_at" "237964_at"
which(rownames(MILLENNIUM100) =="215811_at")#49
M99 <- data.frame(MILLENNIUM100[!(rownames(MILLENNIUM100) %in% diff),])
rownames(M99) <- rownames(MILLENNIUM100)[-49]
M100 <- read.table("M100_mmrf.csv",sep = ",",row.names = 1,header = TRUE,stringsAsFactors = FALSE)
M100$weight <- M99[M100$affy_hg_u133_plus_2,]               
colnames(M100) <- c("probe","gene","id","weight")
write.csv(M100,file = "M100_mmrf.csv")
