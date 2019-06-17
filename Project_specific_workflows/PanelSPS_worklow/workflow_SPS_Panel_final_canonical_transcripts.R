setwd("/home/hpz440/Documents/Panel_SPS_files/TSV_only_sps")

library(gdata)
library(gtools)
library(mgsub)
library(data.table)



#Loading the samples
samples_list <- read.csv("/home/hpz440/Documents/Panel_SPS_files/samples_list_SPS_all_samples_only_sps.csv", header= FALSE, sep= ",", quote="", stringsAsFactors = F) #loading a file with the IDs in the right order
all_data <- list.files(pattern='*.tsv')
data_list <- lapply(all_data, fread, sep="\t", header=TRUE, stringsAsFactors = FALSE, data.table=FALSE) #load the data

for (i in 1:length(data_list)) {  
  data_list[[i]]$sample_ID <-NA
  data_list[[i]]$sample_ID <- as.character(samples_list[i,]) 
 }

data_compilation <- do.call("rbind", data_list) #merge all the data frames in one big data frame 
data_compilation$CHROM <- gsub("chr", " ", data_compilation$CHROM)

#export raw
#write.table(data_compilation, file="all_samples_raw.tsv", sep="\t", row.names = F, na="NA")


#Canonical transcript filter
canonical_transcripts <- read.xls("/home/hpz440/Documents/Panel_SPS/transcriptos_canonicos_panel_SPS.xlsx")
canonical_transcripts <- as.character(canonical_transcripts$Transcript_ID)
canonical_transcripts <- as.data.frame(canonical_transcripts, row.names = NULL, optional = FALSE)
colnames(canonical_transcripts)<-"ANN[*].FEATUREID"

canonical_filter <- merge(data_compilation,canonical_transcripts)
raw_data_1 <- canonical_filter[c(83,2,3,4,5,6,7,8,9,10,1,11,12,13,14:82)]

#Genes panel filter

#gene_panel_filter <- read.xls("/home/hpz440/Documents/Panel_SPS/gene_panel_list.xlsx")
#colnames(gene_panel_filter)<-"ANN[*].GENE"
#GENE_PANEL_FILTER <- merge(data_compilation,gene_panel_filter)
#raw_data_1 <- GENE_PANEL_FILTER[c(83,2,3,4,5,6,7,8,9,10,1,11,12,13,14:82)] #sorting columns


#Allelic frequency filter

names(raw_data_1)[names(raw_data_1)=="AF_raw"] <- "AF_genomAD" 

raw_data_1[, "dbNSFP_ExAC_Adj_AF"] <- as.numeric(as.character(raw_data_1[, "dbNSFP_ExAC_Adj_AF"]))
raw_data_1[, "AF_genomAD"] <- as.numeric(as.character(raw_data_1[, "AF_genomAD"]))

fq_filter <- which (
     (raw_data_1$dbNSFP_ExAC_Adj_AF <= 0.01 | is.na(raw_data_1$dbNSFP_ExAC_Adj_AF))  &
    (raw_data_1$AF_genomAD <= 0.01 | is.na(raw_data_1$AF_genomAD))
)

raw_data_2<-raw_data_1[fq_filter,]

#Coverage filter
split_AD_column <- strsplit(raw_data_2$`GEN[*].AD`,",")
numconversion_split_AD_column<-lapply(split_AD_column,as.numeric)
sum_split_AD_column<-sapply(numconversion_split_AD_column,sum)
raw_data_2[,"GEN[*].DP"]<-sum_split_AD_column  

cv_filter <- which ((raw_data_2$`GEN[*].DP`>= 10)) 
raw_data_3 <- raw_data_2[cv_filter,]


#Genotype Quality Filter
raw_data_3$`GEN[*].GQ` <- as.numeric(as.character(raw_data_3$`GEN[*].GQ`))
gq_filter <- which ((raw_data_3$`GEN[*].GQ`  >= 50)) 
raw_data_4 <- raw_data_3[gq_filter,]


#Variant Frequency Filter
vf_filter <- which(raw_data_4$`GEN[*].VF` >= 0.2)
raw_data_5 <- raw_data_4[vf_filter,]


#Variant impact filter
raw_data_5$dbNSFP_SIFT_pred_v2 <- gsub(".*D.*", "D", raw_data_5$dbNSFP_SIFT_pred)
raw_data_5$dbNSFP_LRT_pred_v2 <-gsub(".*D.*", "D", raw_data_5$dbNSFP_LRT_pred)
raw_data_5$dbNSFP_Polyphen2_HVAR_pred_v2 <- gsub(".*D.*", "D", raw_data_5$dbNSFP_Polyphen2_HVAR_pred)
raw_data_5$dbNSFP_Polyphen2_HVAR_pred_v2 <- gsub(".*P.*", "P", raw_data_5$dbNSFP_Polyphen2_HVAR_pred_v2)

raw_data_5$dbNSFP_phyloP46way_placental <- as.numeric(as.character(raw_data_5$dbNSFP_phyloP46way_placental))
raw_data_5$dbNSFP_SIFT_pred_v2 <- as.character(raw_data_5$dbNSFP_SIFT_pred_v2)
raw_data_5$dbNSFP_Polyphen2_HVAR_pred_v2 <- as.character(raw_data_5$dbNSFP_Polyphen2_HVAR_pred_v2)
raw_data_5$dbNSFP_CADD_phred <- as.numeric(as.character(raw_data_5$dbNSFP_CADD_phred))
raw_data_5$dbNSFP_LRT_pred_v2 <- as.character(raw_data_5$dbNSFP_LRT_pred_v2)
raw_data_5$dbNSFP_MutationTaster_pred <- as.character(raw_data_5$dbNSFP_MutationTaster_pred)

a<-as.numeric((raw_data_5$dbNSFP_phyloP46way_placental)>1.6)
b<-as.numeric((raw_data_5$dbNSFP_SIFT_pred_v2)=="D")
c<-as.numeric((raw_data_5$dbNSFP_Polyphen2_HVAR_pred_v2)=="D" | (raw_data_5$dbNSFP_Polyphen2_HVAR_pred_v2)=="P")
d<-as.numeric((raw_data_5$dbNSFP_CADD_phred)>15.0)
e<-as.numeric((raw_data_5$dbNSFP_LRT_pred_v2)=="D")
g<-as.numeric((raw_data_5$dbNSFP_MutationTaster_pred)=="D" | (raw_data_5$dbNSFP_MutationTaster_pred)=="A")


a[which(is.na(a)==TRUE)]<-0.25
b[which(is.na(b)==TRUE)]<-0.25
c[which(is.na(c)==TRUE)]<-0.25
d[which(is.na(d)==TRUE)]<-0.25
e[which(is.na(e)==TRUE)]<-0.25
g[which(is.na(g)==TRUE)]<-0.25

columns <- list(a, b, c, d, e, g)
raw_data_5$tools_sum <-apply(data.frame(columns),1,sum)

#write.table(raw_data_5, file="SPS_rawdata5_tool_sum.tsv", sep="\t", row.names = F, na="NA")

HIGH_rows <- which(raw_data_5$`ANN[*].IMPACT` == "HIGH")
missense_rows <- grep("missense", raw_data_5$`ANN[*].EFFECT`)


missense_filter <- raw_data_5[missense_rows,]
missense_rows_filtered <- which(missense_filter$tools_sum >=2) 

raw_data_6 <- rbind(missense_filter[missense_rows_filtered,], raw_data_5[HIGH_rows,])
raw_data_6$CHROM <- as.numeric(raw_data_6$CHROM)
raw_data_6$POS <- as.numeric(raw_data_6$POS)
raw_data_6_rows <- order(raw_data_6[,2], raw_data_6[,3])
raw_data_6.1 <- raw_data_6[raw_data_6_rows,]



#Splitting data per sample
x<-list()
for ( i in 1:nrow(samples_list)) {   
  sample_ID <- samples_list[i,]
  data_vs_ID <- which(raw_data_6.1$sample_ID == sample_ID)
  x[[i]] <- raw_data_6.1[data_vs_ID,]
}

#Exporting data
for (i in 1:nrow(samples_list)) {
  write.table(x[[i]], file= paste0(samples_list[i,], "-ID_Panel_SPS_3.0.tsv"), sep="\t", row.names = F, na="NA")
  cat("*")
}

write.table(raw_data_6.1, file="Panel_SPS_all_samples_only_SPS_3.0.tsv", sep="\t", row.names = F, na="NA")

#Allelic frequency filter #second one (<=0.001)
fq_filter <- which (
  (raw_data_6.1$dbNSFP_ExAC_Adj_AF <= 0.001 | is.na(raw_data_6.1$dbNSFP_ExAC_Adj_AF))  &
    (raw_data_6.1$AF_genomAD <= 0.001 | is.na(raw_data_6.1$AF_genomAD))
)

raw_data_7 <-raw_data_6.1[fq_filter,]

#Path prediction tools (Missense filter)  #second one (>=3)
HIGH_rows <- which(raw_data_7$`ANN[*].IMPACT` == "HIGH")
missense_rows <- grep("missense", raw_data_7$`ANN[*].EFFECT`)

missense_filter <- raw_data_7[missense_rows,]
missense_rows_filtered <- which(missense_filter$tools_sum >=3) 

raw_data_8 <- rbind(missense_filter[missense_rows_filtered,], raw_data_7[HIGH_rows,])
raw_data_8$CHROM <- as.numeric(raw_data_8$CHROM)
raw_data_8$POS <- as.numeric(raw_data_8$POS)
raw_data_8_rows <- order(raw_data_8[,2], raw_data_8[,3])
raw_data_8.1 <- raw_data_8[raw_data_8_rows,]

#extract per sample ID 

y<-list()

for ( i in 1:nrow(samples_list)) {  
  sample_ID <- samples_list[i,]
  data_vs_ID <- which(raw_data_8.1$sample_ID == sample_ID)
  y[[i]] <- raw_data_8.1[data_vs_ID,]
}

#Exporting data

for (i in 1:nrow(samples_list)) {
  write.table(y[[i]], file= paste0(samples_list[i,], "-ID_Panel_SPS_6.0.tsv"), sep="\t", row.names=F, na="NA") 
  cat("*")
}

write.table(raw_data_8.1, file="Panel_SPS_all_samples_only_SPS_6.0.tsv", sep="\t", row.names = F, na="NA")


                            

