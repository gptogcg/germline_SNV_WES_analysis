library(gdata)
library(gtools)

### Pongo mis comentarios después de 3 # ok? (Marcos)


#Load the tsv file
#raw_data<-read.table("AA3542.PON_SPS_39samples.vcf_PASS_snpEff_dbNSFP_genomad.tsv",sep="\t",header=T,stringsAsFactors = F)

samples_list <- read.csv("sampleslist.csv", header= FALSE, sep= ",", quote="", stringsAsFactors = F) #loading a file with the IDs in the right order


all_data <- list.files(pattern='*.tsv')
data_list <- lapply(all_data,read.table,sep="\t",header=T, stringsAsFactors = F) #load the data

for (i in 1:length(data_list)) {  
  data_list[[i]]$sample_ID <-NA
  data_list[[i]]$sample_ID <- as.character(samples_list[i,]) #attribute each element from the sample_list to my data_list is tricky here, it was a coincidence it did it in order, but might not happen always
  
  ###It isn't so tricky. You just have to order in the same way all_data vector and sample_list dataframe you loaded. If the name of the file in all_data starts with the sample name you put in sampleslist.csv file everything perfect!! Really well done!
  
#  data_list[[i]]$sample_ID <- as.character(data_list[[i]]$sample_ID) 
  ### It's easier to combine both commands -> see the command before
 }

data_compilation <- do.call("rbind", data_list) #merge all the data frames in one big data frame 

#Canonical transcript filter
canonical_transcripts <- read.xls("/home/hpz440/Documents/Panel_SPS/canonical_transcripts_GRCh37.xls")
canonical_transcripts <- as.character(canonical_transcripts$Transcript.stable.ID)
canonical_transcripts <- as.data.frame(canonical_transcripts, row.names = NULL, optional = FALSE)
colnames(canonical_transcripts)<-"ANN....FEATUREID"

raw_data_1 <- merge(data_compilation,canonical_transcripts) #mucho más rápido

###Ojo que al hacer el merge te queda la columna de ANN....FEATUREID como la primera. Si prefieras mantenerla donde estaba inicialmente tienes que cambiarla de sitio, acuérdate. ;) Y lo mismo pasa con el orden de las filas. Antes estaban ordenadas por posición y cromosoma, ahora están ordenadas por el valor de ANN....FEATUREID. Lo mismo, no es que haya nada mal pero ten en cuenta como tienes los datos ahora y que sepas el porqué.


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
split_AD_column <- strsplit(raw_data_2$GEN.NORMAL..AD,",")
numconversion_split_AD_column<-lapply(split_AD_column,as.numeric)
sum_split_AD_column<-sapply(numconversion_split_AD_column,sum)
raw_data_2[,"GEN.NORMAL..DP"]<-sum_split_AD_column  #I filtered using only the NORMAL field ( didn't want to make it more complicated, I'll probably just have one)

cv_filter <- which ((raw_data_2$GEN.NORMAL..DP >= 10)) 
raw_data_3 <- raw_data_2[cv_filter,]


#Genotype Quality Filter
#raw_data_3[, "GEN.D3204.PBL9..GQ"] <- as.numeric(as.character(raw_data_3[,"GEN.D3204.PBL9..GQ"]))
#gq_filter <- which ((raw_data_3$GEN.D3204.PBL9..GQ  >= 50)) 
#raw_data_4 <- raw_data_3[gq_filter,]
raw_data_4 <- raw_data_3

#Variant impact filter
raw_data_5 <- raw_data_4
raw_data_5$dbNSFP_phyloP46way_placental <- as.numeric(as.character(raw_data_5$dbNSFP_phyloP46way_placental))
raw_data_5$dbNSFP_SIFT_pred <- as.character(raw_data_5$dbNSFP_SIFT_pred)
raw_data_5$dbNSFP_Polyphen2_HVAR_pred <- as.character(raw_data_5$dbNSFP_Polyphen2_HVAR_pred)
raw_data_5$dbNSFP_CADD_phred <- as.numeric(as.character(raw_data_5$dbNSFP_CADD_phred))
raw_data_5$dbNSFP_LRT_pred <- as.character(raw_data_5$dbNSFP_LRT_pred)
raw_data_5$dbNSFP_MutationTaster_pred <- as.character(raw_data_5$dbNSFP_MutationTaster_pred)

a<-as.numeric((raw_data_5$dbNSFP_phyloP46way_placental)>1.6)
b<-as.numeric((raw_data_5$dbNSFP_SIFT_pred)=="D")
c<-as.numeric((raw_data_5$dbNSFP_Polyphen2_HVAR_pred)=="D" | (raw_data_5$dbNSFP_Polyphen2_HVAR_pred)=="P")
d<-as.numeric((raw_data_5$dbNSFP_CADD_phred)>15.0)
e<-as.numeric((raw_data_5$dbNSFP_LRT_pred)=="D")
g<-as.numeric((raw_data_5$dbNSFP_MutationTaster_pred)=="D" | (raw_data_5$dbNSFP_MutationTaster_pred)=="A")

a[which(is.na(a)==TRUE)]<-0.25
b[which(is.na(b)==TRUE)]<-0.25
c[which(is.na(c)==TRUE)]<-0.25
d[which(is.na(d)==TRUE)]<-0.25
e[which(is.na(e)==TRUE)]<-0.25
g[which(is.na(g)==TRUE)]<-0.25

columns <- list(a, b, c, d, e, g)
raw_data_5$tools_sum <-apply(data.frame(columns),1,sum)

HIGH_rows <- which(raw_data_5$ANN....IMPACT == "HIGH")
missense_rows <- grep("missense", raw_data_5$ANN....EFFECT)

missense_filter <- raw_data_5[missense_rows,]
missense_rows_filtered <- which(missense_filter$tools_sum >=2) 

raw_data_6 <- rbind(missense_filter[missense_rows_filtered,], raw_data_5[HIGH_rows,])


#extract per sample ID #want to save each loop result in differents files

#x <- vector("list", length(samples_list))  #no funciona bien, esperaba que hubiese un length de tres

### I don't understand what you're trying to do here. I mean, you just have to initialize x list before starting to save things in it, but this can be done as follows:

x<-list()


for ( i in 1:nrow(samples_list)) {   #mismo con el problema de mi vector vacio tener solo un elemento, aquí terminó por funcionar
  sample_ID <- samples_list[i,]
  data_vs_ID <- which(raw_data_6$sample_ID == sample_ID)
  x[[i]] <- raw_data_6[data_vs_ID,]
}

#Exporting data
#WriteXLS(raw_data_6, ExcelFileName = "../Panel_SPS/Panel_SPS_3.0.xls" , perl = "perl", na ="NA" )
for (i in 1:nrow(samples_list)) {
  #write.table(x[[i]], file= paste0(samples_list[i,], "-ID_Panel_SPS_3.0.tsv"), sep="\t", col.names = NA, na="NA")
  
  ###The important option is to avoid print the row names, not the column names
  
  write.table(x[[i]], file= paste0(samples_list[i,], "-ID_Panel_SPS_3.0.tsv"), sep="\t", row.names = F, na="NA")
  cat("*") ###And this? Do you like *? haha. It's cool!
}

#Allelic frequency filter #second one (<=0.001)
fq_filter <- which (
  (raw_data_6$dbNSFP_ExAC_Adj_AF <= 0.001 | is.na(raw_data_6$dbNSFP_ExAC_Adj_AF))  &
    (raw_data_6$AF_genomAD <= 0.001 | is.na(raw_data_6$AF_genomAD))
)

raw_data_7 <-raw_data_6[fq_filter,]


#Path prediction tools (Missense filter)  #second one (>=3)
HIGH_rows <- which(raw_data_7$ANN....IMPACT == "HIGH")
missense_rows <- grep("missense", raw_data_7$ANN....EFFECT)

missense_filter <- raw_data_7[missense_rows,]
missense_rows_filtered <- which(missense_filter$tools_sum >=3) 

raw_data_8 <- rbind(missense_filter[missense_rows_filtered,], raw_data_7[HIGH_rows,])

#extract per sample ID 

#y <- vector("list", length(samples_list)) 
y<-list()

for ( i in 1:nrow(samples_list)) {  
  sample_ID <- samples_list[i,]
  data_vs_ID <- which(raw_data_8$sample_ID == sample_ID)
  y[[i]] <- raw_data_8[data_vs_ID,]
}

#Exporting data

for (i in 1:nrow(samples_list)) {
  #write.table(x[[i]], file= paste0(samples_list[i,], "-ID_Panel_SPS_6.0.tsv"), sep="\t", col.names = NA, na="NA")
  write.table(y[[i]], file= paste0(samples_list[i,], "-ID_Panel_SPS_6.0.tsv"), sep="\t", row.names=F, na="NA") ###Ojo que no habías cambiado la x por la y!
  cat("*")
}


#Exporting data
#WriteXLS(raw_data_8, ExcelFileName = "../Panel_SPS/Panel_SPS_6.0.xls" )
