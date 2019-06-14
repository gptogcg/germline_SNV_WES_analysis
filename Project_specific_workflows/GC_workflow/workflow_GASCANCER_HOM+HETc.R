library(WriteXLS)
library(gdata)
library(data.table)

dir_input <- ""
dir.create("RESULTS_RECESSIVE")
dir_output <- "./RESULTS_RECESSIVE/"

##################      Loading files needed      ##################      
SAMPLES <- read.table("./LIST_samples_GASCANCER_01.txt", sep="\t", header=T)
input_files <- read.table("files_to_analyze.txt", header=T, sep="\t")


FILE <- as.character(input_files[1,1])
NAME <- as.character(input_files[1,2])
TAB <- fread(paste(dir_input, FILE, sep=""), header=TRUE, sep="\t", dec=".", data.table=FALSE)
large_TAB <- TAB

colnames(TAB)[1]<-"CHROM"
####################################################################


################    NOM?S GERMINALS & GC    ########################
germinals <- which(SAMPLES[,"MOSTRA"]=="germinal" & SAMPLES[,"ESTUDI"]=="GC")
SAMPLES_2 <- SAMPLES[germinals,]

TAB2 <- TAB
rownames(SAMPLES_2) <- 1:dim(SAMPLES_2)[1]
####################################################################


####################### 	CONTROL INDIVIDUS (start)	########################
samples_in_TAB2 <- gsub("_GT", "", colnames(TAB2)[grep("_GT", colnames(TAB))])

list_presents <- list()
for (j in 1:length(SAMPLES_2[,1])){
	ind <- SAMPLES_2[,1][j]
	fam <- SAMPLES_2[,4][j]
	if (length(grep(ind, samples_in_TAB2))==0){
		message <- paste("NO! : ", ind, ", from family ", fam, ", is not here!!", sep="")
		print(message)
		list_presents[[j]] <- paste(message)
	} else {
		list_presents[[j]] <- paste(ind, fam, "YES", sep="\t")
	}
}
list_extra_samples <- list()
for (i in 1:length(samples_in_TAB2)){
	ind <- samples_in_TAB2[i]
	if (length(grep(ind, SAMPLES_2[,1]))==0){
		message_1 <- paste("HEY!, ", ind, ", from family ", fam, ", is here, but it is not in the lists!!", sep="")
		print(message_1)
## columns from the missing individual will be removed ... ##
		ind_cols <- grep(ind, colnames(TAB2))
		TAB2 <- TAB2[,-ind_cols]
		message_2 <- (paste(ind, ", from family ", fam, ", has been eliminated from analysis."))
		print(message_2)
		##############################

		list_extra_samples[[i]] <- paste(message_1)
	}
}
 PRESENTS <- do.call(rbind, list_presents)
 NAME2 <- paste(dir_output, NAME, "__samples_present.txt", sep="") 
 write.table(PRESENTS, NAME2, quote=F, row.names=F, col.names=F, sep="\t")
 EXTRA_samples <- do.call(rbind, list_extra_samples)
 NAME2 <- paste(dir_output, NAME, "__extra_samples.txt", sep="") 
 write.table(EXTRA_samples, NAME2, quote=F, row.names=F, col.names=F, sep="\t")

####################### 	CONTROL INDIVIDUS (end)		########################




########################################################################################
##########|                     	        [1.0]	                        |###############
######################################### START ########################################

EXOMES_1.0 <- TAB2

########################################################################################
###########|			                      	[1.0]				                   |##############
######################################## END ###########################################


########################################################################################
##########|                      [2.0 & 3.0] primers filtres                  |###############
######################################### START ########################################
first_exome_column <- as.numeric(grep("_GT", colnames(EXOMES_1.0))[1])
nonexome_columns <- 1:(first_exome_column-1)


#HIGH/missense filter
HIGH_filter <- grep("HIGH", EXOMES_1.0$Annotation_Impact)
MISSENSE_filter <- grep("missense", EXOMES_1.0$Annotation)

EXOMES_1.0.2 <- EXOMES_1.0[c(HIGH_filter, MISSENSE_filter),]


##Inner frequency filter
cols_GT <- grep("_GT", colnames(EXOMES_1.0.2))
total <- as.numeric(length(cols_GT))
filas_bucle<-dim(EXOMES_1.0.2)[1]

for (fq in 1:filas_bucle){
  print(paste("Row ", fq, " from ", filas_bucle, sep=""))
  ###frequency annotation###
  row <- as.character(EXOMES_1.0.2[fq,cols_GT])
  relative <- as.numeric(length(which(row!="./." & row!="0/0")))
  FQ <- relative/total
  EXOMES_1.0.2$inner_frequency[fq] <- as.numeric(FQ)#fq annotation in every row
  ###frequency annotation###
}

EXOMES_1.0.2[, "inner_frequency"] <- as.numeric(EXOMES_1.0.2[,"inner_frequency"])


inner_frequency_filter <- which(EXOMES_1.0.2$inner_frequency <=0.25) 
EXOMES_1.0.X <- EXOMES_1.0.2[inner_frequency_filter,] #filtered for inner_frequency
###########################  inner freq FILTER (end) ###

EXOMES_1.0.2<-EXOMES_1.0.X







#Coverage and Genotype Quality filter
first_exome_column <- as.numeric(grep("_GT", colnames(EXOMES_1.0.2))[1])
nonexome_columns <- 1:(first_exome_column-1)

n_cols_x_ind <- as.numeric(grep("_GT", colnames(EXOMES_1.0.2))[2])-first_exome_column

all_fams <- list()
for (ss in 1:length(unique(SAMPLES_2[,"FAMILIES"]))){
  fam <- as.character(unique(SAMPLES_2[,"FAMILIES"])[ss])
  inds <- as.character(SAMPLES_2[which(SAMPLES_2[,"FAMILIES"]==fam),1])
  
  
  if (length(inds)>6){
    print("Nom?s est? preparat per fam?les amb un m?xim de 6 persones!!!!")
  } else if (length(inds)<=6){
    all_cols <- list()
    for (xx in 1:length(inds)){
      name <- inds[xx]
      cols <- grep(name, colnames(EXOMES_1.0.2))
      all_cols[[name]] <- cols
    }
    all_cols2 <- do.call(cbind, all_cols)
    INDS_COLS <- as.vector(all_cols2)
    FAM <- EXOMES_1.0.2[,c(nonexome_columns, INDS_COLS)]
    first_exome_column <- as.numeric(grep("_GT", colnames(FAM))[1])
    nonexome_columns <- 1:(first_exome_column-1)
    
    first_GT <- first_exome_column
    first_DP <- as.numeric(grep("_DP", colnames(FAM))[1])
    first_GQ <- as.numeric(grep("_GQ", colnames(FAM))[1])
    
    if (length(inds)==2){
      col_GT_1 <- first_GT
      col_GT_2 <- first_GT + n_cols_x_ind
      col_DP_1 <- first_DP
      col_DP_2 <- first_DP + n_cols_x_ind
      col_GQ_1 <- first_GQ
      col_GQ_2 <- first_GQ + n_cols_x_ind
      
      
      sel <- which( ((FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0")) & (as.numeric(as.character(FAM[,col_DP_1]))>=10 & as.numeric(as.character(FAM[,col_DP_2]))>=10) & (as.numeric(as.character(FAM[,col_GQ_1]))>=50 & as.numeric(as.character(FAM[,col_GQ_2]))>=50)  )
      # _GT and _DP filters!
      
      
    } else if (length(inds)==3){
      col_GT_1 <- first_GT
      col_GT_2 <- first_GT + n_cols_x_ind
      col_GT_3 <- first_GT + (n_cols_x_ind*2)
      col_DP_1 <- first_DP
      col_DP_2 <- first_DP + n_cols_x_ind
      col_DP_3 <- first_DP + (n_cols_x_ind*2)
      col_GQ_1 <- first_GQ
      col_GQ_2 <- first_GQ + n_cols_x_ind
      col_GQ_3 <- first_GQ + (n_cols_x_ind*2)
      
      sel <- which( ((FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0") & (FAM[,col_GT_3]!="./." & FAM[,col_GT_3]!="0/0")) & as.numeric(as.character((FAM[,col_DP_1]))>=10 & as.numeric(as.character(FAM[,col_DP_2]))>=10 & as.numeric(as.character(FAM[,col_DP_3]))>=10) & as.numeric(as.character((FAM[,col_GQ_1]))>=50 & as.numeric(as.character(FAM[,col_GQ_2]))>=50 & as.numeric(as.character(FAM[,col_GQ_3]))>=50) )
      # _GT and _DP filters!
      
      
    } else if (length(inds)==4){
      col_GT_1 <- first_GT
      col_GT_2 <- first_GT + n_cols_x_ind
      col_GT_3 <- first_GT + (n_cols_x_ind*2)
      col_GT_4 <- first_GT + (n_cols_x_ind*3)
      col_DP_1 <- first_DP
      col_DP_2 <- first_DP + n_cols_x_ind
      col_DP_3 <- first_DP + (n_cols_x_ind*2)
      col_DP_4 <- first_DP + (n_cols_x_ind*3)
      col_GQ_1 <- first_GQ
      col_GQ_2 <- first_GQ + n_cols_x_ind
      col_GQ_3 <- first_GQ + (n_cols_x_ind*2)
      col_GQ_4 <- first_GQ + (n_cols_x_ind*3)
      
      
      sel <- which( ((FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0") & (FAM[,col_GT_3]!="./." & FAM[,col_GT_3]!="0/0") & (FAM[,col_GT_4]!="./." & FAM[,col_GT_4]!="0/0")) & (as.numeric(as.character(FAM[,col_DP_1]))>=10 & as.numeric(as.character(FAM[,col_DP_2]))>=10 & as.numeric(as.character(FAM[,col_DP_3]))>=10 & as.numeric(as.character(FAM[,col_DP_4]))>=10) & (as.numeric(as.character(FAM[,col_GQ_1]))>=50 & as.numeric(as.character(FAM[,col_GQ_2]))>=50 & as.numeric(as.character(FAM[,col_GQ_3]))>=50 & as.numeric(as.character(FAM[,col_GQ_4]))>=50) )
      # _GT and _DP filters!
      
      
    } else if (length(inds)==5){
      col_GT_1 <- first_GT
      col_GT_2 <- first_GT + n_cols_x_ind
      col_GT_3 <- first_GT + (n_cols_x_ind*2)
      col_GT_4 <- first_GT + (n_cols_x_ind*3)
      col_GT_5 <- first_GT + (n_cols_x_ind*4)
      col_DP_1 <- first_DP
      col_DP_2 <- first_DP + n_cols_x_ind
      col_DP_3 <- first_DP + (n_cols_x_ind*2)
      col_DP_4 <- first_DP + (n_cols_x_ind*3)
      col_DP_5 <- first_DP + (n_cols_x_ind*4)
      col_GQ_1 <- first_GQ
      col_GQ_2 <- first_GQ + n_cols_x_ind
      col_GQ_3 <- first_GQ + (n_cols_x_ind*2)
      col_GQ_4 <- first_GQ + (n_cols_x_ind*3)
      col_GQ_5 <- first_GQ + (n_cols_x_ind*4)
      
      
      sel <- which( ((FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0") & (FAM[,col_GT_3]!="./." & FAM[,col_GT_3]!="0/0") & (FAM[,col_GT_4]!="./." & FAM[,col_GT_4]!="0/0") & (FAM[,col_GT_5]!="./." & FAM[,col_GT_5]!="0/0")) & (as.numeric(as.character(FAM[,col_DP_1]))>=10 & as.numeric(as.character(FAM[,col_DP_2]))>=10 & as.numeric(as.character(FAM[,col_DP_3]))>=10 & as.numeric(as.character(FAM[,col_DP_4]))>=10 & as.numeric(as.character(FAM[,col_DP_5]))>=10) & (as.numeric(as.character(FAM[,col_GQ_1]))>=50 & as.numeric(as.character(FAM[,col_GQ_2]))>=50 & as.numeric(as.character(FAM[,col_GQ_3]))>=50 & as.numeric(as.character(FAM[,col_GQ_4]))>=50 & as.numeric(as.character(FAM[,col_GQ_5]) )>=50))
      # _GT and _DP filters!
      
      
    } else if (length(inds)==6){
      col_GT_1 <- first_GT
      col_GT_2 <- first_GT + n_cols_x_ind
      col_GT_3 <- first_GT + (n_cols_x_ind*2)
      col_GT_4 <- first_GT + (n_cols_x_ind*3)
      col_GT_5 <- first_GT + (n_cols_x_ind*4)
      col_GT_6 <- first_GT + (n_cols_x_ind*5)
      col_DP_1 <- first_DP
      col_DP_2 <- first_DP + n_cols_x_ind
      col_DP_3 <- first_DP + (n_cols_x_ind*2)
      col_DP_4 <- first_DP + (n_cols_x_ind*3)
      col_DP_5 <- first_DP + (n_cols_x_ind*4)
      col_DP_6 <- first_DP + (n_cols_x_ind*5)
      col_GQ_1 <- first_GQ
      col_GQ_2 <- first_GQ + n_cols_x_ind
      col_GQ_3 <- first_GQ + (n_cols_x_ind*2)
      col_GQ_4 <- first_GQ + (n_cols_x_ind*3)
      col_GQ_5 <- first_GQ + (n_cols_x_ind*4)
      col_GQ_6 <- first_GQ + (n_cols_x_ind*5)
      
      
      sel <- which( ((FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0") & (FAM[,col_GT_3]!="./." & FAM[,col_GT_3]!="0/0") & (FAM[,col_GT_4]!="./." & FAM[,col_GT_4]!="0/0") & (FAM[,col_GT_5]!="./." & FAM[,col_GT_5]!="0/0") & (FAM[,col_GT_6]!="./." & FAM[,col_GT_6]!="0/0")) & (as.numeric(as.character(FAM[,col_DP_1]))>=10 & as.numeric(as.character(FAM[,col_DP_2]))>=10 & as.numeric(as.character(FAM[,col_DP_3]))>=10 & as.numeric(as.character(FAM[,col_DP_4]))>=10 & as.numeric(as.character(FAM[,col_DP_5]))>=10 & as.numeric(as.character(FAM[,col_DP_6]))>=10) & (as.numeric(as.character(FAM[,col_GQ_1]))>=50 & as.numeric(as.character(FAM[,col_GQ_2]))>=50 & as.numeric(as.character(FAM[,col_GQ_3]))>=50 & as.numeric(as.character(FAM[,col_GQ_4]))>=50 & as.numeric(as.character(FAM[,col_GQ_5]))>=50 & as.numeric(as.character(FAM[,col_GQ_6]) )>=50))
      # _GT and _DP filters!
      
      
    } else if (length(inds)==1){
      col_GT_1 <- first_GT
      col_DP_1 <- first_DP
      col_GQ_1 <- first_GQ
      
      sel <- which(((FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0")) & (as.numeric(as.character(FAM[,col_DP_1]))>=10) & (as.numeric(as.character(FAM[,col_GQ_1]))>=50) )
      # _GT and _DP filters!
    }
  }
  all_fams[[paste(fam)]] <- data.frame(sel) #keeping variants rows as dataframe
}
rows_from_fams <- do.call(rbind, all_fams)
ROWS_FROM_FAMS <- unique(as.vector(rows_from_fams))

EXOMES_1.0.3 <- EXOMES_1.0.2[ROWS_FROM_FAMS[,1],]





###########################################
############# COMMAS ISSUE  ###############
###########################################
ExAC_comas<-grep(",",EXOMES_1.0.3[,"dbNSFP_ExAC_Adj_AF"])
EXOMES_exac_comas<-EXOMES_1.0.3[ExAC_comas,]
EXOMES_1.0.4<-EXOMES_1.0.3[-ExAC_comas,]

gnomAD_comas<-grep(",",EXOMES_1.0.4[,"gnomAD_Ex_AF"])
EXOMES_gnomAD_comas<-EXOMES_1.0.4[gnomAD_comas,]
EXOMES_1.0.4<-EXOMES_1.0.4[-gnomAD_comas,]


cadd_comas<-grep(",",EXOMES_1.0.4$dbNSFP_CADD_phred)
EXOMES_cadd_comas<-EXOMES_1.0.4[cadd_comas,]
EXOMES_1.0.4<-EXOMES_1.0.4[-cadd_comas,]

phylop_comas<-grep(",",EXOMES_1.0.4$dbNSFP_phyloP46way_placental)
lrt_comas<-grep(",",EXOMES_1.0.4$dbNSFP_LRT_score)
pp2_comas<-grep(",",EXOMES_1.0.4$dbNSFP_Polyphen2_HVAR_pred)
mt_comas<-grep(",",EXOMES_1.0.4$dbNSFP_MutationTaster_pred)
sift_comas<-grep(",",EXOMES_1.0.4$dbNSFP_SIFT_pred)



EXOMES_1.0_comas<-rbind(EXOMES_exac_comas,EXOMES_gnomAD_comas,EXOMES_cadd_comas)

write.table(EXOMES_1.0_comas,"EXOMES_1.0_comas.txt",sep="\t",quote=F,row.names=F)

# ExAC_comas<-grep(",",EXOMES_1.0[,"dbNSFP_ExAC_Adj_AF"])
# EXOMES_1.0_exac_comas<-EXOMES_1.0[ExAC_comas,]
# 
# comas_1<-strsplit(EXOMES_1.0_exac_comas$dbNSFP_ExAC_Adj_AF,",")
# comas_1_duplicated<-sapply(sapply(sapply(comas_1,duplicated),which),length)
# comas_1_total<-sapply(comas_1,length)
# comas_iguales<-which(comas_1_duplicated+1==comas_1_total)
# 
# EXOMES_1.0_exac_comas_iguales<-EXOMES_1.0_exac_comas[comas_iguales,]

###########################################
############# COMMAS ISSUE  ###############
###########################################



#Germline allelic frequency filter

EXOMES_1.0.4[, "dbNSFP_ExAC_Adj_AF"] <- as.numeric(as.character(EXOMES_1.0.4[,"dbNSFP_ExAC_Adj_AF"]))
EXOMES_1.0.4[, "gnomAD_Ex_AF"] <- as.numeric(as.character(EXOMES_1.0.4[,"gnomAD_Ex_AF"]))


fq_filters <- which(
                (EXOMES_1.0.4$dbNSFP_ExAC_Adj_AF <= 0.05 | is.na(EXOMES_1.0.4$dbNSFP_ExAC_Adj_AF)) &
                (EXOMES_1.0.4$gnomAD_Ex_AF <= 0.05 | is.na(EXOMES_1.0.4$gnomAD_Ex_AF))
              )


data <- EXOMES_1.0.4[fq_filters,] 


### Pathogenicity prediction Tools Filtering

data$dbNSFP_phyloP46way_placental <- as.numeric(as.character(data$dbNSFP_phyloP46way_placental))
data$dbNSFP_SIFT_pred <- as.character(data$dbNSFP_SIFT_pred)
data$dbNSFP_Polyphen2_HVAR_pred <- as.character(data$dbNSFP_Polyphen2_HVAR_pred)
data$dbNSFP_CADD_phred <- as.numeric(as.character(data$dbNSFP_CADD_phred))
data$dbNSFP_LRT_pred <- as.character(data$dbNSFP_LRT_pred)
data$dbNSFP_MutationTaster_pred <- as.character(data$dbNSFP_MutationTaster_pred)

a<-as.numeric((data$dbNSFP_phyloP46way_placental)>1.6)
b<-as.numeric((data$dbNSFP_SIFT_pred)=="D")
c<-as.numeric((data$dbNSFP_Polyphen2_HVAR_pred)=="D" | (data$dbNSFP_Polyphen2_HVAR_pred)=="P")
d<-as.numeric((data$dbNSFP_CADD_phred)>15.0)
e<-as.numeric((data$dbNSFP_LRT_pred)=="D")
g<-as.numeric((data$dbNSFP_MutationTaster_pred)=="D" | (data$dbNSFP_MutationTaster_pred)=="A")

a[which(is.na(a)==TRUE)]<-0.25
b[which(data$dbNSFP_SIFT_pred=="" | data$dbNSFP_SIFT_pred==".")]<-0.25
c[which(data$dbNSFP_Polyphen2_HVAR_pred=="" | data$dbNSFP_Polyphen2_HVAR_pred==".")]<-0.25
d[which(is.na(d)==TRUE)]<-0.25
e[which(data$dbNSFP_LRT_pred=="" | data$dbNSFP_LRT_pred==".")]<-0.25
g[which(data$dbNSFP_MutationTaster_pred=="" | data$dbNSFP_MutationTaster_pred==".")]<-0.25

columns <- list(a, b, c, d, e, g)

data$extra<-apply(data.frame(columns),1,sum)

first_exome_column <- grep("_GT", colnames(data))[1]
data <- data[,c(1:(first_exome_column-1), dim(data)[2], first_exome_column:(dim(data)[2]-1))]#column position change

HIGH_rows <- which(data$Annotation_Impact == "HIGH")

if (length(HIGH_rows)==0){
  df_MISSENSE <- data
} else {
  df_MISSENSE <- data[-HIGH_rows,]
}

sel_extra <- which(df_MISSENSE$extra >= 2)

data <- rbind(df_MISSENSE[sel_extra,], data[HIGH_rows,])


################ pipeline_Maria.R ###################
source("./docs_+_pipeline_MARIA/pipeline_Marcos_7.0.R")

for (i in 1:nrow(data)) {
	data$GeneRIF[i]<-strtrim(data$GeneRIF[i],(30000))[[1]]
	data$Filt.GeneRIF[i]<-strtrim(data$Filt.GeneRIF[i],(30000))[[1]]
}
######################################################





################### gene_families_&_pathways ###################
data3 <- data

###	Gene Families
fam_gf_path<-read.table("./docs_+_pipeline_MARIA/gene_families_&_pathways/Gene_Families_con_genes_hereditarios.txt", header=T, sep="\t")

fam_functfam<-read.table("./docs_+_pipeline_MARIA/gene_families_&_pathways/Functional_Families_con_genes_hereditarios.txt", header=T, sep="\t")

prot_fam<-list()
fun_fam<-list()
par_fam<-list()

for (i in 1:dim(fam_gf_path)[1]){
	prot_fam[[as.character(fam_gf_path$Gene[i])]]<-strsplit(as.character(fam_gf_path$Protein.Family[i]),",")[[1]]
	names(prot_fam[[as.character(fam_gf_path$Gene[i])]])<-rep(as.character(fam_gf_path$Gene[i]),length(prot_fam[[as.character(fam_gf_path$Gene[i])]]))

	fun_fam[[as.character(fam_functfam$Gene[i])]]<-strsplit(as.character(fam_functfam$Functional.Family[i]),",")[[1]]
	names(fun_fam[[as.character(fam_functfam$Gene[i])]])<-rep(as.character(fam_functfam$Gene[i]),length(fun_fam[[as.character(fam_functfam$Gene[i])]]))

	par_fam[[as.character(fam_gf_path$Gene[i])]]<-strsplit(as.character(fam_gf_path$Paralogues[i]),",")[[1]]
	names(par_fam[[as.character(fam_gf_path$Gene[i])]])<-rep(as.character(fam_gf_path$Gene[i]),length(par_fam[[as.character(fam_gf_path$Gene[i])]]))
}

list_prot_fam<-do.call("c",prot_fam)
names(list_prot_fam)<-do.call(cbind,strsplit(names(list_prot_fam),"[.]"))[1,]
list_fun_fam<-do.call("c",fun_fam)
names(list_fun_fam)<-do.call(cbind,strsplit(names(list_fun_fam),"[.]"))[1,]
list_par_fam<-do.call("c",par_fam)
names(list_par_fam)<-do.call(cbind,strsplit(names(list_par_fam),"[.]"))[1,]

#################################################
genes_gf<-list_prot_fam
	
filas<-integer(0)
filas_fam<-integer(0)
for (g in 1:length(genes_gf)){
	filas<-c(filas,which(data3$unique_name==genes_gf[g]))
	filas_fam<-c(filas_fam,rep(names(genes_gf)[g],length(which(data3$unique_name==genes_gf[g]))))
}

data3$Protein.Family<-NA
data3$Protein.Family[filas]<-filas_fam
#################################################
#################################################
genes_gf<-list_fun_fam
	
filas<-integer(0)
filas_fam<-integer(0)
for (g in 1:length(genes_gf)){
	filas<-c(filas,which(data3$unique_name==genes_gf[g]))
	filas_fam<-c(filas_fam,rep(names(genes_gf)[g],length(which(data3$unique_name==genes_gf[g]))))
}

data3$Functional.Family<-NA
data3$Functional.Family[filas]<-filas_fam
#################################################
#################################################	
genes_gf<-list_par_fam
	
filas<-integer(0)
filas_fam<-integer(0)
for (g in 1:length(genes_gf)){
	filas<-c(filas,which(data3$unique_name==genes_gf[g]))
	filas_fam<-c(filas_fam,rep(names(genes_gf)[g],length(which(data3$unique_name==genes_gf[g]))))
}

data3$Paralogues<-NA
data3$Paralogues[filas]<-filas_fam
#################################################	


####	Pathways
pt<-read.table("./docs_+_pipeline_MARIA/gene_families_&_pathways/PATHWAYS_db.txt", header=T, sep="\t")

go_pt<-list()
kegg_pt<-list()
reactome_pt<-list()

for (i in 1:dim(pt)[1]){
	go_pt[[as.character(pt$Pathway[i])]]<-strsplit(as.character(pt$GO_terms[i]),",")[[1]]
	names(go_pt[[as.character(pt$Pathway[i])]])<-rep(as.character(pt$Pathway[i]),length(go_pt[[as.character(pt$Pathway[i])]]))

	kegg_pt[[as.character(pt$Pathway[i])]]<-strsplit(as.character(pt$KEGG_terms[i]),",")[[1]]
	names(kegg_pt[[as.character(pt$Pathway[i])]])<-rep(as.character(pt$Pathway[i]),length(kegg_pt[[as.character(pt$Pathway[i])]]))

	reactome_pt[[as.character(pt$Pathway[i])]]<-strsplit(as.character(pt$Reactome_terms[i]),",")[[1]]
	names(reactome_pt[[as.character(pt$Pathway[i])]])<-rep(as.character(pt$Pathway[i]),length(reactome_pt[[as.character(pt$Pathway[i])]]))
}

go_pt<-go_pt[which(is.na(go_pt)==FALSE)]
kegg_pt<-kegg_pt[which(is.na(kegg_pt)==FALSE)]
reactome_pt<-reactome_pt[which(is.na(reactome_pt)==FALSE)]

cols_pt<-as.data.frame(matrix(NA,dim(data3)[1],3))
colnames(cols_pt)<-c("GO","KEGG","Reactome")

#################################################	
pathways_go<-go_pt

filas<-list()
for (j in 1:length(pathways_go)){
	for (k in 1:length(pathways_go[[j]])){
		filas[[names(pathways_go)[j]]]<-c(filas[[names(pathways_go)[j]]],grep(pathways_go[[j]][k],data3$GeneOntology_BP,ignore.case=TRUE))
	}
	filas[[j]]<-filas[[j]][which(duplicated(filas[[j]])==FALSE)]
	
	cols_pt$GO[filas[[j]]]<-paste(cols_pt$GO[filas[[j]]],names(filas)[j],sep=", ")
}
#################################################	
#################################################	
pathways_kegg<-kegg_pt

filas<-list()
for (j in 1:length(pathways_kegg)){
	for (k in 1:length(pathways_kegg[[j]])){
		filas[[names(pathways_kegg)[j]]]<-c(filas[[names(pathways_kegg)[j]]],grep(pathways_kegg[[j]][k],data3$KEGG,ignore.case=TRUE))
	}
	filas[[j]]<-filas[[j]][which(duplicated(filas[[j]])==FALSE)]
	
	cols_pt$KEGG[filas[[j]]]<-paste(cols_pt$KEGG[filas[[j]]],names(filas)[j],sep=", ")
}
#################################################	
#################################################	
pathways_reactome<-reactome_pt

filas<-list()
for (j in 1:length(pathways_reactome)){
	for (k in 1:length(pathways_reactome[[j]])){
		filas[[names(pathways_reactome)[j]]]<-c(filas[[names(pathways_reactome)[j]]],grep(pathways_reactome[[j]][k],data3$Reactome,ignore.case=TRUE))
	}
	filas[[j]]<-filas[[j]][which(duplicated(filas[[j]])==FALSE)]
	
	cols_pt$Reactome[filas[[j]]]<-paste(cols_pt$Reactome[filas[[j]]],names(filas)[j],sep=", ")
}
#################################################	
for (w in 1:dim(cols_pt)[2]){
	cols_pt[,w]<-gsub("NA, ","",cols_pt[,w])
}

cols_pt_GO<-strsplit(cols_pt[,1],", ")
cols_pt_KEGG<-strsplit(cols_pt[,2],", ")
cols_pt_Reactome<-strsplit(cols_pt[,3],", ")
cols_def<-list()
for (a in 1:dim(data3)[1]){
	cols_def[[a]]<-c(cols_pt_GO[[a]],cols_pt_KEGG[[a]],cols_pt_Reactome[[a]])
	cols_def[[a]]<-cols_def[[a]][which(duplicated(cols_def[[a]])==FALSE)]
	cols_def[[a]]<-cols_def[[a]][which(is.na(cols_def[[a]])==FALSE)]
	cols_def[[a]]<-paste(cols_def[[a]],collapse=", ")
}
data3$Pathway<-do.call("c",cols_def)
################### gene_families_&_pathways ###################


################### 		pLI		 ###################
pli_tab<-fread("./docs_+_pipeline_MARIA/pLI/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt", sep="\t", header=TRUE, dec=".", data.table=FALSE)

data3$pLI<-NA
for (rt in 1:dim(data3)[1]){
	pli_row<-which(pli_tab$gene==data3$Gene_Name[rt])
	if (length(pli_row)>0){	
		data3$pLI[rt]<-pli_tab$pLI[pli_row]
	}
}
################### 		pLi		 ###################


# 
# 	   ##########################################################
# 	################### 	LD known regions	 ###################
# 	   ##########################################################
# 
# ld_tab<-fread("./docs_+_pipeline_MARIA/LD_regions/LIST_Regions.txt", sep="\t", header=TRUE, data.table=FALSE)
# 
# data3$LD_Region<-NA
# data3$LD_Region_bibl<-NA
# 
# for (ld in 1:dim(data3)[1]){
# 
# 	chrom_ld<-data3$CHROM[ld]
# 	pos_ld<-data3$POS[ld]
# 
# 	ld_row<-which(ld_tab$CROM==chrom_ld & ld_tab$INICIO<=pos_ld & ld_tab$FIN>=pos_ld)
# 	if (length(ld_row)>0){	
# 		data3$LD_Region[ld]<-paste(ld_tab$REGION[ld_row],collapse=" ; ")
# 		data3$LD_Region_bibl[ld]<-paste(ld_tab$AUTORES[ld_row],collapse=" ; ")
# 	}
# 	
# }
# 
# 
# 	   ##########################################################
# 	################### 	LD known regions	 ###################
# 	   ##########################################################


########################################################################################
##########|                 [aux] Cols informatives de segregaci?        |##############
######################################## START #########################################

#Programa para poner nuevas columnas con c?digo individuos que tienen la variante y familias a las q pertenecen


DATA <- data3

#Inicializaci?n nuevas columnas de DATA
DATA$Ind_HET <- NA
DATA$Ind_fam_HET <- NA
DATA$correct_SEGGREGATION_HET <- NA

DATA$Ind_HOM <- NA
DATA$Ind_fam_HOM <- NA
DATA$correct_SEGGREGATION_HOM <- NA

#Puesta de las nuevas columnas como las 4 primeras
DATA <- DATA[,c((dim(DATA)[2]-5):dim(DATA)[2], 1:(dim(DATA)[2]-6))]


#Aqu? homozigotos va a ser un vector en el que se almacenar?n los individuos que tengan la variable en homozigosis y heterozigotos los q la tengan en heterozigosis (est? hecho directamente con el n?mero de columna, por lo que habr?a que cambiarlo manualmente si se cambia cualquier otra columna o el n?mero de individuos general) [es lo q quer?a corregir con lo de antes...]

cols_GT <- grep("_GT", colnames(DATA))

for (pp in 1:dim(DATA)[1]){

	####individus
		##heterozygosis
	cols_in_GT_with_het <- which(DATA[pp,cols_GT]!="./." & DATA[pp,cols_GT]!="0/0" & DATA[pp,cols_GT]!="1/1" & DATA[pp,cols_GT]!="2/2" & DATA[pp,cols_GT]!="3/3" & DATA[pp,cols_GT]!="4/4" & DATA[pp,cols_GT]!="5/5" & DATA[pp,cols_GT]!="6/6")
	inds_het <- gsub("_GT", "", colnames(DATA)[cols_GT[cols_in_GT_with_het]])
	if (length(inds_het)>0){
		DATA[pp,"Ind_HET"] <- paste(inds_het, collapse=", ")
	}
	
		##homozygosis
	cols_in_GT_with_hom <- which(DATA[pp,cols_GT]=="1/1" | DATA[pp,cols_GT]=="2/2" | DATA[pp,cols_GT]=="3/3" | DATA[pp,cols_GT]=="4/4" | DATA[pp,cols_GT]=="5/5" | DATA[pp,cols_GT]=="6/6")
	inds_hom <- gsub("_GT", "", colnames(DATA)[cols_GT[cols_in_GT_with_hom]])
	if (length(inds_hom)>0){
		DATA[pp,"Ind_HOM"] <- paste(inds_hom, collapse=", ")
	}
	
	####families
		##heterozygosis	
	if (length(inds_het)>0){
		list_families_het <- list()
		for (vv in 1:length(inds_het)){
			name <- inds_het[vv]
			samples_row <- grep(name, SAMPLES_2[,"SAMPLES"])
			family <- SAMPLES_2[samples_row,"FAMILIES"]
			list_families_het[[vv]] <- paste(family)
		}
		FAMILIES_het <- unique(as.vector(do.call(cbind, list_families_het)))
		#if (length(FAMILIES_het)>0){
		#	DATA[pp,"Ind_fam_HET"] <- paste(FAMILIES_het, collapse=",")
		#}

			#SEGGREGATION
		inds_seggregating <- list()
		families_seggregating <- list()
		for (tt in 1:length(FAMILIES_het)){
			family <- FAMILIES_het[tt]
			family_members <- SAMPLES_2[which(SAMPLES_2[,"FAMILIES"]==family),"SAMPLES"]
			members_found <- list()
			for (rr in 1:length(family_members)){
				if (length(grep(family_members[rr], DATA[pp,"Ind_HET"]))==1){
					members_found[[paste(family_members[rr])]] <- "YES"
				} else {
					members_found[[paste(family_members[rr])]] <- "NO"
				}
			}
			comprovation <- as.vector(do.call(cbind, members_found))
			n_yes <- length(grep("YES", comprovation))
			n_no <- length(grep("NO", comprovation))
			n_members <- length(family_members)
			if (n_yes > 0){
				inds_seggregating[[paste(family)]] <- paste(family,"-",n_yes, "/", n_members, sep="")
			}
			if (n_no==0){
				families_seggregating[[tt]] <- paste(family)
			}
		}
		INDS_with_correct_SEGG <- as.vector(do.call(cbind, inds_seggregating))
		FAMS_with_correct_SEGG <- as.vector(do.call(cbind, families_seggregating))
		if (length(FAMS_with_correct_SEGG)>0){
			DATA[pp,"correct_SEGGREGATION_HET"] <- paste(unique(FAMS_with_correct_SEGG), collapse=",")
		}
		if (length(INDS_with_correct_SEGG)>0) {
			DATA[pp,"Ind_fam_HET"] <- paste(unique(INDS_with_correct_SEGG), collapse=" ; ")
		}
	}
	
		##homozygosis	
	if (length(inds_hom)>0){
		list_families_hom <- list()
		for (vv in 1:length(inds_hom)){
			name <- inds_hom[vv]
			samples_row <- grep(name, SAMPLES_2[,"SAMPLES"])
			family <- SAMPLES_2[samples_row,"FAMILIES"]
			list_families_hom[[vv]] <- paste(family)
		}
		FAMILIES_hom <- unique(as.vector(do.call(cbind, list_families_hom)))
		#if (length(FAMILIES_hom)>0){
		#	DATA[pp,"Ind_fam_HOM"] <- paste(FAMILIES_hom, collapse=",")
		#}

			#SEGGREGATION
		inds_seggregating <- list()
		families_seggregating <- list()
		for (tt in 1:length(FAMILIES_hom)){
			family <- FAMILIES_hom[tt]
			family_members <- SAMPLES_2[which(SAMPLES_2[,"FAMILIES"]==family),"SAMPLES"]
			members_found <- list()
			for (rr in 1:length(family_members)){
				if (length(grep(family_members[rr], DATA[pp,"Ind_HOM"]))==1){
					members_found[[paste(family_members[rr])]] <- "YES"
				} else {
					members_found[[paste(family_members[rr])]] <- "NO"
				}
			}
			comprovation <- as.vector(do.call(cbind, members_found))
			n_yes <- length(grep("YES", comprovation))
			n_no <- length(grep("NO", comprovation))
			n_members <- length(family_members)
			if (n_yes > 0){
				inds_seggregating[[paste(family)]] <- paste(family,"-",n_yes, "/", n_members, sep="")
			}
			if (n_no==0){
				families_seggregating[[tt]] <- paste(family)
			}
		}
		INDS_with_correct_SEGG <- as.vector(do.call(cbind, inds_seggregating))
		FAMS_with_correct_SEGG <- as.vector(do.call(cbind, families_seggregating))
		if (length(FAMS_with_correct_SEGG)>0){
			DATA[pp,"correct_SEGGREGATION_HOM"] <- paste(unique(FAMS_with_correct_SEGG), collapse=" ; ")
		}
		if (length(INDS_with_correct_SEGG)>0) {
			DATA[pp,"Ind_fam_HOM"] <- paste(unique(INDS_with_correct_SEGG), collapse=",")
		}
	}

	print(paste("Line ", pp, " from ", nrow(DATA), ". Third part.", sep=""))
}

EXOMES_3.0 <- DATA
########################################################################################
##########|                 [aux] Cols informatives de segregaci?        |##############
######################################## END ###########################################

########################################################################################
##########|  selecci? variants de gens hereditaris i hits de GWAS  |##############
######################################### START ########################################

HER_GENES <- read.table("./docs_+_pipeline_MARIA/genes_candidatos/Candidate_Genes_GC_21.08.2017.txt", sep="\t", header=F)
GWAS_GENES <- read.table("./docs_+_pipeline_MARIA/genes_candidatos/GWAS_GC_21.08.2017.txt", sep="\t", header=F)
CPG <- read.table("./docs_+_pipeline_MARIA/genes_candidatos/Cancer_predisposing_genes_Nature_2014.txt", sep="\t", header=F)#gene list of Cancer Predisposing Genes by Nature 2014 (review)

###
variants_with_hits_HER <- list()
for (fq in 1:dim(EXOMES_3.0)[1]){
	exomes_genename <- EXOMES_3.0[fq,"unique_name"]
	matching <- which(HER_GENES==exomes_genename)
	if (length(matching)==1){
		variants_with_hits_HER[[fq]] <- fq
	}
}
rows_matched_HER <- do.call(cbind, variants_with_hits_HER)
ROWS_MATCHED_HER <- as.vector(rows_matched_HER)
###
###
variants_with_hits_GWAS <- list()
for (fq in 1:dim(EXOMES_3.0)[1]){
  exomes_genename <- EXOMES_3.0[fq,"unique_name"]
  matching <- which(GWAS_GENES==exomes_genename)
  if (length(matching)==1){
    variants_with_hits_GWAS[[fq]] <- fq
  }
}
rows_matched_GWAS <- do.call(cbind, variants_with_hits_GWAS)
ROWS_MATCHED_GWAS <- as.vector(rows_matched_GWAS)
###
variants_with_hits_CPG <- list()
for (fq in 1:dim(EXOMES_3.0)[1]){
  exomes_genename <- EXOMES_3.0[fq,"unique_name"]
  matching <- which(CPG==exomes_genename)
  if (length(matching)==1){
    variants_with_hits_CPG[[fq]] <- fq
  }
}
rows_matched_CPG <- do.call(cbind, variants_with_hits_CPG)
ROWS_MATCHED_CPG <- as.vector(rows_matched_CPG)


EXOMES_3.0$Hereditary_Gene<-"NO"
EXOMES_3.0$GWAS_Gene<-"NO"
EXOMES_3.0$Cancer_Predisposition_Gene<-"NO"

EXOMES_3.0$Hereditary_Gene[ROWS_MATCHED_HER]<-"YES"
EXOMES_3.0$GWAS_Gene[ROWS_MATCHED_GWAS]<-"YES"
EXOMES_3.0$Cancer_Predisposition_Gene[ROWS_MATCHED_CPG]<-"YES"

########################################################################################
##########|  selecci? variants de gens hereditaris i hits de GWAS  |##############
######################################## END ###########################################



NAME2 <- paste(dir_output, NAME, "__WES_3.0.txt", sep="") 
write.table(EXOMES_3.0, NAME2, quote=F, row.names=F, col.names=T, sep="\t")




#loop to select family by family shared variants 
first_exome_column <- as.numeric(grep("_GT", colnames(EXOMES_3.0))[1])
colnames_tail_exomes <- as.numeric(length(as.numeric(grep("_PL", colnames(EXOMES_3.0)))))
last_exome_column <- as.numeric(grep("_PL", colnames(EXOMES_3.0))[colnames_tail_exomes])
NONE_EXOME_COLUMNS <- 1:(first_exome_column-1)
NONE_EXOME_COLUMNS_TAIL <- (last_exome_column+1):length(colnames(EXOMES_3.0))

n_cols_x_ind <- as.numeric(grep("_GT", colnames(EXOMES_3.0))[2])-first_exome_column


all_fams <- list()
for (ss in 1:length(unique(SAMPLES_2[,"FAMILIES"]))){
  fam <- as.character(unique(SAMPLES_2[,"FAMILIES"])[ss])
  inds <- as.character(SAMPLES_2[which(SAMPLES_2[,"FAMILIES"]==fam),1])
  if (length(inds)>6){
    print("Nom?s est? preparat per fam?les amb un m?xim de 6 persones!!!!")
  } else if (length(inds)<=6){
    all_cols <- list()
    for (xx in 1:length(inds)){
      name <- inds[xx]
      cols <- grep(name, colnames(EXOMES_3.0))
      all_cols[[as.character(name)]] <- cols
    }
    all_cols2 <- do.call(cbind, all_cols)
    INDS_COLS <- as.vector(all_cols2)
    
    FAM <- EXOMES_3.0[,c(NONE_EXOME_COLUMNS, INDS_COLS, NONE_EXOME_COLUMNS_TAIL)]
    first_exome_column <- as.numeric(grep("_GT", colnames(FAM))[1])
    nonexome_columns <- 1:(first_exome_column-1)
    
    first_GT <- first_exome_column
    first_DP <- as.numeric(grep("_DP", colnames(FAM))[1])
    
    if (length(inds)==2){
      col_GT_1 <- first_GT
      col_GT_2 <- first_GT + n_cols_x_ind
      
      sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0"))
      
    } else if (length(inds)==3){
      col_GT_1 <- first_GT
      col_GT_2 <- first_GT + n_cols_x_ind
      col_GT_3 <- first_GT + (n_cols_x_ind*2)
      
      sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0") & (FAM[,col_GT_3]!="./." & FAM[,col_GT_3]!="0/0"))
      
    } else if (length(inds)==4){
      col_GT_1 <- first_GT
      col_GT_2 <- first_GT + n_cols_x_ind
      col_GT_3 <- first_GT + (n_cols_x_ind*2)
      col_GT_4 <- first_GT + (n_cols_x_ind*3)
      
      sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0") & (FAM[,col_GT_3]!="./." & FAM[,col_GT_3]!="0/0") & (FAM[,col_GT_4]!="./." & FAM[,col_GT_4]!="0/0"))
      
    } else if (length(inds)==5){
      col_GT_1 <- first_GT
      col_GT_2 <- first_GT + n_cols_x_ind
      col_GT_3 <- first_GT + (n_cols_x_ind*2)
      col_GT_4 <- first_GT + (n_cols_x_ind*3)
      col_GT_5 <- first_GT + (n_cols_x_ind*4)
      
      sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0") & (FAM[,col_GT_3]!="./." & FAM[,col_GT_3]!="0/0") & (FAM[,col_GT_4]!="./." & FAM[,col_GT_4]!="0/0") & (FAM[,col_GT_5]!="./." & FAM[,col_GT_5]!="0/0"))
      
    } else if (length(inds)==6){
      col_GT_1 <- first_GT
      col_GT_2 <- first_GT + n_cols_x_ind
      col_GT_3 <- first_GT + (n_cols_x_ind*2)
      col_GT_4 <- first_GT + (n_cols_x_ind*3)
      col_GT_5 <- first_GT + (n_cols_x_ind*4)
      col_GT_6 <- first_GT + (n_cols_x_ind*5)
      
      sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0") & (FAM[,col_GT_3]!="./." & FAM[,col_GT_3]!="0/0") & (FAM[,col_GT_4]!="./." & FAM[,col_GT_4]!="0/0") & (FAM[,col_GT_5]!="./." & FAM[,col_GT_5]!="0/0") & (FAM[,col_GT_6]!="./." & FAM[,col_GT_6]!="0/0"))
      
    } else if (length(inds)==1){
      col_GT_1 <- first_GT
      
      sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0"))
      
    }
    all_fams[[fam]] <- FAM[sel,] #keeping variants rows as dataframe
  }
}

NAME2 <- paste(dir_output, NAME, "__WES_3.0_by_family.xls", sep="") 
WriteXLS("all_fams", NAME2)




########################################################################################
##########|                      [2.0 & 3.0] primers filtres                  |###############
######################################## END ###########################################






                               #######################
#-------------------#################################################------------------#
                               #######################





########################################################################################
##########|           [4.0] Termes bibliogr?fics -termes_filtres.xls-      |############
######################################### START ########################################


cols <- c("Filt.Function", "Filt.GeneOntology_BP", "Filt.GeneRIF", "Filt.KEGG", "Filt.Reactome")

all_info <- list()
for(rr in 1:dim(EXOMES_3.0)[1]){
	no_info <- length(which(is.na(EXOMES_3.0[rr,cols]) | EXOMES_3.0[rr,cols]==""))
	if (no_info < 5){
		all_info[[rr]] <- as.numeric(rr)
	}
}
all_info2 <- do.call(cbind, all_info)
ALL_INFO <- as.vector(all_info2)

EXOMES_5.0 <- EXOMES_3.0[ALL_INFO,]

########################################################################################
##########|         [4.0] Termes bibliogr?fics -termes_filtres.xls-      |##############
######################################## END ###########################################




########################################################################################
##########|                 [6.0] Segona tanda de filtres	        |###############
######################################### START ########################################
first_exome_column <- as.numeric(grep("_GT", colnames(EXOMES_5.0))[1])
nonexome_columns <- 1:(first_exome_column-1)

fq_filters <- which(
  (EXOMES_5.0$dbNSFP_ExAC_Adj_AF <= 0.01 | is.na(EXOMES_5.0$dbNSFP_ExAC_Adj_AF)) &
    (EXOMES_5.0$gnomAD_Ex_AF <= 0.01 | is.na(EXOMES_5.0$gnomAD_Ex_AF))
)

EXOMES_5.1 <- EXOMES_5.0[fq_filters,] 

HIGH_rows_2 <- which(EXOMES_5.1$Annotation_Impact == "HIGH")

if (length(HIGH_rows_2)==0){
	df_MISSENSE_2 <- EXOMES_5.1
} else {
	df_MISSENSE_2 <- EXOMES_5.1[-HIGH_rows_2,]
}

sel_6.0 <- which(df_MISSENSE_2$extra >= 3)

EXOMES_6.0 <- rbind(df_MISSENSE_2[sel_6.0,], EXOMES_5.1[HIGH_rows_2,])


NAME2 <- paste(dir_output, NAME, "__WES_6.0.txt", sep="")
#NAME3 <- paste(dir_output, NAME, "__EXOMES_6.0.xlsx", sep="")
write.table(EXOMES_6.0, NAME2, quote=F, row.names=F, col.names=T, sep="\t")
#WriteXLS(EXOMES_6.0, NAME3)
########################################################################################
##########|                 [6.0] Segona tanda de filtres    	         |##############
######################################## END ###########################################






                               #######################
#-------------------#################################################------------------#
                               #######################




#loop to select family by family shared variants 
first_exome_column <- as.numeric(grep("_GT", colnames(EXOMES_6.0))[1])
colnames_tail_exomes <- as.numeric(length(as.numeric(grep("_PL", colnames(EXOMES_6.0)))))
last_exome_column <- as.numeric(grep("_PL", colnames(EXOMES_6.0))[colnames_tail_exomes])
NONE_EXOME_COLUMNS <- 1:(first_exome_column-1)
NONE_EXOME_COLUMNS_TAIL <- (last_exome_column+1):length(colnames(EXOMES_6.0))

n_cols_x_ind <- as.numeric(grep("_GT", colnames(EXOMES_6.0))[2])-first_exome_column


all_fams <- list()
for (ss in 1:length(unique(SAMPLES_2[,"FAMILIES"]))){
	fam <- as.character(unique(SAMPLES_2[,"FAMILIES"])[ss])
	inds <- as.character(SAMPLES_2[which(SAMPLES_2[,"FAMILIES"]==fam),1])
	if (length(inds)>6){
		print("Nom?s est? preparat per fam?les amb un m?xim de 6 persones!!!!")
	} else if (length(inds)<=6){
		all_cols <- list()
		for (xx in 1:length(inds)){
			name <- inds[xx]
			cols <- grep(name, colnames(EXOMES_6.0))
			all_cols[[as.character(name)]] <- cols
		}
		all_cols2 <- do.call(cbind, all_cols)
		INDS_COLS <- as.vector(all_cols2)
		
		FAM <- EXOMES_6.0[,c(NONE_EXOME_COLUMNS, INDS_COLS, NONE_EXOME_COLUMNS_TAIL)]
		first_exome_column <- as.numeric(grep("_GT", colnames(FAM))[1])
		nonexome_columns <- 1:(first_exome_column-1)

		first_GT <- first_exome_column
		first_DP <- as.numeric(grep("_DP", colnames(FAM))[1])
	
		if (length(inds)==2){
			col_GT_1 <- first_GT
			col_GT_2 <- first_GT + n_cols_x_ind
	
			sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0"))
			# _GT and _DP filters!
		
			
		} else if (length(inds)==3){
			col_GT_1 <- first_GT
			col_GT_2 <- first_GT + n_cols_x_ind
			col_GT_3 <- first_GT + (n_cols_x_ind*2)

			sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0") & (FAM[,col_GT_3]!="./." & FAM[,col_GT_3]!="0/0"))
			# _GT and _DP filters!


		} else if (length(inds)==4){
			col_GT_1 <- first_GT
			col_GT_2 <- first_GT + n_cols_x_ind
			col_GT_3 <- first_GT + (n_cols_x_ind*2)
			col_GT_4 <- first_GT + (n_cols_x_ind*3)

			sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0") & (FAM[,col_GT_3]!="./." & FAM[,col_GT_3]!="0/0") & (FAM[,col_GT_4]!="./." & FAM[,col_GT_4]!="0/0"))
			# _GT and _DP filters!

		
		} else if (length(inds)==5){
			col_GT_1 <- first_GT
			col_GT_2 <- first_GT + n_cols_x_ind
			col_GT_3 <- first_GT + (n_cols_x_ind*2)
			col_GT_4 <- first_GT + (n_cols_x_ind*3)
			col_GT_5 <- first_GT + (n_cols_x_ind*4)

			sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0") & (FAM[,col_GT_3]!="./." & FAM[,col_GT_3]!="0/0") & (FAM[,col_GT_4]!="./." & FAM[,col_GT_4]!="0/0") & (FAM[,col_GT_5]!="./." & FAM[,col_GT_5]!="0/0"))
			# _GT and _DP filters!


		} else if (length(inds)==6){
			col_GT_1 <- first_GT
			col_GT_2 <- first_GT + n_cols_x_ind
			col_GT_3 <- first_GT + (n_cols_x_ind*2)
			col_GT_4 <- first_GT + (n_cols_x_ind*3)
			col_GT_5 <- first_GT + (n_cols_x_ind*4)
			col_GT_6 <- first_GT + (n_cols_x_ind*5)

			sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0") & (FAM[,col_GT_3]!="./." & FAM[,col_GT_3]!="0/0") & (FAM[,col_GT_4]!="./." & FAM[,col_GT_4]!="0/0") & (FAM[,col_GT_5]!="./." & FAM[,col_GT_5]!="0/0") & (FAM[,col_GT_6]!="./." & FAM[,col_GT_6]!="0/0"))
			# _GT and _DP filters!


		} else if (length(inds)==1){
			col_GT_1 <- first_GT
	
			sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0"))
			# _GT and _DP filters!
		}
		all_fams[[fam]] <- FAM[sel,] #keeping variants rows as dataframe
	}
}

NAME2 <- paste(dir_output, NAME, "__WES_6.0_by_family.xls", sep="") 
WriteXLS("all_fams", NAME2)


# 
# 
# ############################################################
# #loop to make an excel by genes
# 
# DATA <- EXOMES_6.0
# DATA$unique_name<-factor(DATA$unique_name)
# DATA$Annotation_Impact<-factor(DATA$Annotation_Impact)
# 
# Gene<-levels(DATA$unique_name)
# levels(DATA$Annotation_Impact)<-c("HIGH","MISSENSE")
# 
# df<-data.frame(Gene)
# 
# for (i in 1:length(df$Gene)){
# 	gen<-df$Gene[i]
# 	aux<-which(DATA$unique_name==gen)
# 	
# 	df[i,"Num_Vars"]<-length(aux)
# 	df[i,"Ind_fam_HET"]<-paste(DATA$Ind_fam_HET[aux],collapse=" and ")
# 	df[i,"Ind_fam_HOM"]<-paste(DATA$Ind_fam_HOM[aux],collapse=" and ")
# 	df[i,"Var_Type"]<-paste(DATA$Annotation_Impact[aux],collapse=" and ")
# 	df[i,"HGVS.c"]<-paste(DATA$HGVS.c[aux],collapse=" and ")
# 	df[i,"HGVS.p"]<-paste(DATA$HGVS.p[aux],collapse=" and ")
# 	df[i,"ExAC_AF"]<-paste(DATA$ExAC_AF[aux],collapse=" and ")
# 	df[i,"pLI"]<-DATA[aux[1],"pLI"]
# 	df[i,"LD_Region"]<-DATA[aux[1],"LD_Region"]
# 	df[i,"LD_Region_bibl"]<-DATA[aux[1],"LD_Region_bibl"]
# 	df[i,"Protein.Family"]<-DATA[aux[1],"Protein.Family"]
# 	df[i,"Functional.Family"]<-DATA[aux[1],"Functional.Family"]
# 	df[i,"Paralogues"]<-DATA[aux[1],"Paralogues"]
# 	df[i,"Pathway"]<-DATA[aux[1],"Pathway"]
# 	df[i,"RefSeq"]<-DATA[aux[1],"RefSeq"]
# 	df[i,"Function"]<-DATA[aux[1],"Function"]
# 	df[i,"NCBI_gene"]<-DATA[aux[1],"NCBI_gene"]
# 	df[i,"GeneOntology_BP"]<-DATA[aux[1],"GeneOntology_BP"]
# 	df[i,"GeneRIF"]<-DATA[aux[1],"GeneRIF"]
# 	df[i,"Interactions"]<-DATA[aux[1],"Interactions"]
# 	df[i,"KEGG"]<-DATA[aux[1],"KEGG"]
# 	df[i,"Reactome"]<-DATA[aux[1],"Reactome"]
# 	df[i,"Filt.Function"]<-DATA[aux[1],"Filt.Function"]
# 	df[i,"Filt.GeneOntology_BP"]<-DATA[aux[1],"Filt.GeneOntology_BP"]
# 	df[i,"Filt.GeneRIF"]<-DATA[aux[1],"Filt.GeneRIF"]
# 	df[i,"Filt.KEGG"]<-DATA[aux[1],"Filt.KEGG"]
# 	df[i,"Filt.Reactome"]<-DATA[aux[1],"Filt.Reactome"]
# }
# 
# df<-df[order(-df$Num_Vars),]
# NAME2 <- paste(dir_output, NAME, "__WES_6.0_by_genes.xls", sep="")
# WriteXLS(df,NAME2)
# 
# ############################################################
# 
# 
# 



                               #######################
#-------------------#################################################------------------#
                               #######################











#END#
dim(EXOMES_1.0)

dim(EXOMES_3.0)

dim(EXOMES_5.0)

dim(EXOMES_6.0)


























################################################################################
################################################################################
################################################################################
################################################################################
####################		Second part		########################
################################################################################
################################################################################
################################################################################
################################################################################

TABEXOMES_3.0<-EXOMES_3.0
TABEXOMES_6.0<-EXOMES_6.0

################################################################################
#EXOMES_3.0#####################################################################
################################################################################


TABEXOMES<-TABEXOMES_3.0

aux_comprobacion<-TABEXOMES[,c("CHROM","POS","REF","ALT")]
TABEXOMES2<-TABEXOMES[which(duplicated(aux_comprobacion)==FALSE),]
EXOMES_6.0<-TABEXOMES2


#loop to select family by family shared variants 
first_exome_column <- as.numeric(grep("_GT", colnames(EXOMES_6.0))[1])
colnames_tail_exomes <- as.numeric(length(as.numeric(grep("_PL", colnames(EXOMES_6.0)))))
last_exome_column <- as.numeric(grep("_PL", colnames(EXOMES_6.0))[colnames_tail_exomes])
NONE_EXOME_COLUMNS <- 1:(first_exome_column-1)
NONE_EXOME_COLUMNS_TAIL <- (last_exome_column+1):length(colnames(EXOMES_6.0))

n_cols_x_ind <- as.numeric(grep("_GT", colnames(EXOMES_6.0))[2])-first_exome_column

all_fams <- list()
for (ss in 1:length(unique(SAMPLES_2[,"FAMILIES"]))){
  fam <- as.character(unique(SAMPLES_2[,"FAMILIES"])[ss])
  inds <- as.character(SAMPLES_2[which(SAMPLES_2[,"FAMILIES"]==fam),1])
  if (length(inds)>6){
    print("Només està preparat per famíles amb un màxim de 6 persones!!!!")
  } else if (length(inds)<=6){
    all_cols <- list()
    for (xx in 1:length(inds)){
      name <- inds[xx]
      cols <- grep(name, colnames(EXOMES_6.0))
      all_cols[[as.character(name)]] <- cols
    }
    all_cols2 <- do.call(cbind, all_cols)
    INDS_COLS <- as.vector(all_cols2)
    
    FAM <- EXOMES_6.0[,c(NONE_EXOME_COLUMNS, INDS_COLS, NONE_EXOME_COLUMNS_TAIL)]
    first_exome_column <- as.numeric(grep("_GT", colnames(FAM))[1])
    nonexome_columns <- 1:(first_exome_column-1)
    
    first_GT <- first_exome_column
    first_DP <- as.numeric(grep("_DP", colnames(FAM))[1])
    
    if (length(inds)==2){
      col_GT_1 <- first_GT
      col_GT_2 <- first_GT + n_cols_x_ind
      
      sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0"))
      # _GT and _DP filters!
      
      
    } else if (length(inds)==3){
      col_GT_1 <- first_GT
      col_GT_2 <- first_GT + n_cols_x_ind
      col_GT_3 <- first_GT + (n_cols_x_ind*2)
      
      sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0") & (FAM[,col_GT_3]!="./." & FAM[,col_GT_3]!="0/0"))
      # _GT and _DP filters!
      
      
    } else if (length(inds)==4){
      col_GT_1 <- first_GT
      col_GT_2 <- first_GT + n_cols_x_ind
      col_GT_3 <- first_GT + (n_cols_x_ind*2)
      col_GT_4 <- first_GT + (n_cols_x_ind*3)
      
      sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0") & (FAM[,col_GT_3]!="./." & FAM[,col_GT_3]!="0/0") & (FAM[,col_GT_4]!="./." & FAM[,col_GT_4]!="0/0"))
      # _GT and _DP filters!
      
      
    } else if (length(inds)==5){
      col_GT_1 <- first_GT
      col_GT_2 <- first_GT + n_cols_x_ind
      col_GT_3 <- first_GT + (n_cols_x_ind*2)
      col_GT_4 <- first_GT + (n_cols_x_ind*3)
      col_GT_5 <- first_GT + (n_cols_x_ind*4)
      
      sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0") & (FAM[,col_GT_3]!="./." & FAM[,col_GT_3]!="0/0") & (FAM[,col_GT_4]!="./." & FAM[,col_GT_4]!="0/0") & (FAM[,col_GT_5]!="./." & FAM[,col_GT_5]!="0/0"))
      # _GT and _DP filters!
      
      
    } else if (length(inds)==6){
      col_GT_1 <- first_GT
      col_GT_2 <- first_GT + n_cols_x_ind
      col_GT_3 <- first_GT + (n_cols_x_ind*2)
      col_GT_4 <- first_GT + (n_cols_x_ind*3)
      col_GT_5 <- first_GT + (n_cols_x_ind*4)
      col_GT_6 <- first_GT + (n_cols_x_ind*5)
      
      sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0") & (FAM[,col_GT_3]!="./." & FAM[,col_GT_3]!="0/0") & (FAM[,col_GT_4]!="./." & FAM[,col_GT_4]!="0/0") & (FAM[,col_GT_5]!="./." & FAM[,col_GT_5]!="0/0") & (FAM[,col_GT_6]!="./." & FAM[,col_GT_6]!="0/0"))
      # _GT and _DP filters!
      
      
    } else if (length(inds)==1){
      col_GT_1 <- first_GT
      
      sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0"))
      # _GT and _DP filters!
    }
    all_fams[[fam]] <- FAM[sel,] #keeping variants rows as dataframe
  }
}


familias<-as.character(unique(SAMPLES_2$FAMILIES))

all_fams_2_homo<-list()
all_fams_2_chet<-list()
for (i in 1:length(familias)){
  
  asdf<-all_fams[[i]]
  
  familia<-familias[i]
  
  inds<-as.character(SAMPLES_2$SAMPLES[which(SAMPLES_2$FAMILIES==familia)])
  
  colsGTinds<-paste0(inds,"_GT")
  
  n_colsGTinds<-integer()
  for (j in 1:length(inds)){
    n_colsGTinds<-c(n_colsGTinds,which(colnames(asdf)==colsGTinds[j]))
  }
  
  n_gene_name<-grep("Gene_Name",colnames(asdf))
  
  
  for (k in 1:dim(asdf)[1]){
    aux<-asdf[,c(n_gene_name,n_colsGTinds)]
  }
  
  
  #Compound heterozygosis selection
  dupl_genes_names<-unique(aux$Gene_Name[which(duplicated(aux$Gene_Name))])
  
  total_dupl_genes_names<-integer()
  for (w in 1:length(dupl_genes_names)){
    total_dupl_genes_names<-c(total_dupl_genes_names,which(aux$Gene_Name==dupl_genes_names[w]))
  }
  
  n_selection_comp_het<-total_dupl_genes_names
  
  selection_comp_het<-asdf[n_selection_comp_het,]
  
  
  
  
  
  #Homozygosis selection
  auxl<-list()
  for (l in 1:length(n_colsGTinds)){
    if (length(which(asdf[,n_colsGTinds[l]]=="1/1" | asdf[,n_colsGTinds[l]]=="2/2" | asdf[,n_colsGTinds[l]]=="3/3" | asdf[,n_colsGTinds[l]]=="4/4" | asdf[,n_colsGTinds[l]]=="5/5" | asdf[,n_colsGTinds[l]]=="6/6"))>0){
      auxl[[l]]<-which(asdf[,n_colsGTinds[l]]=="1/1" | asdf[,n_colsGTinds[l]]=="2/2" | asdf[,n_colsGTinds[l]]=="3/3" | asdf[,n_colsGTinds[l]]=="4/4" | asdf[,n_colsGTinds[l]]=="5/5" | asdf[,n_colsGTinds[l]]=="6/6")
    } else {
      auxl[[l]]<-NA
    }
  }
  
  auxlc<-do.call("c",auxl)
  
  n_selection_homo<-as.integer(names(table(auxlc))[which(table(auxlc)==length(inds))])
  
  selection_homo<-asdf[n_selection_homo,]
  
  
  all_fams_2_chet[[familia]]<-selection_comp_het
  all_fams_2_homo[[familia]]<-selection_homo
  
  
  
  
  
  print(i)
  print(familia)
}


NAME2 <- paste(dir_output, NAME, "__","WES_3.0",".Homoz_by_family.xls", sep="") 
WriteXLS(all_fams_2_homo, NAME2)


NAME2 <- paste(dir_output, NAME, "__","WES_3.0",".Comp-Het_by_family.xls", sep="") 
WriteXLS(all_fams_2_chet, NAME2)








################################################################################
#EXOMES_6.0#####################################################################
################################################################################

TABEXOMES<-TABEXOMES_6.0
aux_comprobacion<-TABEXOMES[,c("CHROM","POS","REF","ALT")]
TABEXOMES2<-TABEXOMES[which(duplicated(aux_comprobacion)==FALSE),]
EXOMES_6.0<-TABEXOMES2


#loop to select family by family shared variants 
first_exome_column <- as.numeric(grep("_GT", colnames(EXOMES_6.0))[1])
colnames_tail_exomes <- as.numeric(length(as.numeric(grep("_PL", colnames(EXOMES_6.0)))))
last_exome_column <- as.numeric(grep("_PL", colnames(EXOMES_6.0))[colnames_tail_exomes])
NONE_EXOME_COLUMNS <- 1:(first_exome_column-1)
NONE_EXOME_COLUMNS_TAIL <- (last_exome_column+1):length(colnames(EXOMES_6.0))

n_cols_x_ind <- as.numeric(grep("_GT", colnames(EXOMES_6.0))[2])-first_exome_column

all_fams <- list()
for (ss in 1:length(unique(SAMPLES_2[,"FAMILIES"]))){
  fam <- as.character(unique(SAMPLES_2[,"FAMILIES"])[ss])
  inds <- as.character(SAMPLES_2[which(SAMPLES_2[,"FAMILIES"]==fam),1])
  if (length(inds)>6){
    print("Només està preparat per famíles amb un màxim de 6 persones!!!!")
  } else if (length(inds)<=6){
    all_cols <- list()
    for (xx in 1:length(inds)){
      name <- inds[xx]
      cols <- grep(name, colnames(EXOMES_6.0))
      all_cols[[as.character(name)]] <- cols
    }
    all_cols2 <- do.call(cbind, all_cols)
    INDS_COLS <- as.vector(all_cols2)
    
    FAM <- EXOMES_6.0[,c(NONE_EXOME_COLUMNS, INDS_COLS, NONE_EXOME_COLUMNS_TAIL)]
    first_exome_column <- as.numeric(grep("_GT", colnames(FAM))[1])
    nonexome_columns <- 1:(first_exome_column-1)
    
    first_GT <- first_exome_column
    first_DP <- as.numeric(grep("_DP", colnames(FAM))[1])
    
    if (length(inds)==2){
      col_GT_1 <- first_GT
      col_GT_2 <- first_GT + n_cols_x_ind
      
      sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0"))
      # _GT and _DP filters!
      
      
    } else if (length(inds)==3){
      col_GT_1 <- first_GT
      col_GT_2 <- first_GT + n_cols_x_ind
      col_GT_3 <- first_GT + (n_cols_x_ind*2)
      
      sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0") & (FAM[,col_GT_3]!="./." & FAM[,col_GT_3]!="0/0"))
      # _GT and _DP filters!
      
      
    } else if (length(inds)==4){
      col_GT_1 <- first_GT
      col_GT_2 <- first_GT + n_cols_x_ind
      col_GT_3 <- first_GT + (n_cols_x_ind*2)
      col_GT_4 <- first_GT + (n_cols_x_ind*3)
      
      sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0") & (FAM[,col_GT_3]!="./." & FAM[,col_GT_3]!="0/0") & (FAM[,col_GT_4]!="./." & FAM[,col_GT_4]!="0/0"))
      # _GT and _DP filters!
      
      
    } else if (length(inds)==5){
      col_GT_1 <- first_GT
      col_GT_2 <- first_GT + n_cols_x_ind
      col_GT_3 <- first_GT + (n_cols_x_ind*2)
      col_GT_4 <- first_GT + (n_cols_x_ind*3)
      col_GT_5 <- first_GT + (n_cols_x_ind*4)
      
      sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0") & (FAM[,col_GT_3]!="./." & FAM[,col_GT_3]!="0/0") & (FAM[,col_GT_4]!="./." & FAM[,col_GT_4]!="0/0") & (FAM[,col_GT_5]!="./." & FAM[,col_GT_5]!="0/0"))
      # _GT and _DP filters!
      
      
    } else if (length(inds)==6){
      col_GT_1 <- first_GT
      col_GT_2 <- first_GT + n_cols_x_ind
      col_GT_3 <- first_GT + (n_cols_x_ind*2)
      col_GT_4 <- first_GT + (n_cols_x_ind*3)
      col_GT_5 <- first_GT + (n_cols_x_ind*4)
      col_GT_6 <- first_GT + (n_cols_x_ind*5)
      
      sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0") & (FAM[,col_GT_2]!="./." & FAM[,col_GT_2]!="0/0") & (FAM[,col_GT_3]!="./." & FAM[,col_GT_3]!="0/0") & (FAM[,col_GT_4]!="./." & FAM[,col_GT_4]!="0/0") & (FAM[,col_GT_5]!="./." & FAM[,col_GT_5]!="0/0") & (FAM[,col_GT_6]!="./." & FAM[,col_GT_6]!="0/0"))
      # _GT and _DP filters!
      
      
    } else if (length(inds)==1){
      col_GT_1 <- first_GT
      
      sel <- which( (FAM[,col_GT_1]!="./." & FAM[,col_GT_1]!="0/0"))
      # _GT and _DP filters!
    }
    all_fams[[fam]] <- FAM[sel,] #keeping variants rows as dataframe
  }
}


familias<-as.character(unique(SAMPLES_2$FAMILIES))

all_fams_2_homo<-list()
all_fams_2_chet<-list()
for (i in 1:length(familias)){
  
  asdf<-all_fams[[i]]
  
  familia<-familias[i]
  
  inds<-as.character(SAMPLES_2$SAMPLES[which(SAMPLES_2$FAMILIES==familia)])
  
  colsGTinds<-paste0(inds,"_GT")
  
  n_colsGTinds<-integer()
  for (j in 1:length(inds)){
    n_colsGTinds<-c(n_colsGTinds,which(colnames(asdf)==colsGTinds[j]))
  }
  
  n_gene_name<-grep("Gene_Name",colnames(asdf))
  
  
  for (k in 1:dim(asdf)[1]){
    aux<-asdf[,c(n_gene_name,n_colsGTinds)]
  }
  
  
  #Compound heterozygosis selection
  dupl_genes_names<-unique(aux$Gene_Name[which(duplicated(aux$Gene_Name))])
  
  total_dupl_genes_names<-integer()
  for (w in 1:length(dupl_genes_names)){
    total_dupl_genes_names<-c(total_dupl_genes_names,which(aux$Gene_Name==dupl_genes_names[w]))
  }
  
  n_selection_comp_het<-total_dupl_genes_names
  
  selection_comp_het<-asdf[n_selection_comp_het,]
  
  
  
  
  
  #Homozygosis selection
  auxl<-list()
  for (l in 1:length(n_colsGTinds)){
    if (length(which(asdf[,n_colsGTinds[l]]=="1/1" | asdf[,n_colsGTinds[l]]=="2/2" | asdf[,n_colsGTinds[l]]=="3/3" | asdf[,n_colsGTinds[l]]=="4/4" | asdf[,n_colsGTinds[l]]=="5/5" | asdf[,n_colsGTinds[l]]=="6/6"))>0){
      auxl[[l]]<-which(asdf[,n_colsGTinds[l]]=="1/1" | asdf[,n_colsGTinds[l]]=="2/2" | asdf[,n_colsGTinds[l]]=="3/3" | asdf[,n_colsGTinds[l]]=="4/4" | asdf[,n_colsGTinds[l]]=="5/5" | asdf[,n_colsGTinds[l]]=="6/6")
    } else {
      auxl[[l]]<-NA
    }
  }
  
  auxlc<-do.call("c",auxl)
  
  n_selection_homo<-as.integer(names(table(auxlc))[which(table(auxlc)==length(inds))])
  
  selection_homo<-asdf[n_selection_homo,]
  
  
  all_fams_2_chet[[familia]]<-selection_comp_het
  all_fams_2_homo[[familia]]<-selection_homo
  
  
  
  
  
  print(i)
  print(familia)
}


NAME2 <- paste(dir_output, NAME, "__","WES_6.0",".Homoz_by_family.xls", sep="") 
WriteXLS(all_fams_2_homo, NAME2)


NAME2 <- paste(dir_output, NAME, "__","WES_6.0",".Comp-Het_by_family.xls", sep="") 
WriteXLS(all_fams_2_chet, NAME2)












