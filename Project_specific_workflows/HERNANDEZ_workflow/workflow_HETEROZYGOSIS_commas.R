#!/usr/bin/R
#	AUTHOR: Marcos Díaz-Gay
#	EMAIL: diaz2@clinic.cat
#	LAB: Genetic Predisposition to Gastrointestinal Cancer Group (Sergi Castellvi's Lab)
#	DESCRIPTION: Comma issue handling for workflow for germline SNV analysis of different samples provided in different tsv files


# Loading libraries needed      
library(data.table)


# Loading files needed      
input_files <- read.table("files_to_analyze.txt", header=T, sep="\t", stringsAsFactors=F)
direct <- read.table("directories.txt", header=T, sep="\t", stringsAsFactors=F)


# Directory creation      
dir_input <- direct$Input_dir
dir_output <- direct$Output_dir
dir_workflow <- direct$Workflow_dir
dir_annotation <- direct$Annotation_dir

## Sample loop

for (sample in 1:nrow(input_files)){


  NAME <- as.character(input_files[sample,"NAME"])

  # Loading commas issue input file
  asdf<-fread(paste0(dir_output,"/",NAME,"__WES_1.0_comas.txt"),sep="\t",header=T,data.table = F)

  filt_cols<-c("dbNSFP_ExAC_Adj_AF","gnomAD_Ex_AF","dbNSFP_phyloP46way_placental","dbNSFP_SIFT_pred","dbNSFP_Polyphen2_HVAR_pred","dbNSFP_CADD_phred","dbNSFP_LRT_pred","dbNSFP_MutationTaster_pred")


  # Germline allelic frequency filter

  #ExAC
  exac_split<-strsplit(as.character(asdf[,filt_cols[1]]),",")
  exac_split_filt_1<-grep("TRUE",sapply(lapply(exac_split,as.numeric),"<",0.01))
  exac_split_filt_2<-grep("TRUE",sapply(lapply(exac_split,as.numeric),is.na))
  exac_split_filt_3<-grep("0",sapply(lapply(exac_split,as.numeric),length))
  exac_split_filt<-sort(c(exac_split_filt_1,exac_split_filt_2,exac_split_filt_3))
  exac_split_filt<-exac_split_filt[which(duplicated(exac_split_filt)==FALSE)]

  asdf<-asdf[exac_split_filt,]

  #gnomAD
  exac_split<-strsplit(as.character(asdf[,filt_cols[2]]),",")
  exac_split_filt_1<-grep("TRUE",sapply(lapply(exac_split,as.numeric),"<",0.01))
  exac_split_filt_2<-grep("TRUE",sapply(lapply(exac_split,as.numeric),is.na))
  exac_split_filt_3<-grep("0",sapply(lapply(exac_split,as.numeric),length))
  exac_split_filt<-sort(c(exac_split_filt_1,exac_split_filt_2,exac_split_filt_3))
  exac_split_filt<-exac_split_filt[which(duplicated(exac_split_filt)==FALSE)]

  asdf<-asdf[exac_split_filt,]


  # Pathogenicity prediction Tools annotation

  ##phyloP
  phylop_split<-strsplit(as.character(asdf[,filt_cols[3]]),",")
  phylop_split_na<-grep("0",sapply(sapply(lapply(phylop_split,as.numeric),">",1.6),length))
  phylop_split_filt<-grep("TRUE",lapply(lapply(phylop_split,as.numeric),">",1.6))
  a<-rep(0,length(phylop_split))
  a[phylop_split_na]<-0.25
  a[phylop_split_filt]<-1


  #SIFT
  sift_split<-strsplit(as.character(asdf[,filt_cols[4]]),",")
  sift_split_na_1<-which(sapply(sift_split,length)==0)
  sift_split_na_2<-grep("[.]",sift_split)
  sift_split_na<-sort(c(sift_split_na_1,sift_split_na_2))
  sift_split_filt<-grep("D",sift_split)
  b<-rep(0,length(sift_split))
  b[sift_split_na]<-0.25
  b[sift_split_filt]<-1


  #Polyphen2_HVAR
  pp2_split<-strsplit(as.character(asdf[,filt_cols[5]]),",")
  pp2_split_na_1<-which(sapply(pp2_split,length)==0)
  pp2_split_na_2<-grep("[.]",pp2_split)
  pp2_split_na<-sort(c(pp2_split_na_1,pp2_split_na_2))
  pp2_split_filt_1<-grep("D",pp2_split)
  pp2_split_filt_2<-grep("P",pp2_split)
  pp2_split_filt<-sort(c(pp2_split_filt_1,pp2_split_filt_2))
  pp2_split_filt<-pp2_split_filt[which(duplicated(pp2_split_filt)==FALSE)]

  c<-rep(0,length(pp2_split))
  c[pp2_split_na]<-0.25
  c[pp2_split_filt]<-1


  #CADD
  cadd_split<-strsplit(as.character(asdf[,filt_cols[6]]),",")
  cadd_split_na<-grep("0",sapply(sapply(lapply(cadd_split,as.numeric),">",15.0),length))
  cadd_split_filt<-grep("TRUE",lapply(lapply(cadd_split,as.numeric),">",15.0))
  d<-rep(0,length(cadd_split))
  d[cadd_split_na]<-0.25
  d[cadd_split_filt]<-1


  #LRT
  lrt_split<-strsplit(as.character(asdf[,filt_cols[7]]),",")
  lrt_split_na_1<-which(sapply(lrt_split,length)==0)
  lrt_split_na_2<-grep("[.]",lrt_split)
  lrt_split_na<-sort(c(lrt_split_na_1,lrt_split_na_2))
  lrt_split_filt<-grep("D",lrt_split)
  e<-rep(0,length(lrt_split))
  e[lrt_split_na]<-0.25
  e[lrt_split_filt]<-1


  #MutationTaster
  mutationtaster_split<-strsplit(as.character(asdf[,filt_cols[8]]),",")
  mutationtaster_split_na_1<-which(sapply(mutationtaster_split,length)==0)
  mutationtaster_split_na_2<-grep("[.]",mutationtaster_split)
  mutationtaster_split_na<-sort(c(mutationtaster_split_na_1,mutationtaster_split_na_2))
  mutationtaster_split_filt_1<-grep("D",mutationtaster_split)
  mutationtaster_split_filt_2<-grep("A",mutationtaster_split)
  mutationtaster_split_filt<-sort(c(mutationtaster_split_filt_1,mutationtaster_split_filt_2))
  mutationtaster_split_filt<-mutationtaster_split_filt[which(duplicated(mutationtaster_split_filt)==FALSE)]
  g<-rep(0,length(mutationtaster_split))
  g[mutationtaster_split_na]<-0.25
  g[mutationtaster_split_filt]<-1

  ###
  columns <- list(a, b, c, d, e, g)
  data<-asdf
  data$extra<-apply(data.frame(columns),1,sum)


  # HIGH / missense with >=2 tools filtering
  HIGH_rows <- which(data$`ANN[*].IMPACT` == "HIGH")
  df_MISSENSE <- data[grep("missense", data$`ANN[*].EFFECT`),]

  sel_extra <- which(df_MISSENSE$extra >= 2)

  data <- rbind(df_MISSENSE[sel_extra,], data[HIGH_rows,])


  # Annotation pipeline
  source(paste0(dir_workflow,"/SNV_germline_WES_analysis/Annotation/pipeline_Marcos_7.0.R"))

  # Trimming GeneRIF results
  for (i in 1:nrow(data)) {
    data$GeneRIF[i]<-strtrim(data$GeneRIF[i],(30000))[[1]]
    data$Filt.GeneRIF[i]<-strtrim(data$Filt.GeneRIF[i],(30000))[[1]]
  }


    ################### gene_families_&_pathways ###################
    data3 <- data

    ###	Gene Families
    fam_gf_path<-read.table(paste0(dir_annotation,"/gene_families_&_pathways_CRC/Gene_Families_con_genes_hereditarios.txt"), header=T, sep="\t")

    fam_functfam<-read.table(paste0(dir_annotation,"/gene_families_&_pathways_CRC/Functional_Families_con_genes_hereditarios.txt"), header=T, sep="\t")

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
    pt<-read.table(paste0(dir_annotation,"/gene_families_&_pathways_CRC/PATHWAYS_db.txt"), header=T, sep="\t")

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
    pli_tab<-fread(paste0(dir_annotation,"/pLI/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt"), sep="\t", header=TRUE, dec=".", data.table=FALSE)

    data3$pLI<-NA
    for (rt in 1:dim(data3)[1]){
      pli_row<-which(pli_tab$gene==data3[rt,"ANN[*].GENE"])
      if (length(pli_row)>0){	
        data3$pLI[rt]<-pli_tab$pLI[pli_row]
      }
    }
    ################### 		pLi		 ###################



    ################### 	LD known regions	 ###################
    ld_tab<-fread(paste0(dir_annotation,"/LD_regions_CRC/LIST_Regions.txt"), sep="\t", header=TRUE, data.table=FALSE)
    
    data3$LD_Region<-NA
    data3$LD_Region_bibl<-NA

    for (ld in 1:dim(data3)[1]){
    
      chrom_ld<-data3$CHROM[ld]
      pos_ld<-data3$POS[ld]
      
      ld_row<-which(ld_tab$CROM==chrom_ld & ld_tab$INICIO<=pos_ld & ld_tab$FIN>=pos_ld)
      if (length(ld_row)>0){	
        data3$LD_Region[ld]<-paste(ld_tab$REGION[ld_row],collapse=" ; ")
        data3$LD_Region_bibl[ld]<-paste(ld_tab$AUTORES[ld_row],collapse=" ; ")
      }
    }
    ################### 	LD known regions	 ###################


    ########################################################################################
    ##########|  selecció variants de gens hereditaris i hits de GWAS  |##############
    ######################################### START ########################################

    HER_GENES <- read.table(paste0(dir_annotation,"/genes_candidatos_CRC/genes_candidatos_23.09.2016.txt"), sep="\t", header=F)
    CPG <- read.table(paste0(dir_annotation,"/genes_candidatos_CRC/Cancer_predisposing_genes_Nature_2014.txt"), sep="\t", header=F)#gene list of Cancer Predisposing Genes by Nature 2014 (review)

    ###
    EXOMES_3.0<-data3

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
    EXOMES_3.0$Cancer_Predisposition_Gene<-"NO"

    EXOMES_3.0$Hereditary_Gene[ROWS_MATCHED_HER]<-"YES"
    EXOMES_3.0$Cancer_Predisposition_Gene[ROWS_MATCHED_CPG]<-"YES"

    ########################################################################################
    ##########|  selecció variants de gens hereditaris i hits de GWAS  |##############
    ######################################## END ###########################################


  NAME2 <- paste(dir_output, NAME, "__WES_3.0_commas.txt", sep="") 
  write.table(EXOMES_3.0, NAME2, quote=F, row.names=F, col.names=T, sep="\t")




  ########################################################################################
  ##########|           [4.0] Termes bibliogràfics -termes_filtres.xls-      |############
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
  ##########|         [4.0] Termes bibliogràfics -termes_filtres.xls-      |##############
  ######################################## END ###########################################





  ########################################################################################
  ##########|                 [6.0] Segona tanda de filtres	        |###############
  ######################################### START ########################################




  ###########################################################################################################
  ### Germline allelic frequency filter

  asdf<-EXOMES_5.0

  #ExAC
  exac_split<-strsplit(asdf[,filt_cols[1]],",")
  exac_split_filt_1<-grep("TRUE",sapply(lapply(exac_split,as.numeric),"<",0.001))
  exac_split_filt_2<-grep("TRUE",sapply(lapply(exac_split,as.numeric),is.na))
  exac_split_filt_3<-grep("0",sapply(lapply(exac_split,as.numeric),length))
  exac_split_filt<-sort(c(exac_split_filt_1,exac_split_filt_2,exac_split_filt_3))
  exac_split_filt<-exac_split_filt[which(duplicated(exac_split_filt)==FALSE)]

  asdf<-asdf[exac_split_filt,]

  #gnomAD
  exac_split<-strsplit(asdf[,filt_cols[2]],",")
  exac_split_filt_1<-grep("TRUE",sapply(lapply(exac_split,as.numeric),"<",0.001))
  exac_split_filt_2<-grep("TRUE",sapply(lapply(exac_split,as.numeric),is.na))
  exac_split_filt_3<-grep("0",sapply(lapply(exac_split,as.numeric),length))
  exac_split_filt<-sort(c(exac_split_filt_1,exac_split_filt_2,exac_split_filt_3))
  exac_split_filt<-exac_split_filt[which(duplicated(exac_split_filt)==FALSE)]

  asdf<-asdf[exac_split_filt,]

  ###########################################################################################################


  EXOMES_5.1 <- asdf

  HIGH_rows_2 <- which(EXOMES_5.1$`ANN[*].IMPACT` == "HIGH")

  if (length(HIGH_rows_2)==0){
    df_MISSENSE_2 <- EXOMES_5.1
  } else {
    df_MISSENSE_2 <- EXOMES_5.1[-HIGH_rows_2,]
  }

  sel_6.0 <- which(df_MISSENSE_2$extra >= 3)

  EXOMES_6.0 <- rbind(df_MISSENSE_2[sel_6.0,], EXOMES_5.1[HIGH_rows_2,])


  NAME2 <- paste(dir_output, NAME, "__WES_6.0_commas.txt", sep="")
  write.table(EXOMES_6.0, NAME2, quote=F, row.names=F, col.names=T, sep="\t")

  ########################################################################################
  ##########|                 [6.0] Segona tanda de filtres    	         |##############
  ######################################## END ###########################################


  #END#
  dim(EXOMES_3.0)
  dim(EXOMES_5.0)
  dim(EXOMES_6.0)

}












