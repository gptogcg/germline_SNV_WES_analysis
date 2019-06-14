#Variant DB annotation pipeline (Maria Vila-Casadesus / Sebastia Franch-Exposito / Marcos Diaz-Gay)


######################################################################
	##########   DBs reading FIRST PART (start)   ##########
######################################################################

library("data.table")
library("annotate")
library("GO.db")
library("org.Hs.eg.db")
library("gdata")
library("KEGG.db")
library("reactome.db")

ruta <- paste0(dir_annotation,"/DBs_annotation/")

#HGNC_synonyms -> HGNC synonyms for previous symbols used for CNAG variant annotators
HGNC_synonyms<-fread(paste(ruta,"HGNC_synonyms.txt", sep=""), header=TRUE, sep="\t", data.table=F)
colnames(HGNC_synonyms)<-c("HGNC","Previous_Symbols","Synonyms")
HGNC_synonyms$true_synonyms<-paste0(HGNC_synonyms$Previous_Symbols,", ",HGNC_synonyms$Synonyms)



#CONVNORM -> HGNC ID vs RefSeq transcript
HGNC_RefSeq<-fread(paste(ruta, "HGNC_RefSeq.txt", sep=""),header=TRUE,sep="\t",data.table=F)
colnames(HGNC_RefSeq)<-c("HGNC","RefSeq")
HGNC_RefSeq[,"HGNC"]<-as.character(HGNC_RefSeq[,"HGNC"])
HGNC_RefSeq[,"RefSeq"]<-as.character(HGNC_RefSeq[,"RefSeq"])


#CONVENTR -> HGNC ID vs NCBI_gene ID
HGNC_entrez<-fread(paste(ruta, "HGNC_entrez.txt", sep=""),header=TRUE,sep="\t",data.table=F)
colnames(HGNC_entrez)<-c("HGNC","NCBI_gene")
HGNC_entrez[,"NCBI_gene"]<-as.character(HGNC_entrez[,"NCBI_gene"])
HGNC_entrez[,"HGNC"]<-as.character(HGNC_entrez[,"HGNC"])


#SUMMARY.REF -> RefSeq transcript ID and function description
summary.ref<-fread(paste(ruta, "refSeqSummaryfilt_noalmohad.txt", sep=""),sep="\t", header=FALSE, data.table=F)
colnames(summary.ref)<-c("RefSeq","Status","Function")
summary.ref[,"RefSeq"]<-as.character(summary.ref[,"RefSeq"])
summary.ref[,"Function"]<-as.character(summary.ref[,"Function"])


#GOID_GOTERM
GOID_GOTERM<-as.data.frame(GOTERM)
colnames(GOID_GOTERM)[2]<-"go_id2"
GOID_GOTERM[,"go_id"]<-as.character(GOID_GOTERM[,"go_id"])
GOID_GOTERM<-GOID_GOTERM[which(GOID_GOTERM$Ontology=="BP"),]


#GOID_entrez
GOID_entrez<-as.data.frame(org.Hs.egGO)
GOID_entrez[,"gene_id"]<-as.character(GOID_entrez[,"gene_id"])
GOID_entrez[,"go_id"]<-as.character(GOID_entrez[,"go_id"])
GOID_entrez<-GOID_entrez[which(GOID_entrez$Ontology=="BP"),]


#GENERIFS -> NCBI_gene ID and GeneRIF text
generifs<-fread(paste(ruta, "generifs_basic", sep=""),sep="\t",header=TRUE,,data.table=F)
colnames(generifs)<-c("Tax_ID","Gene_ID","PubMed_ID","last_update_timestamp","GeneRIF_text")
generifs[,"Gene_ID"]<-as.character(generifs[,"Gene_ID"])
generifs[,"GeneRIF_text"]<-as.character(generifs[,"GeneRIF_text"])


#GENE2GO -> NCBI_gene ID and GO term and category
gene2go<-fread(paste(ruta, "gene2go",sep=""), header=TRUE, sep="\t", data.table=F)
gene2go[,"GeneID"]<-as.character(gene2go[,"GeneID"])
gene2go[,"Category"]<-as.character(gene2go[,"Category"])
gene2go[,"GO_term"]<-as.character(gene2go[,"GO_term"])
gene2go<-gene2go[which(gene2go$Category=="Process"),]


#INTERACTIONS
interactions<-fread(paste(ruta,"interactions",sep="") ,header=TRUE,sep="\t", data.table=F)
interactions[,"tax_id"]<-as.character(interactions[,"tax_id"])
interactions[,"gene_id"]<-as.character(interactions[,"gene_id"])
interactions[,"interactant_id"]<-as.character(interactions[,"interactant_id"])


#KEGG_entrez -> KEGG ID and NCBI_gene ID
KEGG_entrez<-as.data.frame(KEGGEXTID2PATHID)
KEGG_entrez[,"pathway_id"]<-as.character(KEGG_entrez[,"pathway_id"])
KEGG_entrez[,"gene_or_orf_id"]<-as.character(KEGG_entrez[,"gene_or_orf_id"])

#KEGG_path -> KEGG ID and KEGG pathway
KEGG_path<-as.data.frame(KEGGPATHID2NAME)
KEGG_path[,"path_id"]<-paste0("hsa",as.character(KEGG_path[,"path_id"]))
KEGG_path[,"path_name"]<-as.character(KEGG_path[,"path_name"])


#REACT_entrez -> Reactome ID and NCBI_gene ID
REACT_entrez<-as.data.frame(reactomeEXTID2PATHID)
REACT_entrez[,"DB_ID"]<-as.character(REACT_entrez[,"DB_ID"])
REACT_entrez[,"gene_id"]<-as.character(REACT_entrez[,"gene_id"])


#REACT_path -> Reactome ID and Reactome ID
REACT_path<-as.data.frame(reactomePATHID2NAME)
REACT_path[,"DB_ID"]<-as.character(REACT_path[,"DB_ID"])
REACT_path<-REACT_path[grep("Homo sapiens:",as.character(REACT_path$path_name)),]
REACT_path[,"path_name"]<-as.character(REACT_path[,"path_name"])
REACT_path$path_name<-sapply(strsplit(REACT_path[,"path_name"],"Homo sapiens: "),"[",2)


#HGNC_OMIM -> HGNC ID and MIM ID
HGNC_OMIM<-fread(paste(ruta,"HGNC_OMIM.txt",sep=""), header=TRUE, sep="\t", data.table=F)
colnames(HGNC_OMIM)<-c("HGNC","OMIM_ID")
HGNC_OMIM[,"HGNC"]<-as.character(HGNC_OMIM[,"HGNC"])
HGNC_OMIM[,"OMIM_ID"]<-as.character(HGNC_OMIM[,"OMIM_ID"])


#OMIM_morbidmap -> Phenotype ID and MIM ID
OMIM_morbidmap<-fread(paste(ruta,"OMIM_morbidmap.txt",sep=""), header=TRUE, sep="\t", data.table=F)
colnames(OMIM_morbidmap)<-c("Phenotype","Gene_Symbols","MIM_Number","Cyto_Location")
OMIM_morbidmap[,"Phenotype"]<-as.character(OMIM_morbidmap[,"Phenotype"])
OMIM_morbidmap[,"MIM_Number"]<-as.character(OMIM_morbidmap[,"MIM_Number"])


#HGNC_ProteinAtlas_normal -> Protein Atlas Normal Tissue vs HGNC
HGNC_ProteinAtlas_normal<-fread(paste(ruta,"HGNC_ProteinAtlas_normal_tissue.csv",sep=""),header=TRUE, sep=",", data.table=F)
colnames(HGNC_ProteinAtlas_normal)<-c("Gene","Gene_name","Tissue","Cell_type","Level","Reliability")
HGNC_ProteinAtlas_normal<-HGNC_ProteinAtlas_normal[which(HGNC_ProteinAtlas_normal$Tissue=="colon" | HGNC_ProteinAtlas_normal$Tissue=="rectum"),]
HGNC_ProteinAtlas_normal[,"Gene_name"]<-as.character(HGNC_ProteinAtlas_normal[,"Gene_name"])
HGNC_ProteinAtlas_normal[,"Level"]<-as.character(HGNC_ProteinAtlas_normal[,"Level"])


#HGNC_ProteinAtlas_cancer -> Protein Atlas Cancer Tissue vs HGNC
HGNC_ProteinAtlas_cancer<-fread(paste(ruta,"HGNC_ProteinAtlas_cancer.csv",sep=""),header=TRUE, sep=",", data.table=F)
colnames(HGNC_ProteinAtlas_cancer)<-c("Gene","Gene_name","Tumor","Level","Count_patients","Total_patients")
HGNC_ProteinAtlas_cancer<-HGNC_ProteinAtlas_cancer[which(HGNC_ProteinAtlas_cancer$Count_patients!=0),]
HGNC_ProteinAtlas_cancer$Level_mod<-paste0(HGNC_ProteinAtlas_cancer$Level,"-",HGNC_ProteinAtlas_cancer$Count_patients)
HGNC_ProteinAtlas_cancer<-HGNC_ProteinAtlas_cancer[which(HGNC_ProteinAtlas_cancer$Tumor=="colorectal cancer"),]
HGNC_ProteinAtlas_cancer[,"Gene_name"]<-as.character(HGNC_ProteinAtlas_cancer[,"Gene_name"])
HGNC_ProteinAtlas_cancer[,"Level_mod"]<-as.character(HGNC_ProteinAtlas_cancer[,"Level_mod"])

#HGNC_ProteinAtlas_subcellular_location -> Protein Atlas Subcellular Location vs HGNC
HGNC_ProteinAtlas_sl<-fread(paste(ruta,"HGNC_ProteinAtlas_subcellular_location.csv",sep=""), header=TRUE, sep=",", data.table=F)
colnames(HGNC_ProteinAtlas_sl)<-c("Gene", "Gene_name", "Reliability", "Validated", "Supported", "Approved", "Uncertain", "Cell-to-cell_variation", "Cell-to-cell_variation_intensity","Cell_cycle_dependency","GO_id")
HGNC_ProteinAtlas_sl[,"Gene_name"]<-as.character(HGNC_ProteinAtlas_sl[,"Gene_name"])
HGNC_ProteinAtlas_sl[,"Validated"]<-as.character(HGNC_ProteinAtlas_sl[,"Validated"])
HGNC_ProteinAtlas_sl[,"Supported"]<-as.character(HGNC_ProteinAtlas_sl[,"Supported"])
HGNC_ProteinAtlas_sl[,"Approved"]<-as.character(HGNC_ProteinAtlas_sl[,"Approved"])
HGNC_ProteinAtlas_sl[,"Uncertain"]<-as.character(HGNC_ProteinAtlas_sl[,"Uncertain"])

######################################################################
	##########   DBs reading FIRST PART (end)   ##########
######################################################################

#New columns initialization
data$unique_name<-NA
data$RefSeq<-NA
data$Function<-NA
data$NCBI_gene<-NA
data$GeneOntology_BP<-NA
data$GeneRIF<-NA
data$GeneRIF_PubMedID<-NA
data$Interactions<-NA
data$KEGG<-NA
data$Reactome<-NA
data$OMIM_ID<-NA
data$OMIM_Phenotype<-NA
data$Expression_normal<-NA
data$Expression_cancer<-NA
data$Subcellular_location<-NA




#Creation UNIQUE_NAME previous data
data[,"ANN[*].GENE"]<-as.character(data[,"ANN[*].GENE"])
noms<-data[,"ANN[*].GENE"]
noms.spl<-sapply(strsplit(noms,";"),"[",1)

#############################################################
	##########   LOOP FIRST PART (start)   ##########
#############################################################


for (i in 1:nrow(data)) {

#Annotation unique_name / RefSeq /Function / NCBI_gene
	gen<-noms.spl[i]
	if (length(grep(paste0("\\b",gen,"\\b"),HGNC_synonyms$HGNC))>0){
		data$unique_name[i]<-gen
	} else {	
		selSYN<-grep(paste0("\\b",gen,"\\b"),HGNC_synonyms$true_synonyms)
		if (length(selSYN)>0){
			data$unique_name[i]<-HGNC_synonyms$HGNC[selSYN][1]
		} else {
			data$unique_name[i]<-gen
		}
	}
	sel<-which(HGNC_RefSeq[,"HGNC"]==data$unique_name[i])[1]		
	data$RefSeq[i]<-HGNC_RefSeq[sel,"RefSeq"]
	sel2<-which(summary.ref[,"RefSeq"]==data$RefSeq[i])[1]
	data$Function[i]<-summary.ref[sel2,"Function"]
	sel3<-which(HGNC_entrez[,"HGNC"]==data$unique_name[i])[1]		
	data$NCBI_gene[i]<-HGNC_entrez[sel3,"NCBI_gene"]


#GeneOntology
	go_ids<-sort(GOID_entrez$go_id[which(GOID_entrez$gene_id==data$NCBI_gene[i])])
	auxGO1<-character()
	if (length(go_ids)>0){	
		gos<-integer()
		for (z in 1:length(go_ids)) {
			gos[z]<-which(GOID_GOTERM$go_id==go_ids[z])[1]			
		}		
		auxGO1<-GOID_GOTERM$Term[gos]
	}

	auxGO2<-sort(gene2go$GO_term[which(gene2go$GeneID==data$NCBI_gene[i])])
	
	auxGO1_2<-unique(c(auxGO1,auxGO2))
	
	data$GeneOntology_BP[i]<-paste(auxGO1_2,collapse=";")


#GeneRIF
	sel5<-which(generifs$Gene_ID==data$NCBI_gene[i])
	data$GeneRIF[i]<-paste(generifs$GeneRIF_text[sel5],collapse=";")
	data$GeneRIF_PubMedID[i]<-paste(generifs$PubMed_ID[sel5],collapse=";")


#Interactions
	sel6<-which(interactions$gene_id==data$NCBI_gene[i])
	id.interactants<-unique(interactions$interactant_id[sel6])
	selentr<-integer()	
	for (f in 1:length(id.interactants)) {
		selentr<-c(selentr,which(HGNC_entrez$NCBI_gene==id.interactants[f]))
	}
	noms.interactants<-sort(HGNC_entrez$HGNC[selentr])
	data$Interactions[i]<-paste(noms.interactants,collapse=";")


#KEGG
	sel7<-which(KEGG_entrez$gene_or_orf_id==data$NCBI_gene[i])
	paths<-KEGG_entrez$pathway_id[sel7]
	selkegg<-integer()	
	for (k in 1:length(paths)){
		selkegg<-c(selkegg,which(KEGG_path$path_id==paths[k]))
	}
	noms.kegg<-sort(KEGG_path$path_name[selkegg])
	data$KEGG[i]<-paste(noms.kegg,collapse=";")


#Reactome
	sel8<-which(REACT_entrez$gene_id==data$NCBI_gene[i])
	paths_r<-REACT_entrez$DB_ID[sel8]
	selreact<-integer()
	for (r in 1:length(paths_r)){
		selreact<-c(selreact,which(REACT_path$DB_ID==paths_r[r]))
	}
	noms.react<-sort(REACT_path$path_name[selreact])
	data$Reactome[i]<-paste(noms.react,collapse=";")


#OMIM
	selOMIM<-which(HGNC_OMIM$HGNC==data$unique_name[i])[1]
	data$OMIM_ID[i]<-HGNC_OMIM$OMIM_ID[selOMIM]

	selOMIM_p<-which(OMIM_morbidmap$MIM_Number==data$OMIM_ID[i])
	data$OMIM_Phenotype[i]<-paste(OMIM_morbidmap$Phenotype[selOMIM_p],collapse=";")

#ProteinAtlas normal
	selPAnormal<-which(HGNC_ProteinAtlas_normal$Gene_name==data$unique_name[i])
	data$Expression_normal[i]<-paste(HGNC_ProteinAtlas_normal$Level[selPAnormal],collapse=";")


#ProteinAtlas cancer
	selPAcancer<-which(HGNC_ProteinAtlas_cancer$Gene_name==data$unique_name[i])
	data$Expression_cancer[i]<-paste(HGNC_ProteinAtlas_cancer$Level_mod[selPAcancer],collapse=";")


#ProteinAtlas subcellular location
	selPAsl<-which(HGNC_ProteinAtlas_sl$Gene_name==data$unique_name[i])
	auxPAsl<-c(HGNC_ProteinAtlas_sl$Validated[selPAsl],HGNC_ProteinAtlas_sl$Supported[selPAsl],HGNC_ProteinAtlas_sl$Approved[selPAsl],HGNC_ProteinAtlas_sl$Uncertain[selPAsl])
	auxPAsl<-auxPAsl[which(auxPAsl!="")]

	data$Subcellular_location[i]<-paste(auxPAsl,collapse=";")


	print(paste("Line ", i, " from ", nrow(data), ". First part.", sep=""))
}

#############################################################
	##########   LOOP FIRST PART (end)   ##########
#############################################################



################################################################################
	##########   Filtering Terms reading SECOND PART (start)   ##########
################################################################################
GO<-as.character(read.xls(paste(ruta, "filtering_terms_CRC.xls", sep="") ,sheet="GO",header=FALSE)[,1])

KEGG<-as.character(read.xls(paste(ruta, "filtering_terms_CRC.xls", sep="") ,sheet="KEGG",header=FALSE)[,1])

REACTOME<-as.character(read.xls(paste(ruta, "filtering_terms_CRC.xls", sep=""),sheet="REACTOME",header=FALSE)[,1])

GO_rm<-as.character(read.xls(paste(ruta, "filtering_terms_CRC.xls", sep=""),sheet="GO_rm",header=FALSE)[,1])

bib<-as.character(read.xls(paste(ruta, "filtering_terms_CRC.xls", sep="") ,sheet="bib",header=FALSE)[,1])
################################################################################
	##########   LFiltering Terms reading SECOND PART (end)   ##########
################################################################################


data$GeneOntology_BP<-as.character(data$GeneOntology_BP)
data$KEGG<-as.character(data$KEGG)
data$Reactome<-as.character(data$Reactome)
data$Function<-as.character(data$Function)
data$GeneRIF<-as.character(data$GeneRIF)


data$Function[which(data$Function=="")]<-NA
data$GeneRIF[which(data$GeneRIF=="")]<-NA
data$KEGG[which(data$KEGG=="")]<-NA
data$Reactome[which(data$Reactome=="")]<-NA

#New columns initialization
data$Filt.Function<-NA
data$Filt.GeneOntology_BP<-NA
data$Filt.GeneRIF<-NA
data$Filt.KEGG<-NA
data$Filt.Reactome<-NA


#############################################################
	##########   LOOP SECOND PART (start)   ##########
#############################################################
talls_go<-strsplit(data$GeneOntology_BP,";")
talls_kegg<-strsplit(data$KEGG,";")
talls_react<-strsplit(data$Reactome,";")
talls_generif<-strsplit(data$GeneRIF,";")


for (i in 1:nrow(data)) {

#GO
	tall.noms<-talls_go[[i]]
	resGO<-character()
	for (w in 1:length(GO)){
		term<-GO[w]
		aux<-tall.noms[grep(term,tall.noms)]
		if(length(aux)>0){
			for (z in 1:length(aux)){
				if ((aux[z] %in% GO_rm)==FALSE){
					resGO<-c(resGO,aux[z])
				}		
			}
		}
	}

	data$Filt.GeneOntology_BP[i]<-paste0(resGO,collapse=";")


#KEGG
	tall.noms<-talls_kegg[[i]]
	resKEGG<-character()
	for (k in 1:length(KEGG)){
		term<-KEGG[k]
		resKEGG<-c(resKEGG,tall.noms[grep(term,tall.noms)])
	}

	data$Filt.KEGG[i]<-paste0(resKEGG,collapse=";")


#Reactome
	tall.noms<-talls_react[[i]]
	resREACT<-character()
	for (r in 1:length(REACTOME)){
		term<-REACTOME[r]
		resREACT<-c(resREACT,tall.noms[grep(term,tall.noms)])
	}

	data$Filt.Reactome[i]<-paste0(resREACT,collapse=";")


#GeneRIF
	tall.noms<-talls_generif[[i]]
	resGENERIF<-character()
	for (g in 1:length(bib)){
		term<-bib[g]
		resGENERIF<-c(resGENERIF,tall.noms[grep(term,tall.noms)])
	}

	data$Filt.GeneRIF[i]<-paste0(resGENERIF,collapse=";")


#Function
	resFUNCTTION<-character()
	for (f in 1:length(bib)){
		term<-bib[f]
		if (length(grep(term,data$Function[i]))>0){
			resFUNCTTION<-c(resFUNCTTION,term)
		}	
	}

	data$Filt.Function[i]<-paste(resFUNCTTION,collapse=";")
	

	print(paste("Line ", i, " from ", nrow(data), ". Second part.", sep=""))

}

#############################################################
	##########   LOOP SECOND PART (end)   ##########
#############################################################

