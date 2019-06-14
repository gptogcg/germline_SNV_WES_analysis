# Germline SNV - WES Analysis
Germline SNV - WES analysis pipelines by PGaCCR lab.

_Folder structure changed February 26th 2019 but most workflows not modified already_


Two different models of inheritance are considered: autosomal dominant (heterozygosity) and recessive (homozygosity or compound heterozygosity). You can see an overview of the filtering criteria used in the [workflow diagrams](https://github.com/marcos-diazg/SNV_germline_WES_analysis/tree/master/Workflow_diagram).

Different workflows have also been developed according to the structure of the input files. One file per sample is allowed, as well as one big file containing the information of the germline variants of all the samples of the cohort. In all cases, tab separated text files are required to run these pipelines.

In all cases, an annotation folder will be required. This folder is not directly provided in GitHub because of size issue, but it should contain the following folders and files:
* DBs_annotation
  * filtering_terms_CRC.xls
  * HGNC_ProteinAtlas_cancer.csv
  * interactions
  * gene2go
  * HGNC_ProteinAtlas_normal_tissue.csv
  * OMIM_morbidmap.txt
  * generifs_basic
  * HGNC_ProteinAtlas_subcellular_location.csv
  * refSeqSummaryfilt_noalmohad.txt
  * HGNC_entrez.txt
  * HGNC_RefSeq.txt
  * uniprot-all.tab
  * HGNC_OMIM.txt
  * HGNC_synonyms.txt
* gene_families_&_pathways
  * Functional_Families_con_genes_hereditarios.txt
  * Gene_Families_con_genes_hereditarios.txt
  * PATHWAYS_db.txt
* genes_candidatos
  * genes_candidatos_23.09.2016.txt
  * Cancer_predisposing_genes_Nature_2014.txt
* LD_regions
  * LIST_Regions.txt
* pLI
  * fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt


## One file per sample
In this case two files will be required to run the pipeline:
1. __Directories file__ (directories.txt): Tabulated text file containing the absolute paths of the directories used for running the pipeline (input directory, output/results directory, workflow directory -the home directory by default- and annotation directory). It should contain 4 columns with the following mandatory headers: _Input_dir, Output_dir, Workflow_dir_ and _Annotation_dir_
2. __Files to analyze file__ (file_to_analyze.txt): Tabulated tect file containing in two columns the names of the files and samples to analysze. It should contain 2 columns with the following mandatory headers: _FILE_ and _NAME_.
