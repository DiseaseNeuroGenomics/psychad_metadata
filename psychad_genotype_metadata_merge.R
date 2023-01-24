########################################################################################################
# General functions and settings
########################################################################################################

{
  library(data.table)
  library(plyr)
  library(ggplot2)
  library(ggsci)
  library(ggpubr)
  library(plotly)
  library(synapser)
  library(UpSetR)
  library(stringr)
  synLogin('','')
  
  outDir = "/sc/arion/projects/roussp01a/jaro/project_psychAD/visualization/"
  
  ####################################
  # Minerva paths
  CLINICAL_METADATA = "/sc/arion/projects/psychAD/NPS-AD/freeze1_rc/metadata/syn26527784_latest.csv"                              # Minerva path for clinical metadata
  SNPARRAY_FINAL_VCF = "/sc/arion/projects/psychAD/genotypes/SNParray_PsychAD/genotypes.vcf.gz"                                   # Minerva path for "final" VCF file (containing only valid samples)
  SNPARRAY_FINAL_PLINK = "/sc/arion/projects/psychAD/genotypes/SNParray_PsychAD/genotypes"                                        # Minerva path for "final" plink file (containing only valid samples)
  ANEUPLOIDY_LIST = "/sc/arion/projects/roussp01a/karen/NPS_AD_imputed/phenotype_files/aneuploidy_samples_list.tsv"                            # List of samples with aneuploidy (by Karen)
  CALLED_ANCESTRY_METATADA = "/sc/arion/projects/roussp01a/karen/NPS_AD_imputed/phenotype_files/ancestry_assignment_multinom_3_PCs_final.tsv"  # Ancestry called from metadata (by Karen)
  CALLED_SEX_METADATA = "/sc/arion/projects/roussp01a/karen/NPS_AD_imputed/phenotype_files/sex_info.tsv"                                       # Sex called by plink (by Karen)
  CALLED_APOE_COMBINATIONS = "/sc/arion/projects/roussp01a/jaro/project_psychAD/gtcheck/datafiles/APOE_assignment_NPS_AD_microglia_phase2.txt" # Ancestry called by clustering approach from George's team (by Karen)
  FINAL_METADATA_PATH = "/sc/arion/projects/psychAD/NPS-AD/freeze1_rc/metadata/syn50919280_latest.csv"
  
  ####################################
  # Synapse paths
  SYNAPSE = list(root = "syn31827550")
}

{
  ## Load clinical metadata (containing SNParray IDs)
  clinical_metadata = read.csv(CLINICAL_METADATA)
  rownames(clinical_metadata) = clinical_metadata$SubID
  
  ## Get list of valid SNParray samples (i.e. all samples from "final" plink)
  valid_snparray = read.csv(paste0(SNPARRAY_FINAL_PLINK, ".fam"), sep="\t", header=F)
  valid_snparray$SubID = clinical_metadata[match(valid_snparray$V2, clinical_metadata$SNParray_PsychAD),"SubID"]
  valid_snparray = valid_snparray[,c("V2", "SubID")]
  colnames(valid_snparray) = c("ID", "SubID")
  print(paste0(">> SNParray :: Number of SNP-array samples surviving QC: ", nrow(valid_snparray)))
  
  ## Load imputed genotypes
  apoeFromSnparray = read.csv(CALLED_APOE_COMBINATIONS, sep="\t")
  rownames(apoeFromSnparray) = apoeFromSnparray$Formatted_ID
  apoeFromSnparray = apoeFromSnparray[valid_snparray$ID,]
  print(paste0(">> Numbers of samples with ", paste0(names(table(apoeFromSnparray$apo_alleles)), collapse=" ; ") , " ancestry: ", paste0(table(apoeFromSnparray$apo_alleles), collapse="/")))
  
  ## Load imputed ancestry
  ancestryFromSnpArray = read.csv(CALLED_ANCESTRY_METATADA, sep="\t")
  rownames(ancestryFromSnpArray) = ancestryFromSnpArray$Formatted_ID
  ancestryFromSnpArray = ancestryFromSnpArray[valid_snparray$ID,]
  colnames(ancestryFromSnpArray) = paste0("imp_anc_", colnames(ancestryFromSnpArray))
  print(paste0(">> Numbers of samples with ", paste0(names(table(ancestryFromSnpArray$superpop)), collapse="/") , " ancestry: ", paste0(table(apoeFromSnpArray$superpop), collapse="/")))
  
  ## Load imputed sex
  sexFromSnpArray = read.csv(CALLED_SEX_METADATA, sep="\t")
  colnames(sexFromSnpArray) = gsub("SubID", "Formatted_ID", gsub("^F$", "imp_sex_score", gsub("SubID", "Formatted_ID", gsub("^SNPSEX$", "imp_sex", colnames(sexFromSnpArray)))))
  sexFromSnpArray$imp_sex = gsub("^M$", "Male", gsub("^F$", "Female", sexFromSnpArray$imp_sex))
  rownames(sexFromSnpArray) = sexFromSnpArray$Formatted_ID
  sexFromSnpArray = sexFromSnpArray[valid_snparray$ID,]
  print(paste0(">> Numbers of samples with ", paste0(names(table(sexFromSnpArray$imp_sex)), collapse="/") , " ancestry: ", paste0(table(sexFromSnpArray$imp_sex), collapse="/")))
  
  metadata = do.call("cbind.data.frame", list(valid_snparray,
                                        sexFromSnpArray,
                                        ancestryFromSnpArray,
                                        apoeFromSnparray))
  metadata = df[,!endsWith(colnames(df), "Formatted_ID")]
}

########################################################################################################
# Save donor-level metadata to Synapse & Minerva
########################################################################################################

{
  # Save metadata to Minerva
  write.csv(metadata, file = file.path(outDir, "genotype_metadata.csv"), row.names=F)
  write.csv(metadata, file = FINAL_METADATA_PATH, row.names=F)
  
  # Save metadata to Synapse
  file = synStore(File(path=file.path(outDir, "genotype_metadata.csv"), parent=SYNAPSE$root), used="https://github.com/DiseaseNeuroGenomics/psychad_metadata/blob/main/psychad_genotype_metadata_merge.R")
  
  # Copy merge & querying code for metadata to public github
  cmd_copy_metadata = paste0("cp ", unlist(PSYCHAD_METADATA_CODES), " ", PSYCHAD_METADATA_GITHUB_MINERVA)
  sapply(cmd_copy_metadata, system)
  
  # Save genotype convertor
  write.csv(metadata[,c("SubID", "SNParray_HBBC", "SNParray_CommonMind", "WGS_CommonMind", "WGS_RUSH", "WGS_Ampad", "SNParray_Microglia", "SNParray_PsychAD", "ADSP_SampleId")], file=file.path(outDir, "NPSAD_gt_converter_Jaro.csv"), quote=F, row.names=F)
}
