########################################################################################################
# General functions and settings
########################################################################################################

{
  library(data.table)
  library(plyr)
  library(ggplot2)
  library(ggsci)
  library(plotly)
  library(synapser)
  library(UpSetR)
  library(stringr)
  #synLogin('your_login','your_password')
  synLogin('bendl','Ejhe8uca!')
  
  ####################################
  # Helper functions
  mpdf = function(x,width=7,height=7, ...)eval.parent(substitute({ pdf(paste0(outDir,"/plot_",make.names(x),".pdf"),useDingbats=F,width=width,height=height) })) # outDir must be defined as a global var
  #outDir = "/sc/arion/projects/roussp01a/jaro/project_psychAD/visualization/"
  outDir = "~/Desktop/project_psychAD/"
  dir.create(file.path(outDir, "plotly"), recursive=T)########################################################################################################

  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
  ####################################
  # Synapse paths
  SYNAPSE = list(root = "syn22399913",
                 individual_files = "syn22801001",
                 mssm_clinical_metadata = "syn26387767", 
                 hbcc_clinical_metadata = "syn23016097",
                 rush_clinical_metadata_1 = "syn35839674",
                 rush_clinical_metadata_2 = "syn26559760",
                 rush_clinical_metadata_3 = "syn26559745",
                 rush_clinical_metadata_4 = "syn27555634")
  
  ####################################
  # Minerva paths
  #RANDOMIZED_SAMPLES_LIST = "/sc/arion/projects/roussp01a/jaro/project_psychAD/psychAD-MSSM/cellranger_MULTIseq_pipeline/legacy/NPS_AD_rand_for_Shan.csv"   # Currently not used; list of samples that went through randomization step
  FIXSHEET = "/sc/arion/projects/roussp01a/jaro/project_psychAD/psychAD-MSSM/cellranger_MULTIseq_pipeline/legacy/MSSM_HBCC_RUSH_clinical_metadata_combo.csv" # Helper spreadsheet in which I manually set boolean flag (0/1) for all diagnosis (originally, it was plaintext)
  HARRY_UPDATED_LIST = "/sc/arion/projects/roussp01a/jaro/project_psychAD/psychAD-MSSM/Panos-NPS-Clinicals-1-12-2022.csv"    # Updated clinical metadata from MSBB that contains new variables such as ApoE
  GENOTYPE_METADATA_PATH = "/sc/arion/projects/roussp01a/jaro/project_psychAD/metadata/formatted_phenotypes_NPS_AD.tsv"      # Metadata for PsychAD SNParray; contains also default conversion to SubNum
  CUSTOM_AFTER_CHTCHECK_FIXES_MISC = "/sc/arion/projects/roussp01a/jaro/project_psychAD/metadata/fixes_mainTable.csv"        # Fixed of assignments of external genotypes to donor IDs 
  CUSTOM_AFTER_CHTCHECK_FIXES_SNRNASEQ = "/sc/arion/projects/roussp01a/jaro/project_psychAD/metadata/fixes_snRNAseq_new.csv" # assignment of snRNAseq samples to donors' SubIDs
  BLACKLISTED_SNPARRAYS = "/sc/arion/projects/roussp01a/jaro/project_psychAD/metadata/blacklisted_SNParray_new.csv"   # blacklisted SNParray samples (for PsychAD), usually due to not matching irrecoverable identity   
  BLACKLISTED_SNRNASEQ = "/sc/arion/projects/roussp01a/jaro/project_psychAD/metadata/blacklisted_snRNAseq_new.csv"    # blacklisted snRNAseq samples, usually due to gt contamination or completely unknown identity
  DEMUX_STATS = "/sc/arion/projects/roussp01a/jaro/project_psychAD/metadata/Allpools_w_rerun370_demux_stats.tsv"      # setting & output of vireo demultiplex (per-pool donors, genotypes, number of cells)
  FINAL_METADATA_PATH = "/sc/arion/projects/psychAD/NPS-AD/freeze1_rc/metadata/syn26527784_latest.csv"         # Minerva path for donor-level clinical metadata that is going to be compiled by this script
  FINAL_ALLINFO_PATH = "/sc/arion/projects/psychAD/NPS-AD/freeze1_rc/metadata/sample_to_donor_convertor.csv"   # Minerva path for snRNAseq_ID <-> SubID convertor that is going to be compiled by this script
  SNRNASEQ_GTCHECK_PATH = "/sc/arion/projects/psychAD/gtcheck_2/results/allInfo.csv"     # list of samples processed by the first vireo run (some of them were re-calculated by the second vireo run)
  SNRNASEQ_GTCHECK2_PATH = "/sc/arion/projects/psychAD/gtcheck_4/results/allInfo.csv"    # list of samples processed by the second vire run
  PRS_FROM_DEEPIKA = "/sc/arion/projects/roussp01a/deepika/PRS/psychAD/PRSCS_SCORES/EUR.PRS.scores.tsv"    # PRS estimates from Deepika
  CDR_DOMAINS = "/sc/arion/projects/roussp01a/jaro/project_psychAD/metadata/VA_CDR-Dom-10-2022.csv"        # CDR subdomains for MSSM samples that were additionally provided by Harry 
  CONVERTOR_TO_SUBID = "/sc/arion/projects/roussp01a/jaro/project_psychAD/metadata/convertor_to_SubID.csv" # Convertor of SubIDs to external genotypes (used as additional source of fixes for the initial dataset)
  PSYCHAD_METADATA_GITHUB_MINERVA = "/sc/arion/projects/roussp01a/jaro/github/psychad_metadata/"   # Minerva clone of github code for creating (this file) & querying clinical metadata (we will copy scripts there)
  
  # Paths to the code from Jaro's private bitbucket. We will copy code from merging & querying clinical metadata from there to his clone of psychad-metadata public github  
  PSYCHAD_METADATA_CODES = list(
    metadata_create = "/sc/arion/projects/roussp01a/jaro/atacseq_ad/NYGC_AD_R01/other_projects/psychad_clinical_metadata_merge.R",
    metadata_querying = "/sc/arion/projects/roussp01a/jaro/atacseq_ad/NYGC_AD_R01/other_projects/psychad_clinical_metadata_querying.R"
  )
  
  # Description of external genotypes
  dsetForComparison = list(
    "SNParray_HBBC" = list(name="SNParray_HBBC", 
                           genoPath = "/sc/arion/projects/roussp01a/jaro/project_psychAD/gtcheck/datafiles/hbcc_imputed_genotypes_hg38.annotated.vcf.gz", 
                           plink = "/sc/arion/projects/roussp01a/jaro/project_psychAD/gtcheck/datafiles/hbcc_imputed_genotypes_hg38.annotated", 
                           king = "/sc/arion/projects/roussp01a/jaro/project_psychAD/gtcheck/datafiles/hbcc_imputed_genotypes_hg38.annotated_king.kin", 
                           metadata = "/sc/arion/projects/roussp01a/jaro/project_psychAD/metadata/CMC_Human_SNP_mergedMetadata.csv", 
                           metadata_sampleMatchColumn = "SNP_report.Genotyping_Sample_ID", 
                           metadata_subjectMatchColumn = "Individual_ID",
                           convertor_subjectMatchColumn = "CMC_individual_ID"),
    "SNParray_CommonMind" = list(name="SNParray_CommonMind", 
                                 genoPath = "/sc/arion/projects/roussp01a/jaro/project_psychAD/gtcheck/datafiles/cmc_imputed_genotypes_hg38.annotated.vcf.gz", 
                                 plink = "/sc/arion/projects/roussp01a/jaro/project_psychAD/gtcheck/datafiles/cmc_imputed_genotypes_hg38.annotated", 
                                 king = "/sc/arion/projects/roussp01a/jaro/project_psychAD/gtcheck/datafiles/cmc_imputed_genotypes_hg38.annotated_king.kin", 
                                 metadata = "/sc/arion/projects/roussp01a/jaro/project_psychAD/metadata/CMC_Human_SNP_mergedMetadata.csv", 
                                 metadata_sampleMatchColumn = "SNP_report.Genotyping_Sample_ID", 
                                 metadata_subjectMatchColumn = "Individual_ID",
                                 convertor_subjectMatchColumn = "CMC_individual_ID"),
    "WGS_CommonMind" = list(name="WGS_CommonMind", 
                            genoPath = "/sc/arion/projects/roussp01a/jaro/project_psychAD/gtcheck/datafiles/CommonMind_WGS.vcf.gz", 
                            plink = "/sc/arion/projects/roussp01a/jaro/project_psychAD/gtcheck/datafiles/CommonMind_WGS", 
                            king = "/sc/arion/projects/roussp01a/jaro/project_psychAD/gtcheck/datafiles/CommonMind_WGS_king.kin",
                            metadata = "/sc/arion/projects/CommonMind/CMC-WGS-hg38/combo_by_deepika/total_CMC_Human_WGS_mergedMetadata.csv", 
                            metadata_sampleMatchColumn = "WGS_report.WGS_ID", #"Sample_DNA_ID", 
                            metadata_subjectMatchColumn = "Individual_ID",
                            convertor_subjectMatchColumn = "CMC_individual_ID"),
    "WGS_Ampad" = list(name="WGS_Ampad", 
                   genoPath = "/sc/arion/projects/roussp01a/jaro/project_psychAD/gtcheck/datafiles/AMPAD.vcf.gz", 
                   plink = "/sc/arion/projects/roussp01a/jaro/project_psychAD/gtcheck/datafiles/AMPAD", 
                   king = "/sc/arion/projects/roussp01a/jaro/project_psychAD/gtcheck/datafiles/AMPAD_king.kin",
                   metadata = "/sc/arion/projects/CommonMind/epigenAD/alz_analyses/meta-files/msbb.id.mapping.csv", 
                   metadata_full = "/sc/arion/projects/roussp01a/jaro/project_psychAD/metadata/MSBB_individual_metadata.csv",
                   metadata_sampleMatchColumn = "wgsID", 
                   metadata_subjectMatchColumn = "SynapseBrainID",
                   convertor_subjectMatchColumn = "AMPAD_msbb_individualID"),
    "SNParray_Microglia" = list(name = "SNParray_Microglia", 
                       #genoPath = "/sc/arion/projects/Microglia/genotyping/Combined_Phase1/Phase1_TOPMed_imputation/phase1_2_merge/Phase1.2.microglia.merged.TOPMed.dbSNP_v155.r2_0.3.CHR1..22.vcf.gz", 
                       genoPath = "/sc/arion/projects/Microglia/genotyping/Combined_Phase1and2/Microglia_Phase1and2_322ind_TopMed_r2_0p3_FixedNames.vcf.gz", 
                       plink = "/sc/arion/projects/psychAD/gtcheck/step2/files/Microglia_Phase1and2_322ind_TopMed_r2_0p3_FixedNames", 
                       king = "/sc/arion/projects/roussp01a/jaro/project_psychAD/gtcheck/datafiles/Microglia_king.kin",
                       metadata = "/sc/arion/projects/roussp01a/jaro/project_psychAD/gtcheck/datafiles/Microglia.csv", # TODO: check with Roman that this version is the most current 
                       metadata_sampleMatchColumn = "sample_name", 
                       metadata_subjectMatchColumn = "sample_name",
                       convertor_subjectMatchColumn = "Microglia_ID"),
    "WGS_RUSH" = list(name="WGS_RUSH", 
                      genoPath = "/sc/arion/projects/psychAD/genotypes/WGS_RUSH/mergeRush2.vcf.gz", 
                      plink = "/sc/arion/projects/psychAD/genotypes/WGS_RUSH/mergeRush2", 
                      king = "/sc/arion/projects/psychAD/genotypes/WGS_RUSH/mergeRush2_king.kin",
                      metadata = "/sc/arion/projects/psychAD/genotypes/WGS_RUSH/metadata/metadata.csv", 
                      metadata_sampleMatchColumn = "Sample_ID",
                      metadata_subjectMatchColumn = "projid",
                      convertor_subjectMatchColumn = "RUSH_ID"),
    "ADSP_WGS" = list(name="ADSP_WGS", 
                      genoPath = "/sc/arion/projects/psychAD/genotypes/ADSP/vcfs/gcad.preview.compact_filtered.r3.wgs.16906.GATK.2020.05.26.biallelic.genotypes.only_rosmap_mssm.vcf.gz", 
                      plink = "/sc/arion/projects/psychAD/genotypes/ADSP/vcfs/gcad.preview.compact_filtered.r3.wgs.16906.GATK.2020.05.26.biallelic.genotypes.only_rosmap_mssm", 
                      king = "/sc/arion/projects/psychAD/genotypes/ADSP/vcfs/gcad.preview.compact_filtered.r3.wgs.16906.GATK.2020.05.26.biallelic.genotypes.chr1_only_rosmap_mssm_king.kin", # FIXME: probably comparison with PsychAD SNParray
                      metadata = "/sc/arion/projects/psychAD/genotypes/ADSP/ADSP_WGS_17K_MetaData.csv", 
                      metadata_sampleMatchColumn = "SampleId",
                      metadata_subjectMatchColumn = "SampleId",
                      convertor_subjectMatchColumn = "ADSP_SampleId")
    #"SNParray_PsychAD" = list(name="SNParray_PsychAD", 
    #                          genoPath = "/sc/arion/projects/roussp01a/jaro/project_psychAD/gtcheck/datafiles/genotypes.vcf.gz", 
    #                          plink = "/sc/arion/projects/roussp01a/jaro/project_psychAD/gtcheck/datafiles/genotypes", 
    #                          king = "/sc/arion/projects/roussp01a/jaro/project_psychAD/gtcheck/datafiles/genotypes.kin",
    #                          metadata = "/sc/arion/projects/psychAD/gtcheck/results/allInfo.csv", 
    #                          #metadata = "/sc/arion/projects/psychAD/NPS-AD/freeze1_rc/metadata/syn26527784_latest.csv", 
    #                          metadata_sampleMatchColumn = "ID", #"", SNParray_PsychAD
    #                          metadata_subjectMatchColumn = "SubID", #"SubID",
    #                          convertor_subjectMatchColumn = "SubID"), #"CMC_individual_ID")
  )
  
  # Diagnosis classification
  DIAGNOSIS_CLASSIFICATION = list(
    "neurodegenerative" = c("AD", "MCI", "Dementia", "PD", "PD_uncertain_plus_encephalitic", "DLBD", "FTD", "ALS", "Others_Neurodegenerative"),
    "neurological" = c("MS", "PSP", "Epilepsy", "Seizures", "Tumor", "Migraine_headaches", "Head_Injury", "Vascular", "Others_Neurological"), 
    "neuropsychiatric" = c("SCZ", "MDD", "BD_unspecific", "BD_I", "BD_II", "PTSD", "ADHD", "OCD", "Tardive_Dyskinesia_Neuroleptic_induced", "Schizoaffective_bipolar", 
                           "Schizoaffective_depressive", "Anorexia", "Bulimia", "Anxiety", "Binge_Purge", "Eating_disorder", "Others_Neuropsychiatric"),
    "metabolic" = c("Diabetes_mellitus_unspecified", "TD_I", "TD_II"))
  DIAGNOSIS_BENIGN = c("Anxiety", "Migraine_headaches", DIAGNOSIS_CLASSIFICATION$metabolic)
  
  # SNParray IDs from Karen indicating samples with chromosome aneuploidy (note that only samples from PsychAD SNParray were checked)
  SEX_CHR_ANEUPLOIDY = c("214654", "173395", "120572", "214934", "44560", "192838", "34770", "182904", "201271", "214975", "191374", "214696", "214932", "208726")
  
  # Vocabulary for MSBB / neuropsychiatric symptoms
  MSSM_NPS = list(
    "MoodDysCurValue" = "Dysphoria",
    "DecIntCurValue" = "Anhedonia",
    "WtLossCurValue" = "Weight loss",
    "DecAppCurValue" = "Decreased appetite",
    "WtGainCurValue" = "Weight gain",
    "EarlyInsomCurValue" = "Early insomnia",
    "MidInsomCurValue" = "Middle insomnia",
    "LateInsomCurValue" = "Late insomnia",
    "HypersomCurValue" = "Hypersomnia",
    "PsychoAgiCurValue" = "Psychomotor agitation",
    "PsychoRetardCurValue" = "Psychomotor retardation",
    "FatCurValue" = "Fatigue",
    "WorthCurValue" = "Psychomotor retardation",
    "DelCurValue" = "Delusional worthlessness",
    "RumCurValue" = "Ruminations",
    "ThoughtDeathCurValue" = "Suicidal ideations",
    "Last2WkCurValue" = "Depression (current to 2 weeks)",
    "SixMoCurValue" = "Depression (current to 6 months)",
    "LifeCurValue" = "Depression (lifetime)"
  )
  
  # Lists of additional brain-specific neuropathological metrics available only for MSSM donors
  ADDIT_NEUROPAT_METRICS = list("Hipp" = c("HippoPlaquesValue", "HippoPlaquesWCoresValue", "HippoTanglesValue", "HippoLewyValue"),
    "Entor" = c("EntorPlaquesValue", "EntorPlaquesWCoresValue", "EntorTanglesValue", "EntorLewyValue"),
    "Amyg" = c("AmygPlaquesValue", "AmygPlaquesWCoresValue", "AmygTanglesValue", "AmygLewyValue"),	
    "Mid" = c("MidPlaquesValue", "MidPlaquesWCoresValue", "MidGliosisValue", "MidLewyValue"),
    "Sup" = c("SupPlaquesValue", "SupPlaquesWCoresValue", "SupGliosisValue", "SupLewyValue"),
    "Inf" = c("InfPlaquesValue", "InfPlaquesWCoresValue", "InfGliosisValue", "InfLewyValue"),
    "Occi" = c("OcciPlaquesValue", "OcciPlaquesWCoresValue", "OcciGliosisValue", "OcciLewyValue"),
    "Nucleus" = c("LewyNucleusValue"),
    "Locus" = c("LewyLocusValue"),
    "DorsalV" = c("LewyDorsalVValue"))
  
  # Vocabulary for MSBB / CERAD stages
  CERAD_DICT = list(
    "1" = "Normal brain", 
    "2" = "Definite Alzheimer's disease",
    "3" = "Probable Alzheimer's disease",
    "4" = "Possible Alzheimer's disease",
    "4a" = "Possible Alzheimer's disease",
    "4b" = "Possible Alzheimer's disease",
    "4c" = "Possible Alzheimer's disease",
    "5" = "Definite Parkinson's disease",
    "5a" = "Definite Parkinson's disease",
    "5b" = "Definite Parkinson's disease",
    "6" = "Uncertain Parkinson's disease",
    "7a" = "Vascular disease: infarct(s) only",
    "7b" = "Vascular disease: lacunar state only",
    "7c" = "Vascular disease: combined lacunar and large infractions",
    "7d" = "Vascular disease: Binswanger's disease",
    "7e" = "Vascular disease: Hemorrhage",
    "7f" = "Vascular disease: other",
    "8" = "Pick's disease (with Pick bodies)",
    "9" = "Lobar atrophy (without Pick bodies)",
    "10" = "CJD, spongiform encephalopathy",
    "11a" = "Down syndrome: Clinical Down's syndrome with AD pathology",
    "11b" = "Down syndrome: Clinical Down's syndrome without AD pathology",
    "12" = "Huntington disease",
    "13" = "Leukonecephaopathy",
    "14a" = "Tumor: primary",
    "14b" = "Tumor: secondary",
    "15" = "AIDS encephalopathy due to HIV infection",
    "16a" = "Others: other",
    "16b" = "Others: other",
    "16c" = "Others: other",
    "17" = "Lewy body with no AD changes",
    "18" = "Lewy body with AD changes"
  )
  
  # Vocabulary for RADC / cognitive diagnosis
  RUSH_cogdx = list(
    "1" = "NCI: No cognitive impairment (No impaired domains)",
    "2" = "MCI: Mild cognitive impairment (One impaired domain) and NO other cause of CI",
    "3" =  "MCI: Mild cognitive impairment (One impaired domain) AND another cause of CI",
    "4" = "AD: Alzheimer’s dementia and NO other cause of CI (NINCDS PROB AD)",
    "5" = "AD: Alzheimer’s dementia AND another cause of CI (NINCDS POSS AD)",
    "6" = "Other dementia: Other primary cause of dementia"
  )
  
  # Vocabulary for RADC / NIA-Reagan diagnosis of AD
  RUSH_ad_reagan = list(
    "1" = "AD present by NIA-Reagan pathology criteria (high or intermediate likelihood)",
    "0" = "AD not present by NIA-Reagan pathology criteria (low likelihood or no AD)"
  )
  
  # Vocabulary for RADC / chronic medical conditions and risk factors
  RUSH_med_con_sum_bl = list(
    "1" = "Hypertension",
    "2" = "Diabetes",
    "3" = "Heart disease",
    "4" = "Cancer",
    "5" = "Thyroid disease",
    "6" = "Head injury with loss of consciousness",
    "7" = "Stroke"
  )
  
  # Vocabulary for RADC / Lewy Body disease diagnosis
  RUSH_dlbdx = list(
    "0" = "not present",
    "1" = "nigral-predominant",
    "2" = "limbic-type",
    "3" = "neocortical-type"
  )
  
  # Vocabulary for RADC / Clinical Parkinson's disease
  RUSH_r_pd = list(
    "1" = "Highly Probable",
    "2" = "Probable",
    "3" = "Possible",
    "4" = "Not Present"
  )
  
  # Vocabulary for RADC / Cerebral atherosclerosis
  RUSH_cvda_4gp2 = list(
    "0" = "None",
    "1" = "Mild Small amounts in up to several arteries (typically less than 25% vessel involvement) without significant occlusion", 
    "2" = "Moderate In up to half of all visualized major arteries, with less than 50% occlusion of any single vessel", 
    "3" = "Severe In more than half of all visualized arteries, and/or more than 75% occlusion of one or more vessels")
  
  # Vocabulary for RADC / Cerebral amyloid angiopathy
  RUSH_caa_4gp = list(
    "0" = "None",
    "1" = "Mild",
    "2" = "Moderate",
    "3" = "Severe"
  )
  
  # Vocabulary for RADC / Presence of one or more gross chronic infarcts
  RUSH_ci_num2_gct = list(
    "0" = "No gross chronic Infarctions", 
    "1" = "One or more gross chronic infarctions (regardless of location)"
  )
  
  # Vocabulary for RADC / Parkinsoniasm
  RUSH_ci_num2_mct = list(
    "0" = "No chronic microinfarcts",
    "1" = "One or more chronic microinfarcts (regardless of location)"
  )
  
  # Vocabulary for RADC / Parkinsoniasm
  RUSH_dxpark = list(
    "1" = "Yes",
    "2" = "No",
    "8" = "Refusal",
    "9" = "Do not know"
  )
  
  # Vocabulary for RADC / ALS status
  RUSH_als = list( # not defined in codebook and no values provided for currently delivered samples (~ all) so who knows ...
    "1" = "Yes",
    "2" = "No"
  )
  
  # Vocabulary for RADC / Stroke status
  RUSH_cogdx_stroke = list("1" = "Yes", "2" = "No")
  
  # Vocabulary for RADC / Dementia classification
  DEMENTIA = list(
    "demceph" = list("1" = "Dementia due to communicating hydrocephalus", "0" = ""),
    "demcort" = list("1" = "Corticobasal Degeneration", "0" = ""),
    "demmotor" = list("1" = "Dementia due to Motor Neuron Disease", "0" = ""), 
    "demppa" = list("1" = "Dementia due to progressive aphasia", "0" = ""), 
    "demother" = list("1" = "Dementia due to Other condition", "0" = ""), 
    "demdelir" = list("1" = "Dementia due to Delirium", "0" = ""), 
    "demdep" = list("1" = "Dementia due to Depression", "0" = "")
  )
  
  # Vocabulary for RADC / Dementia classification (II)
  RUSH_demothsp = list("1" = "Dementia due to Other - specify condition", "0" = "")
  RUSH_demvasc = list("1" = "Dementia due to cerebrovascular disease", "0" = "")
  RUSH_dempark = list("1" = "Dementia due to Parkinson's disease", "0" = "")
  RUSH_demlewy = list("1" = "Dementia due to Lewy Body disease - Variant", "0" = "")
  RUSH_dempicks = list("1" = "Dementia due to Pick's - frontal lobe dementia", "0" = "")
  RUSH_dempsp = list("1" = "Dementia due to Progressive Supranuclear Palsy", "0" = "")
}

########################################################################################################
# Metadata load & fix
########################################################################################################

{
  #############################
  # Load MSSM clinical metadata & fix format of the data
  MSSM_clinical_meta = read.csv(synGet(entity=SYNAPSE$mssm_clinical_metadata)$path, header = T)
  colnames(MSSM_clinical_meta) = gsub("OcciLweyValue", "OcciLewyValue", colnames(MSSM_clinical_meta))
  
  for(clmn in names(MSSM_NPS)) {
    MSSM_clinical_meta[,paste0("nps_", clmn)] = ifelse(MSSM_clinical_meta[,clmn] %in% c(1,2), as.logical(gsub(1, TRUE, gsub(2, FALSE, MSSM_clinical_meta[,clmn]))), NA)
  }
  
  for(clmn in gsub("Cur", "Hx", names(MSSM_NPS))) {
    MSSM_clinical_meta[,paste0("nps_", clmn)] = ifelse(MSSM_clinical_meta[,clmn] %in% c(1,2), as.logical(gsub(1, TRUE, gsub(2, FALSE, MSSM_clinical_meta[,clmn]))), NA)
  }
  
  MSSM_clinical_meta$CERAD_1_text = unlist(sapply(MSSM_clinical_meta$CERAD_1, function(x) {
    if((nchar(x) > 0) & (x %in% names(CERAD_DICT))) { 
      return(CERAD_DICT[[x]])
    }
    return("")
  }))
  
  MSSM_clinical_meta$CERAD_2_text = unlist(sapply(MSSM_clinical_meta$CERAD_2, function(x) {
    if((nchar(x) > 0) & (x %in% names(CERAD_DICT))) { 
      return(CERAD_DICT[[x]])
    }
    return("")
  }))
  
  MSSM_clinical_meta$CERAD = sapply(1:nrow(MSSM_clinical_meta), function(i) {
    cerad1 = as.numeric(gsub("[abc]", "", MSSM_clinical_meta$CERAD_1)[i])
    cerad2 = as.numeric(gsub("[abc]", "", MSSM_clinical_meta$CERAD_2)[i])
    if(!is.na(cerad1) & (cerad1 %in% c(1:4))) {
      return(cerad1)
    } else if(!is.na(cerad2) & (cerad2 %in% c(1:4))) {
      return(cerad2)
    } else {
      return(NA)
    }
  })
  
  #############################
  # Load HBCC clinical metadata & fix format of the data
  HBCC_clinical_meta = read.csv(synGet(entity=SYNAPSE$rush_clinical_metadata_1)$path, header = T)
  HBCC_clinical_meta$Dx = paste0("Diagnosis: ", HBCC_clinical_meta$primaryDiagnosisDetail, " | death description: ", HBCC_clinical_meta$DescDeath)
  
  #############################
  # Load RUSH clinical metadata & fix format of the data
  RUSH_redundant_columns = c("scaled_to", "set1..Box", "set1.Position", "set1.Region", "set1.Shelf", "set1.Row", "set1.Intra.box", "set2..Box", "set2.Position", "set2.Region", "set2.Freezer", "set2.Shelf", "set2.Row", "set2.Intra.box")
  RUSH_clinical_meta_original = read.csv(synGet(entity=SYNAPSE$hbcc_clinical_metadata)$path, header = T)
  RUSH_clinical_meta_original_fixes = read.csv(synGet(entity=SYNAPSE$rush_clinical_metadata_2)$path, header = T)
  columns_to_be_fixed = intersect(colnames(RUSH_clinical_meta_original), colnames(RUSH_clinical_meta_original_fixes))
  RUSH_clinical_meta_original[match(RUSH_clinical_meta_original_fixes$projid, RUSH_clinical_meta_original$projid),columns_to_be_fixed] = RUSH_clinical_meta_original_fixes[,columns_to_be_fixed]
  columns_to_be_added = setdiff(colnames(RUSH_clinical_meta_original_fixes), colnames(RUSH_clinical_meta_original))
  RUSH_clinical_meta_original[,columns_to_be_added] = NA
  RUSH_clinical_meta_original[match(RUSH_clinical_meta_original_fixes$projid, RUSH_clinical_meta_original$projid),columns_to_be_added] = RUSH_clinical_meta_original_fixes[,columns_to_be_added]
  RUSH_clinical_meta_original = RUSH_clinical_meta_original[,!colnames(RUSH_clinical_meta_original) %in% RUSH_redundant_columns]
  
  RUSH_clinical_meta_additional = read.csv(synGet(entity=SYNAPSE$rush_clinical_metadata_3)$path, header = T)
  RUSH_clinical_meta_additional_2 = read.csv(synGet(entity=SYNAPSE$rush_clinical_metadata_4)$path, header = T)
  RUSH_intersectCols = intersect(intersect(colnames(RUSH_clinical_meta_original), colnames(RUSH_clinical_meta_additional)), colnames(RUSH_clinical_meta_additional_2))
  RUSH_unionCols = union(union(colnames(RUSH_clinical_meta_original), colnames(RUSH_clinical_meta_additional)), colnames(RUSH_clinical_meta_additional_2))
  
  RUSH_clinical_meta_original[,setdiff(RUSH_unionCols, colnames(RUSH_clinical_meta_original))] = NA
  RUSH_clinical_meta_additional[,setdiff(RUSH_unionCols, colnames(RUSH_clinical_meta_additional))] = NA
  RUSH_clinical_meta_additional_2[,setdiff(RUSH_unionCols, colnames(RUSH_clinical_meta_additional_2))] = NA
  RUSH_clinical_meta = data.frame(rbind(RUSH_clinical_meta_original[!(RUSH_clinical_meta_original$projid %in% RUSH_clinical_meta_additional$projid) & !(RUSH_clinical_meta_original$projid %in% RUSH_clinical_meta_additional_2$projid),RUSH_unionCols],
                                        #RUSH_clinical_meta_additional[!(RUSH_clinical_meta_additional$projid %in% RUSH_clinical_meta_additional_2$projid),RUSH_unionCols],
                                        RUSH_clinical_meta_additional_2[,RUSH_unionCols]))
  dim(RUSH_clinical_meta)
  RUSH_clinical_meta[RUSH_clinical_meta == "NULL"] = NA
  
  RUSH_clinical_meta$cogdx_text = unlist(sapply(RUSH_clinical_meta$cogdx, function(x) {
    if((nchar(x) > 0) & (x %in% names(RUSH_cogdx))) { 
      return(RUSH_cogdx[[x]])
    }
    return("")
  }))
  RUSH_clinical_meta$ad_reagan_text = unlist(sapply(as.character(RUSH_clinical_meta$ad_reagan), function(x) {
    if((nchar(x) > 0) & (x %in% names(RUSH_ad_reagan))) { 
      return(RUSH_ad_reagan[[x]])
    }
    return("")
  }))
  RUSH_clinical_meta$med_con_sum_bl_text = unlist(sapply(RUSH_clinical_meta$med_con_sum_bl, function(x) {
    if((nchar(x) > 0) & (x %in% names(RUSH_med_con_sum_bl))) { 
      return(RUSH_med_con_sum_bl[[x]])
    }
    return("")
  }))
  RUSH_clinical_meta$dlbdx_text = unlist(sapply(as.character(RUSH_clinical_meta$dlbdx), function(x) {
    if((nchar(x) > 0) & (x %in% names(RUSH_dlbdx))) { 
      return(RUSH_dlbdx[[x]])
    }
    return("")
  }))
  
  RUSH_clinical_meta$caa_4gp_text = unlist(sapply(as.character(RUSH_clinical_meta$caa_4gp), function(x) {
    if((nchar(x) > 0) & (x %in% names(RUSH_caa_4gp))) { 
      return(RUSH_caa_4gp[[x]])
    }
    return("")
  }))
  
  RUSH_clinical_meta$cvda_4gp2_text = unlist(sapply(as.character(RUSH_clinical_meta$cvda_4gp2), function(x) {
    if((nchar(x) > 0) & (x %in% names(RUSH_cvda_4gp2))) { 
      return(RUSH_cvda_4gp2[[x]])
    }
    return("")
  }))
  
  RUSH_clinical_meta$ci_num2_gct_text = unlist(sapply(as.character(RUSH_clinical_meta$ci_num2_gct), function(x) {
    if((nchar(x) > 0) & (x %in% names(RUSH_ci_num2_gct))) { 
      return(RUSH_ci_num2_gct[[x]])
    }
    return("")
  }))
  
  RUSH_clinical_meta$ci_num2_mct_text = unlist(sapply(as.character(RUSH_clinical_meta$ci_num2_mct), function(x) {
    if((nchar(x) > 0) & (x %in% names(RUSH_ci_num2_mct))) { 
      return(RUSH_ci_num2_mct[[x]])
    }
    return("")
  }))
  
  RUSH_clinical_meta$dxpark_text = unlist(sapply(as.character(RUSH_clinical_meta$dxpark), function(x) {
    if((nchar(x) > 0) & (x %in% names(RUSH_dxpark))) { 
      return(RUSH_dxpark[[x]])
    }
    return("")
  }))
  
  RUSH_clinical_meta$cogdx_stroke_text = unlist(sapply(as.character(RUSH_clinical_meta$cogdx_stroke), function(x) {
    if((nchar(x) > 0) & (x %in% names(RUSH_cogdx_stroke))) { 
      return(RUSH_cogdx_stroke[[x]])
    }
    return("")
  }))
}

########################################################################################################
# Metadata merging across institutions
########################################################################################################

{
  # Make plaintext AD column; for MSSM, it's from Harry's Dx column; for HBCC, it's from HBCC "primary diagnosis detail" / "description of death" columns; for RAD, it's from multiple diagnosis-related columns
  HBCC_clinical_meta$Dx = paste0("Diagnosis: ", HBCC_clinical_meta$primaryDiagnosisDetail, " | death description: ", HBCC_clinical_meta$DescDeath)
  MSSM_clinical_meta$Dx = paste0("Harry's plaintext: ", MSSM_clinical_meta$Final.Dx, " | CERAD_1: ", MSSM_clinical_meta$CERAD_1_text, " | CERAD_2: ", MSSM_clinical_meta$CERAD_2_text)
  RUSH_clinical_meta$Dx = paste0("cogdx: ", RUSH_clinical_meta$cogdx, " | ad_reagan: ", RUSH_clinical_meta$ad_reagan_text, " | dlbdx: ", RUSH_clinical_meta$dlbdx_text, " | med.con: ", RUSH_clinical_meta$med_con_sum_bl_text, " | stroke: ", RUSH_clinical_meta$cogdx_stroke_text)
  
  # Create SubID-like values: 'M' prefix for MSSM / "SubNum"; 'H' prefix for HBCC / "Brain ID"; 'R' prefix for RADC / "projid"
  MSSM_clinical_meta$SubNum = MSSM_clinical_meta$SubID
  HBCC_clinical_meta$Brain.ID = paste0('H', HBCC_clinical_meta$Brain.ID)
  RUSH_clinical_meta$projid = paste0('R', RUSH_clinical_meta$projid)
  
  # Shared clinical info: Age, Sex, Ethnicity, Education, PMI, Dx, pH, CERAD. Note: MSSM has no education info; RUSH has no pH; HBCC has no CERAD and no BRAAK.
  MSSM_shared = data.frame(SubID = MSSM_clinical_meta$SubID, Brain_bank = "MSSM", Age = MSSM_clinical_meta$Age, Sex = MSSM_clinical_meta$Sex, Ethnicity = MSSM_clinical_meta$Race, Education = NA, PMI = MSSM_clinical_meta$PMI/60, Dx = MSSM_clinical_meta$Dx, pH = MSSM_clinical_meta$pH, CERAD = MSSM_clinical_meta$CERAD, BRAAK = MSSM_clinical_meta$Braak.Stage, PLAQUE = MSSM_clinical_meta$Plq_Mn)
  MSSM_shared$CERAD = with(MSSM_shared, ifelse(CERAD == "1", 1, ifelse(CERAD == "2", 4, ifelse(CERAD == "3", 3, ifelse(CERAD == "4", 2, NA)))))
  MSSM_shared = MSSM_shared[,!colnames(MSSM_shared) %in% "CERAD2"]
  MSSM_unique = MSSM_clinical_meta[,!colnames(MSSM_clinical_meta) %in% c("Age", "Sex", "Race", "PMI", "Final.Dx", "pH", "Braak.Stage")]

  # HBCC traits
  HBCC_shared = data.frame(SubID = HBCC_clinical_meta$Brain.ID, Brain_bank = "HBCC", Age = HBCC_clinical_meta$ageOfDeath, Sex = HBCC_clinical_meta$Reported.Gender, Ethnicity = HBCC_clinical_meta$Ethnicity, Education = HBCC_clinical_meta$Years.of.Education, PMI = HBCC_clinical_meta$PMI..in.hours., Dx = HBCC_clinical_meta$Dx, pH = HBCC_clinical_meta$pH, CERAD = NA, BRAAK = NA, PLAQUE = NA)
  HBCC_shared$Ethnicity = with(HBCC_shared, ifelse(Ethnicity == "Caucasian", "White", ifelse(Ethnicity == "African-American", "Black", ifelse(Ethnicity == "Hispanic", "Hispanic", ifelse(Ethnicity == "Asian", "Asian", NA)))))
  HBCC_unique = HBCC_clinical_meta[,!colnames(HBCC_clinical_meta) %in% c("ageOfDeath", "Reported.Gender", "Ethnicity", "Education", "PMI..in.hours.", "Dx", "pH")]
  colnames(HBCC_unique) = gsub("Brain.ID","SubID",colnames(HBCC_unique))
  
  # RUSH traits
  RUSH_shared = data.frame(SubID = RUSH_clinical_meta$projid, Brain_bank = "RUSH", Age = RUSH_clinical_meta$age_death, Sex = RUSH_clinical_meta$msex, Ethnicity = RUSH_clinical_meta$race7, Education = RUSH_clinical_meta$educ, PMI = RUSH_clinical_meta$pmi, Dx = RUSH_clinical_meta$Dx, pH = NA, CERAD = RUSH_clinical_meta$ceradsc, BRAAK = RUSH_clinical_meta$braaksc, PLAQUE = RUSH_clinical_meta$plaq_n)
  RUSH_shared$Sex2 = with(RUSH_shared, ifelse(Sex == 0, "Female", ifelse(Sex == 1, "Male", NA)))
  RUSH_shared$Ethnicity2 = with(RUSH_shared, ifelse(Ethnicity == 1, "White", ifelse(Ethnicity == 2, "Black", ifelse(Ethnicity == 3, "American Indian or Alaska Native", ifelse(Ethnicity == 5, "Asian", ifelse(Ethnicity == 7, "Unknown", NA))))))
  RUSH_shared$CERAD2 = with(RUSH_shared, ifelse(CERAD == "4", 1, ifelse(CERAD == "3", 2, ifelse(CERAD == "2", 3, ifelse(CERAD == "1", 4, NA)))))
  RUSH_shared$Sex = RUSH_shared$Sex2
  RUSH_shared$Ethnicity = RUSH_shared$Ethnicity2
  RUSH_shared$Dx2 = NA
  RUSH_shared$CERAD = RUSH_shared$CERAD2
  RUSH_shared = RUSH_shared[,!colnames(RUSH_shared) %in% c("Sex2","Ethnicity2","Dx2","CERAD2","braaksc")]
  RUSH_unique = RUSH_clinical_meta[,!colnames(RUSH_clinical_meta) %in% c("age_death","msex","race7","educ","pmi","ceradsc","braaksc", "Dx")]
  colnames(RUSH_unique) = gsub("projid","SubID",colnames(RUSH_unique))
  
  ##############ffixsheet###############
  # Load & apply fixsheet (manual assignment of donors to diagnosis)
  fixsheet = read.csv(FIXSHEET)
  fixsheet = fixsheet[,c("SubID", "Comment", unlist(DIAGNOSIS_CLASSIFICATION))]

  for(cname in colnames(fixsheet)[!colnames(fixsheet) %in% c("SubID", "Comment")]) {
    fixsheet[,cname] = as.numeric(ifelse(is.na(fixsheet[,cname]), F, T))
  }
  rownames(fixsheet) = fixsheet$SubID
  fixsheet$Comment[is.na(fixsheet$Comment)] = ""
  
  #############################
  # Combine MSSM, HBCC and RUSH shared traits
  shared_combo = rbind(unique(MSSM_shared), unique(HBCC_shared), unique(RUSH_shared))
  combo1 = join(shared_combo, unique(MSSM_unique), by = "SubID")
  combo2 = join(combo1, unique(HBCC_unique), by = "SubID")
  combo3 = join(combo2, unique(RUSH_unique), by = "SubID")
  MSSM_HBCC_RUSH_clinical_metadata_combo = combo3
  MSSM_HBCC_RUSH_clinical_metadata_combo = MSSM_HBCC_RUSH_clinical_metadata_combo[!duplicated(MSSM_HBCC_RUSH_clinical_metadata_combo$SubID),]
  MSSM_HBCC_RUSH_clinical_metadata_combo = MSSM_HBCC_RUSH_clinical_metadata_combo[,!colnames(MSSM_HBCC_RUSH_clinical_metadata_combo) %in% c("Institution", "race")]
  
  metadata = cbind.data.frame(MSSM_HBCC_RUSH_clinical_metadata_combo, fixsheet[MSSM_HBCC_RUSH_clinical_metadata_combo$SubID,])
  metadata$is_problematic = sapply(metadata$Comment, function(x) nchar(x) > 0)
  rownames(metadata) = metadata$SubID
  
  #############################
  # Load & apply Harry's updated info (containing also ApoE, pH and few other columns)
  updatedHarry = read.csv(HARRY_UPDATED_LIST)
  
  metadata$BRAAK_old = metadata$BRAAK
  metadata$BRAAK_PD = NA
  
  metadata[metadata$Brain_bank=="MSSM", "BRAAK"] = as.integer(gsub("-1", NA, gsub("9", NA, updatedHarry[match(metadata[metadata$Brain_bank=="MSSM", "SubID"],updatedHarry$SubID), "BRAAK_Alz"])))
  metadata[which((metadata$Brain_bank=="MSSM") & is.na(metadata$BRAAK) & (metadata$CERAD==1)),"BRAAK"] = 0
  metadata[metadata$Brain_bank=="MSSM", "BRAAK_PD"] = updatedHarry[match(metadata[metadata$Brain_bank=="MSSM", "SubID"],updatedHarry$SubID), "BRAAK_Park"]
  metadata[metadata$Brain_bank=="MSSM", "pH"] = updatedHarry[match(metadata[metadata$Brain_bank=="MSSM", "SubID"],updatedHarry$SubID), "pH"]
  metadata[metadata$Brain_bank=="MSSM", "apoe_genotype"] = gsub("/", "", updatedHarry[match(metadata[metadata$Brain_bank=="MSSM", "SubID"],updatedHarry$SubID), "ApoE"])
  
  #############################
  # Add link to PsychAD SNParray
  GENOTYPE_METADATA = read.csv(GENOTYPE_METADATA_PATH, sep="\t")
  GENOTYPE_METADATA = GENOTYPE_METADATA[!is.na(GENOTYPE_METADATA$SubNum),]
  GENOTYPE_METADATA$SubID = sapply(1:nrow(GENOTYPE_METADATA), function(i) paste0(substring(GENOTYPE_METADATA$cohort[i], 1, 1), GENOTYPE_METADATA$SubNum[i]))
  metadata$SNParray_PsychAD = GENOTYPE_METADATA[match(metadata$SubID, GENOTYPE_METADATA$SubID),"Formatted_ID"]
  
  #############################
  # Make additional fixes for RADC samples
  metadata[metadata$Brain_bank=="RUSH", "PD"] = sapply(1:sum(metadata$Brain_bank=="RUSH"), function(i) { 
    ifelse((metadata[metadata$Brain_bank=="RUSH","dxpark_text"][i] == "Yes"), 1, 0 ) })
  metadata[metadata$Brain_bank=="RUSH", "Vascular"] = sapply(1:sum(metadata$Brain_bank=="RUSH"), function(i) { 
    ifelse(((metadata[metadata$Brain_bank=="RUSH", "cogdx_stroke_text"][i] == "Yes") | 
              (metadata[metadata$Brain_bank=="RUSH", "demvasc"][i] == "1")), 1, 0 ) })
  metadata[metadata$Brain_bank=="RUSH", "Diabetes_mellitus_unspecified"] = sapply(1:sum(metadata$Brain_bank=="RUSH"), function(i) { 
    ifelse(((metadata[metadata$Brain_bank=="RUSH", "med_con_sum_bl"][i] == "2")), 1, 0 ) })
  metadata[metadata$Brain_bank=="RUSH", "Head_Injury"] = sapply(1:sum(metadata$Brain_bank=="RUSH"), function(i) { 
    ifelse(((metadata[metadata$Brain_bank=="RUSH", "headinjrloc_bl"][i] == "1")), 1, 0 ) })
  metadata[metadata$Brain_bank=="RUSH", "PSP"] = sapply(1:sum(metadata$Brain_bank=="RUSH"), function(i) { 
    ifelse((metadata[metadata$Brain_bank=="RUSH", "dempsp"][i] == "1"), 1, 0 ) })
  metadata[metadata$Brain_bank=="RUSH", "FTD"] = sapply(1:sum(metadata$Brain_bank=="RUSH"), function(i) { 
    ifelse((metadata[metadata$Brain_bank=="RUSH", "dempicks"][i] == "1"), 1, 0 ) })
  metadata[metadata$Brain_bank=="RUSH", "DLBD"] = sapply(1:sum(metadata$Brain_bank=="RUSH"), function(i) { 
    ifelse((metadata[metadata$Brain_bank=="RUSH", "demlewy"][i] == "1"), 1, 0 ) })
  metadata[metadata$Brain_bank=="RUSH", "ALS"] = sapply(1:sum(metadata$Brain_bank=="RUSH"), function(i) { 
    ifelse((metadata[metadata$Brain_bank=="RUSH", "als"][i] == "1"), 1, 0 ) }) # codebook doesn't define als so i am just guessing ... 
  
  metadata[metadata$Brain_bank=="RUSH", "Dx_text"] = sapply(which(metadata$Brain_bank=="RUSH"), function(i)
    paste0(na.omit(sapply(c("demceph", "demcort", "demmotor", "demppa", "demother", "demdelir", "demdep"), function(x) {
      #print(paste0(i, "-", x))
      z = metadata[i,x]
      if(is.na(z) | is.null(z) | z==0)
        return(NA) else {
          out = DEMENTIA[[x]][[as.character(z)]]
          if(x == "demother") {
            out = paste0(out, ":", metadata[i,"demothsp"])
          }
        }
    })), collapse = "; "))
  
  metadata[which(grepl("brain injury", metadata$demothsp)), "Head_Injury"] = 1
  
  metadata[metadata$Brain_bank=="RUSH", "DLBD"] = sapply(1:sum(metadata$Brain_bank=="RUSH"), function(i) {
    ifelse((metadata[metadata$Brain_bank=="RUSH", "demlewy"][i] == "1"), 1, 0 ) })
  
  metadata$AD = NA
  metadata$MCI = NA
  metadata$Dementia = NA
  metadata[metadata$Brain_bank=="RUSH", "MCI"] = sapply(1:sum(metadata$Brain_bank=="RUSH"), function(i) { 
    ifelse((metadata[metadata$Brain_bank=="RUSH", "cogdx"][i] %in% c("2", "3")), 1, ifelse(is.na(metadata[metadata$Brain_bank=="RUSH", "cogdx"][i]), NA, 0) ) })
  
  metadata[metadata$Brain_bank=="RUSH", "Dementia"] = sapply(1:sum(metadata$Brain_bank=="RUSH"), function(i) { 
    ifelse((metadata[metadata$Brain_bank=="RUSH", "cogdx"][i] %in% c("4", "5", "6")), 1, ifelse(is.na(metadata[metadata$Brain_bank=="RUSH", "cogdx"][i]), NA, 0) ) })
  
  metadata[metadata$Brain_bank=="RUSH", "AD"] = sapply(1:sum(metadata$Brain_bank=="RUSH"), function(i) { 
    ifelse(((metadata[metadata$Brain_bank=="RUSH", "CERAD"][i] %in% c(2, 3, 4)) & (metadata[metadata$Brain_bank=="RUSH", "BRAAK"][i] %in% c(3, 4, 5, 6)) & (metadata[metadata$Brain_bank=="RUSH","MCI"][i] | metadata[metadata$Brain_bank=="RUSH","Dementia"][i])) | 
             ((metadata[metadata$Brain_bank=="RUSH", "CERAD"][i] %in% c(2, 3, 4)) & (metadata[metadata$Brain_bank=="RUSH", "BRAAK"][i] %in% c(3, 4, 5, 6)) & (is.na(metadata[metadata$Brain_bank=="RUSH","MCI"][i]) & is.na(metadata[metadata$Brain_bank=="RUSH","Dementia"][i]))) & 
             (metadata[metadata$Brain_bank=="RUSH", "PD"][i] == 0), 1, 0) })
  
  metadata[metadata$Brain_bank=="MSSM", "MCI"] = sapply(1:sum(metadata$Brain_bank=="MSSM"), function(i) { 
    ifelse((metadata[metadata$Brain_bank=="MSSM", "CDRScore"][i] %in% c("0.5")), 1, ifelse(is.na(metadata[metadata$Brain_bank=="MSSM", "CDRScore"][i]), NA, 0) ) })
  
  metadata[metadata$Brain_bank=="MSSM", "Dementia"] = sapply(1:sum(metadata$Brain_bank=="MSSM"), function(i) { 
    ifelse((metadata[metadata$Brain_bank=="MSSM", "CDRScore"][i] %in% c("1", "2", "3", "4", "5")), 1, ifelse(is.na(metadata[metadata$Brain_bank=="MSSM", "CDRScore"][i]), NA, 0) ) })
  
  metadata[metadata$Brain_bank=="MSSM", "AD"] = sapply(1:sum(metadata$Brain_bank=="MSSM"), function(i) { 
    ifelse((((metadata[metadata$Brain_bank=="MSSM", "CERAD"][i] %in% c(2, 3, 4)) & (metadata[metadata$Brain_bank=="MSSM", "BRAAK"][i] %in% c(3, 4, 5, 6)) & (metadata[metadata$Brain_bank=="MSSM", "CDRScore"][i] %in% c(0.5, 1, 2, 3, 4, 5))) | 
             (is.na(metadata[metadata$Brain_bank=="MSSM", "CERAD"][i]) & (metadata[metadata$Brain_bank=="MSSM", "BRAAK"][i] %in% c(3, 4, 5, 6)) & (metadata[metadata$Brain_bank=="MSSM", "CDRScore"][i] %in% c(0.5, 1, 2, 3, 4, 5))) | 
             ((metadata[metadata$Brain_bank=="MSSM", "CERAD"][i] %in% c(2, 3, 4)) & is.na(metadata[metadata$Brain_bank=="MSSM", "BRAAK"][i]) & (metadata[metadata$Brain_bank=="MSSM", "CDRScore"][i] %in% c(0.5, 1, 2, 3, 4, 5))) | 
             ((metadata[metadata$Brain_bank=="MSSM", "CERAD"][i] %in% c(2, 3, 4)) & (metadata[metadata$Brain_bank=="MSSM", "BRAAK"][i] %in% c(3, 4, 5, 6)) & is.na(metadata[metadata$Brain_bank=="MSSM", "CDRScore"][i])))
           & (metadata[metadata$Brain_bank=="MSSM", "PD"][i] == 0) & (metadata[metadata$Brain_bank=="MSSM", "PD_uncertain_plus_encephalitic"][i] == 0),
           1, 0) })
}

########################################################################################################
# Additional custom changes that resulted from gtcheck (after psychAD_gtcheck1/2/3/4.R)
########################################################################################################

{
  # Remove blacklisted SNParrays 
  blacklisted_snparray = read.csv(BLACKLISTED_SNPARRAYS)
  metadata$SNParray_PsychAD = ifelse(metadata$SNParray_PsychAD %in% blacklisted_snparray$ID, NA, metadata$SNParray_PsychAD)
  
  # Make custom changes, mostly assigning different SNParrays (or CMC/AMPAD/HBCC/Microglia IDs) to donors (SubIDs)
  convertToSubID = read.csv(CONVERTOR_TO_SUBID, stringsAsFactors=F)
  convertToSubID$RUSH_ID = str_pad(convertToSubID$RUSH_ID, 8, pad = "0")
  rownames(convertToSubID) = convertToSubID$SubID
  metadata[,c("CMC_individual_ID", "AMPAD_msbb_individualID", "Microglia_ID", "RUSH_ID", "ADSP_SampleId")] = NA
  for(i in 1:nrow(convertToSubID)) {
    if(!rownames(convertToSubID)[i] %in% metadata$SubID) {
      next
    }
    for(gtfix in c("CMC_individual_ID", "AMPAD_msbb_individualID", "Microglia_ID", "RUSH_ID", "ADSP_SampleId")) {
      metadata[convertToSubID$SubID[i],gtfix] = convertToSubID[i,gtfix]
    }
  }
  
  fixesTable = read.csv(CUSTOM_AFTER_CHTCHECK_FIXES_MISC)
  for(i in 1:nrow(fixesTable)) {
    fx = fixesTable[i,]
    if(fx$SubID %in% metadata$SubID) {
      if(fx$columnName %in% c("ADSP_SampleId", "AMPAD_msbb_individualID", "CMC_individual_ID", "Microglia_ID", "RUSH_ID", "SNParray_PsychAD")) {
        metadata[which(metadata[,fx$columnName]==fx$value),fx$columnName] = NA
      }
      metadata[fx$SubID, fx$columnName] = ifelse(fx$value == "", NA, fx$value)
    }
  }
}

########################################################################################################
# Add links to snRNAseq libraries
########################################################################################################

{
  allInfo1 = read.csv(SNRNASEQ_GTCHECK_PATH)
  rownames(allInfo1) = allInfo1$Sample_ID
  allInfo2 = read.csv(SNRNASEQ_GTCHECK2_PATH)
  rownames(allInfo2) = allInfo2$Sample_ID
  
  table(rownames(allInfo2) %in% rownames(allInfo1))
  sharedCols = intersect(colnames(allInfo1), colnames(allInfo2))
  allInfo = rbind.data.frame(allInfo1[!allInfo1$poolName %in% allInfo2$poolName,sharedCols], allInfo2[,sharedCols])
  
  fixesTable2 = read.csv(CUSTOM_AFTER_CHTCHECK_FIXES_SNRNASEQ)
  for(i in 1:nrow(fixesTable2)) {
    fx = fixesTable2[i,]
    allInfo[fx$ID,c("SubNum", "SubID")] = c(fx$SubNum, fx$SubID)
  }
  
  # Remove blacklisted snRNAseq libraries
  blacklisted_snrnaseq = read.csv(BLACKLISTED_SNRNASEQ)
  allInfo = allInfo[!(rownames(allInfo) %in% blacklisted_snrnaseq$ID),]
  allInfo = allInfo[sapply(allInfo$SubID, function(x) { ifelse(startsWith(x, prefix="donor"), F, T) } ),]
  
  metadata$snRNAseq_ID = sapply(metadata$SubID, function(subID) {
    paste0(allInfo[(allInfo$SubID == subID),"Sample_ID"], collapse=";")
  })
  metadata$snRNAseq_ID_count = sapply(metadata$SubID, function(subID) {
    nrow(allInfo[(allInfo$SubID == subID),])
  })
  
  #metadata = metadata[metadata$SubID %in% fixsheet$SubID,]
  
  # Remove donors without any snRNAseq library
  metadata = metadata[metadata$snRNAseq_ID_count > 0,]
}

unresolvedMssm = metadata[!(metadata$SubID %in% metadataRestrictive$SubID) & (metadata$Brain_bank == "MSSM"),]
unresolvedMssm = unresolvedMssm[!is.na(unresolvedMssm$SNParray_PsychAD),]
unresolvedMssm$problematic_fromSnpArray = ifelse(unresolvedMssm$SNParray_PsychAD %in% problematicSnpArray, T, NA)
main = unresolvedMssm

GENOTYPE_FILE_SEX_MISMATCHES = "/sc/arion/projects/roussp01a/jaro/project_psychAD/metadata/snparray_from_karen_sex_mismatch_samples_list.tsv" # provided by Karen :: List of samples with a mismatch between reported sex and genetically imputed sex from genotype inf
GENOTYPE_FILE_SEX_GTCALLED = "/sc/arion/projects/roussp01a/jaro/project_psychAD/metadata/MA000115_Kleopoulos_PlatesP1-P13_annotated_all_CHR_updated_imputed_sex.sexcheck.csv"
GENOTYPE_FILE_SEX_GTCALLED_FROM_KAREN = "/sc/arion/projects/roussp01a/jaro/project_psychAD/metadata/snparray_sex_from_karen.tsv"
GENOTYPE_FILE_PROBLEMATIC_MISSINGNESS = "/sc/arion/projects/roussp01a/jaro/project_psychAD/metadata/snparray_from_karen_mind_0.05_het_3SD_samples.tsv" # provided by Karen :: List of samples with high (greater than 5%) missingness and/or excess (+/- 3SD) heterozygosity: /sc/arion/projects/roussp01a/karen/NPS_AD_imputed/phenotype_files/mind_0.05_het_3SD_samples.tsv --> I will remove these samples from the final version of QC'd/imputed data once the mismatches are rectified
GENOTYPE_FILE_APOE_COMBINATIONS = "/sc/arion/projects/roussp01a/jaro/project_psychAD/gtcheck/datafiles/APOE_assignment_NPS_AD_microglia_phase2.txt"
GENOTYPE_FILE_ANEUPLOIDY = "/sc/arion/projects/roussp01a/jaro/project_psychAD/gtcheck/datafiles/genotypes_aneuploidy_samples_list.tsv"
GENOTYPE_FILE_ETHNICITY = "/sc/arion/projects/roussp01a/jaro/project_psychAD/gtcheck/datafiles/genotypes_ancestry_assignment_multinom_3_PCs.tsv"

sexGtcalledArray = read.csv(GENOTYPE_FILE_SEX_GTCALLED_FROM_KAREN, sep="\t")
sexGtcalledArray$SubID = sexGtcalledArray$SubID
rownames(sexGtcalledArray) = sexGtcalledArray$SubID
main$sex_fromGenotype = sapply(sexGtcalledArray[main$SNParray_PsychAD,"SNPSEX"], function(sex) {
  ifelse(sex == "F", "Female", ifelse(sex == "M", "Male", "Unassigned"))
})
main$check_sex = main$sex_fromGenotype == main$Sex

sexGtcalledArray = read.csv(GENOTYPE_FILE_SEX_GTCALLED_FROM_KAREN, sep="\t")
sexGtcalledArray$SubID = sexGtcalledArray$SubID
rownames(sexGtcalledArray) = sexGtcalledArray$SubID
main$sex_fromGenotype = sapply(sexGtcalledArray[main$SNParray_PsychAD,"SNPSEX"], function(sex) {
  ifelse(sex == "F", "Female", ifelse(sex == "M", "Male", "Unassigned"))
})
main$check_sex = main$sex_fromGenotype == main$Sex

apoeFromSnpArray = read.csv(GENOTYPE_FILE_APOE_COMBINATIONS, sep="\t")
main$ApoE_gt = as.character(main$apoe_genotype)
main$ApoE_gt_fromSnpArray = apoeFromSnpArray[match(main$SNParray_PsychAD, apoeFromSnpArray$Formatted_ID), "apo_alleles"]
main$ApoE_gt_fromSnpArray = gsub("[e/]", "", main$ApoE_gt_fromSnpArray)
main$MetaMatch_ApoeE_snpArray = sapply(1:nrow(main), function(i) {
  if(is.na(main[i,"ApoE_gt"]) | is.na(main[i,"ApoE_gt_fromSnpArray"])) {
    NA
  } else if(main[i,"ApoE_gt"] == main[i,"ApoE_gt_fromSnpArray"]) {
    T
  } else if(((main[i,"ApoE_gt"] == "13") | (main[i,"ApoE_gt"] == "24")) & (main[i,"ApoE_gt_fromSnpArray"] == "13_or_24")) {
    T
  } else {
    F
  }
})
main$check_apoe = main$ApoE_gt_fromSnpArray == main$apoe_genotype

## SNParray :: Ethnicity
ethnicitySnpArray = read.csv(GENOTYPE_FILE_ETHNICITY, sep="\t")
main$Ethnicity_fromSnpArray = ethnicitySnpArray[match(main$SNParray_PsychAD, ethnicitySnpArray$Formatted_ID), "superpop"]
print(">> SNParray :: Ethnicity")
print(table(main$Ethnicity_fromSnpArray))
main$Ethnicity_standardized = gsub("White", "EUR", gsub("Black", "AFR", gsub("Hispanic", "AMR", main$Ethnicity)))
main$Ethnicity_standardized[main$Ethnicity_standardized %in% c("Other", "Unknown", "Asian", "American Indian or Alaska Native")] = NA

main$MC_Demo_SNParray_Ethnicity = main$Ethnicity_standardized == main$Ethnicity_fromSnpArray

ampad_wgs = read.csv("/sc/arion/projects/roussp01a/jaro/project_psychAD/gtcheck_4//results/king__Ampad_PsychAD_SNParray.txt.kin", sep="\t")
gtcheck_snRNAseq_CommonMind_SNParray = read.csv("/sc/arion/projects/roussp01a/jaro/project_psychAD/gtcheck_4//results/king__CommonMind_SNParray_PsychAD_SNParray.txt.kin", sep="\t")
microglia_snparray = read.csv("/sc/arion/projects/roussp01a/jaro/project_psychAD/gtcheck_4//results/king__Microglia_PsychAD_SNParray.txt.kin", sep="\t")
commonmind_wgs = read.csv("/sc/arion/projects/roussp01a/jaro/project_psychAD/gtcheck_4//results/king__CommonMind_WGS_PsychAD_SNParray.txt.kin", sep="\t")
adsp_wgs = read.csv("/sc/arion/projects/roussp01a/jaro/project_psychAD/gtcheck_4//results/king__PsychAD_SNParray_ADSP_WGS.txt.kin", sep="\t")
rush_wgs = read.csv("/sc/arion/projects/roussp01a/jaro/project_psychAD/gtcheck_4//results/king__Ps", sep="\t")

####

ampad_wgs_copy = ampad_wgs
ampad_wgs_copy$ID1 = ampad_wgs$ID2
ampad_wgs_copy$ID2 = ampad_wgs$ID1
ampad_wgs_merged = rbind(ampad_wgs, ampad_wgs_copy)

metadataExtern = EXTERNAL_METADATA$WGS_Ampad
ampad_wgs_merged = ampad_wgs_merged[(ampad_wgs_merged$ID1 %in% GENOTYPE_METADATA$Formatted_ID) & !(ampad_wgs_merged$ID2 %in% GENOTYPE_METADATA$Formatted_ID),]
dim(ampad_wgs_merged)
ampad_wgs_merged = ampad_wgs_merged[ampad_wgs_merged$ID1 %in% main$SNParray_PsychAD,]
View(ampad_wgs_merged)
write.csv(ampad_wgs_merged[order(ampad_wgs_merged$Kinship, decreasing=T),][1:1000,], file="~/Desktop/ampad.csv")

####

gtcheck_snRNAseq_CommonMind_SNParray_copy = gtcheck_snRNAseq_CommonMind_SNParray
gtcheck_snRNAseq_CommonMind_SNParray_copy$ID1 = gtcheck_snRNAseq_CommonMind_SNParray$ID2
gtcheck_snRNAseq_CommonMind_SNParray_copy$ID2 = gtcheck_snRNAseq_CommonMind_SNParray$ID1
gtcheck_snRNAseq_CommonMind_SNParray_merged = rbind(gtcheck_snRNAseq_CommonMind_SNParray, gtcheck_snRNAseq_CommonMind_SNParray_copy)

metadataExtern = EXTERNAL_METADATA$SNParray_CommonMind
gtcheck_snRNAseq_CommonMind_SNParray_merged = gtcheck_snRNAseq_CommonMind_SNParray_merged[(gtcheck_snRNAseq_CommonMind_SNParray_merged$ID1 %in% GENOTYPE_METADATA$Formatted_ID) & !(gtcheck_snRNAseq_CommonMind_SNParray_merged$ID2 %in% GENOTYPE_METADATA$Formatted_ID),]
gtcheck_snRNAseq_CommonMind_SNParray_merged = gtcheck_snRNAseq_CommonMind_SNParray_merged[gtcheck_snRNAseq_CommonMind_SNParray_merged$ID1 %in% main$SNParray_PsychAD,]
gtcheck_snRNAseq_CommonMind_SNParray_merged$ID2_CMC_ID = metadataExtern[match(gtcheck_snRNAseq_CommonMind_SNParray_merged$ID2, metadataExtern$SNP_report.Genotyping_Sample_ID),"Individual_ID"]
dim(gtcheck_snRNAseq_CommonMind_SNParray_merged)
write.csv(gtcheck_snRNAseq_CommonMind_SNParray_merged[order(gtcheck_snRNAseq_CommonMind_SNParray_merged$Kinship, decreasing=T),][1:1000,], file="~/Desktop/commonmind.csv")

####

gtcheck_snRNAseq_CommonMind_SNParray_copy = commonmind_wgs
gtcheck_snRNAseq_CommonMind_SNParray_copy$ID1 = commonmind_wgs$ID2
gtcheck_snRNAseq_CommonMind_SNParray_copy$ID2 = commonmind_wgs$ID1
gtcheck_snRNAseq_CommonMind_SNParray_merged = rbind(gtcheck_snRNAseq_CommonMind_SNParray, gtcheck_snRNAseq_CommonMind_SNParray_copy)

metadataExtern = EXTERNAL_METADATA$WGS_CommonMind
gtcheck_snRNAseq_CommonMind_SNParray_merged = gtcheck_snRNAseq_CommonMind_SNParray_merged[(gtcheck_snRNAseq_CommonMind_SNParray_merged$ID1 %in% GENOTYPE_METADATA$Formatted_ID) & !(gtcheck_snRNAseq_CommonMind_SNParray_merged$ID2 %in% GENOTYPE_METADATA$Formatted_ID),]
gtcheck_snRNAseq_CommonMind_SNParray_merged = gtcheck_snRNAseq_CommonMind_SNParray_merged[gtcheck_snRNAseq_CommonMind_SNParray_merged$ID1 %in% main$SNParray_PsychAD,]
gtcheck_snRNAseq_CommonMind_SNParray_merged$ID2_CMC_ID = metadataExtern[match(gtcheck_snRNAseq_CommonMind_SNParray_merged$ID2, metadataExtern$SNP_report.Genotyping_Sample_ID),"Individual_ID"]
dim(gtcheck_snRNAseq_CommonMind_SNParray_merged)
write.csv(gtcheck_snRNAseq_CommonMind_SNParray_merged[order(gtcheck_snRNAseq_CommonMind_SNParray_merged$Kinship, decreasing=T),][1:1000,], file="~/Desktop/commonmind_wgs.csv")

####

adsp_wgs_copy = adsp_wgs
adsp_wgs_copy$ID1 = adsp_wgs$ID2
adsp_wgs_copy$ID2 = adsp_wgs$ID1
adsp_wgs_merged = rbind(adsp_wgs, adsp_wgs_copy)

metadataExtern = EXTERNAL_METADATA$ADSP_WGS
adsp_wgs_merged = adsp_wgs_merged[(adsp_wgs_merged$ID1 %in% GENOTYPE_METADATA$Formatted_ID) & !(adsp_wgs_merged$ID2 %in% GENOTYPE_METADATA$Formatted_ID),]
dim(adsp_wgs_merged)
adsp_wgs_merged = adsp_wgs_merged[adsp_wgs_merged$ID1 %in% main$SNParray_PsychAD,]
View(adsp_wgs_merged)
write.csv(adsp_wgs_merged[order(adsp_wgs_merged$Kinship, decreasing=T),][1:1000,], file="~/Desktop/adsp.csv")

####

mglia = microglia_snparray
mglia_copy = microglia_snparray
mglia_copy$ID1 = mglia$ID2
mglia_copy$ID2 = mglia$ID1
mglia_merged = rbind(mglia, mglia_copy)

metadataExtern = EXTERNAL_METADATA$mglia
mglia_merged = mglia_merged[(mglia_merged$ID1 %in% GENOTYPE_METADATA$Formatted_ID) & !(mglia_merged$ID2 %in% GENOTYPE_METADATA$Formatted_ID),]
dim(mglia_merged)
mglia_merged = mglia_merged[mglia_merged$ID1 %in% main$SNParray_PsychAD,]
View(mglia_merged)
write.csv(mglia_merged[order(mglia_merged$Kinship, decreasing=T),][1:1000,], file="~/Desktop/mglia.csv")

####

rush_wgs_copy = rush_wgs
rush_wgs_copy$ID1 = rush_wgs$ID2
rush_wgs_copy$ID2 = rush_wgs$ID1
rush_wgs_merged = rbind(rush_wgs, rush_wgs_copy)

metadataExtern = EXTERNAL_METADATA$WGS_RUSH
rush_wgs_merged = rush_wgs_merged[(rush_wgs_merged$ID1 %in% GENOTYPE_METADATA$Formatted_ID) & !(rush_wgs_merged$ID2 %in% GENOTYPE_METADATA$Formatted_ID),]
dim(rush_wgs_merged)
rush_wgs_merged = rush_wgs_merged[rush_wgs_merged$ID1 %in% main$SNParray_PsychAD,]
View(rush_wgs_merged)
write.csv(rush_wgs_merged[order(rush_wgs_merged$Kinship, decreasing=T),][1:1000,], file="~/Desktop/mglia.csv")


########################################################################################################
# Add metadata from PRS
########################################################################################################

{
  # Add PRS score to metadata df
  prsDeepika = read.csv(PRS_FROM_DEEPIKA, sep="\t")
  colnames(prsDeepika) = paste0("prs_", colnames(prsDeepika))
  metadata = cbind.data.frame(metadata, prsDeepika[match(metadata$SNParray_PsychAD, prsDeepika$prs_sample_id),2:ncol(prsDeepika)])
  
  selcols = c("BRAAK", "CERAD", "CDRScore", "PLAQUE", paste0("nps_", c(names(MSSM_NPS), gsub("Cur", "Hx", names(MSSM_NPS))) ))
  
  prsCorrelation = do.call("rbind.data.frame", lapply(selcols, function(pheno) {
    sapply(colnames(prsDeepika)[2:ncol(prsDeepika)], function(x) {
      i = !is.na(metadata[,x]) & !is.na(metadata[,pheno])
      phenoVector = metadata[i,pheno]
      if(is.logical(metadata[i,pheno])) {
        phenoVector = ifelse(metadata[i,pheno], 1, 0)
      }
      cor.test(metadata[i,x], phenoVector, method="pearson")$estimate
    })
  }))
  rownames(prsCorrelation) = selcols
  colnames(prsCorrelation) = gsub("prs_", "", colnames(prsDeepika)[2:ncol(prsDeepika)])

  prsCorrelation$phenotype = rownames(prsCorrelation)
  
  prsCorrelationDf = reshape2::melt(prsCorrelation)
  prsCorrelationDf$phenotype = gsub("nps_", "", prsCorrelationDf$phenotype)
  prsCorrelationDf$phenotype = ordered(prsCorrelationDf$phenotype, levels=c(unique(prsCorrelationDf$phenotype)))
  colnames(prsCorrelationDf) = gsub("value", "Spearman", colnames(prsCorrelationDf))
  
  prsCorrelation = ggplot(prsCorrelationDf, aes(variable, phenotype, fill=Spearman)) + 
    geom_tile() + scale_fill_gradient2(low="navy", mid="white", high="red", midpoint=0) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("GWAS") + ylab("Phenotype")
  mpdf("prsCorrelation", width=15, height=6); print(prsCorrelation); dev.off()
  prsCorrelation
  htmlwidgets::saveWidget(as_widget(ggplotly(prsCorrelation, width=2000, height=800)), paste0(outDir, "prsCorrelation.html"))
}

########################################################################################################
# Add CDR domains for MSSM samples
########################################################################################################

{
  # Load & add CDR subdomains to the main metadata
  cdrDomains = read.csv(CDR_DOMAINS)
  colnames(cdrDomains)[2:ncol(cdrDomains)] = paste0("CDR_", gsub("Home.Hobies", "HomeHobbies", colnames(cdrDomains)[2:ncol(cdrDomains)] ))
  metadata = cbind.data.frame(metadata, cdrDomains[match(metadata$BB, cdrDomains$BB), 2:ncol(cdrDomains)])
  
  # Plot correlations among CDR subdomains
  cdrDomainsFull = cbind.data.frame(cdrDomains, metadata[match(cdrDomains$BB, metadata$BB),"CDRScore"])
  colnames(cdrDomainsFull) = c(colnames(cdrDomains), "CDRScore")
  cdrDomainsFull = cdrDomainsFull[!is.na(cdrDomainsFull$CDRScore),]
  corrMatrix = cdrDomainsFull[,2:ncol(cdrDomainsFull)]
  corrMatrixOut = data.frame(sapply(1:ncol(corrMatrix), function(x) sapply(1:ncol(corrMatrix), function(y) { cor.test(as.numeric(corrMatrix[,x]), as.numeric(corrMatrix[,y]), method="spearman")$estimate } )))
  colnames(corrMatrixOut) = colnames(cdrDomainsFull[,2:ncol(cdrDomainsFull)])
  rownames(corrMatrixOut) = colnames(corrMatrixOut)
  corrMatrixOut = get_upper_tri(corrMatrixOut)
  corrMatrixOut$ID = rownames(corrMatrixOut)
  corrMatrixOut = melt(corrMatrixOut, id.vars="ID", factorsAsStrings=T)
  corrMatrixOut$ID = ordered(corrMatrixOut$ID, levels=colnames(cdrDomainsFull[,2:ncol(cdrDomainsFull)]))
  
  corr = ggplot(corrMatrixOut, aes(variable, ID, fill = value)) +
    geom_tile() + scale_fill_material("red") +
    coord_equal() + theme_bw() + 
    theme(axis.title.x=element_blank(), axis.text.x=element_text(colour = "black", angle = 45, hjust = 1), axis.text.y=element_text(colour = "black"), axis.ticks.x=element_blank(), 
          axis.title.y=element_blank(),  axis.ticks.y=element_blank())
  corr = corr + geom_text(aes(variable, ID, label = round(value, 2)), color = "black", size = 3)
  mpdf("cdrCorrelation", width=15, height=6); print(corr); dev.off()
}

########################################################################################################
# Add links to genotypes & metadata derived from genotypes
########################################################################################################

{
  metadata$Sex_chr_aneuploidy = sapply(metadata$SNParray_PsychAD, function(snparray_id) { 
    if(is.na(snparray_id)) { NA } else if(snparray_id %in% SEX_CHR_ANEUPLOIDY) { T } else { F }
  })
  
  for(dset in dsetForComparison) {
    samples = read.table(paste0(dset$plink, ".fam"))[,2]
    external_metadata = read.csv(dset$metadata)
    
    metadata[,dset$name] = sapply(metadata[,dset$convertor_subjectMatchColumn], function(x) {
      matchedVal = NA
      if(!is.na(x)) {
        matchedVal = external_metadata[match(x, external_metadata[,dset$metadata_subjectMatchColumn])[1], dset$metadata_sampleMatchColumn]
        matchedVal = ifelse(matchedVal %in% samples, matchedVal, NA)
      }
      if((dset$name == "WGS_RUSH") & !is.na(matchedVal)) {
        matchedVal = ifelse(grepl(matchedVal, pattern="_"), matchedVal, paste0(matchedVal, "_", matchedVal))
      }
      matchedVal
    })
    print(paste0("> Mapped for ", dset$name, ": ", sum(!is.na(metadata[,dset$name]))))
  }
  
  gtCols = c("SNParray_HBBC", "SNParray_CommonMind", "WGS_Ampad", "SNParray_Microglia", "WGS_RUSH", "SNParray_PsychAD", "ADSP_WGS") # "WGS_CommonMind"
  gtAvailability = lapply(gtCols, function(dsetName) {
    metadata[!is.na(metadata[,dsetName]),"SubID"]
  })
  names(gtAvailability) = gtCols
  mpdf("gtAvailability", width=15, height=7); print(upset(fromList(gtAvailability), order.by = "freq", nsets = length(gtAvailability))); dev.off()
  print(upset(fromList(gtAvailability), order.by = "freq", nsets = length(gtAvailability)))
  
  gtCols = c("SNParray_PsychAD", "SNParray_HBBC", "SNParray_CommonMind", "ADSP_SampleId", "WGS_RUSH", "WGS_CommonMind", "WGS_Ampad", "SNParray_Microglia")
  gtColsDef = c("SNParray_PsychAD", "SNParray_HBBC", "SNParray_CommonMind","ADSP_SampleId", "WGS_RUSH", "WGS_CommonMind", "WGS_Ampad", "SNParray_Microglia")
  
  metadata$primary_genotype = sapply(1:nrow(metadata), function(i) {
    dset = gtCols[which(!is.na(metadata[i,gtCols]))[1]]
    ifelse(is.na(dset), NA, metadata[i,dset])
  })
  metadata$primary_genotype_dset = sapply(1:nrow(metadata), function(i) {
    gtCols[which(!is.na(metadata[i,gtCols]))[1]]
  })
  sum(is.na(metadata$primary_genotype))
}

########################################################################################################
# Save donor-level metadata to Synapse & Minerva
########################################################################################################

{
  # Create list of variables that will be part of main metadata clinical file
  prioritized = c("SubID", "Brain_bank", "Age", "Sex", "Sex_chr_aneuploidy", "Ethnicity", "Dx", "pH", "PMI", "CERAD", "BRAAK_AD", "BRAAK_PD", "CDRScore", "Plq_Mn", "ApoE_gt", "CMC_individual_ID", "AMPAD_msbb_individualID", "snRNAseq_ID", unlist(DIAGNOSIS_CLASSIFICATION, use.names=F))
  prioritized = c(prioritized, "SNParray_HBBC", "SNParray_CommonMind", "WGS_CommonMind", "WGS_RUSH", "WGS_Ampad", "SNParray_Microglia", "SNParray_PsychAD", "ADSP_SampleId")
  prioritized = c(prioritized, paste0("nps_", names(MSSM_NPS)),  paste0("nps_", gsub("Cur", "Hx", names(MSSM_NPS))))
  prioritized = c(prioritized, colnames(prsDeepika)[2:length(colnames(prsDeepika))])
  prioritized = c(prioritized, as.character(unlist(ADDIT_NEUROPAT_METRICS)))
  prioritized = c(prioritized, colnames(cdrDomains)[2:ncol(cdrDomains)])
  
  # Save metadata to Minerva
  colnames(metadata) = gsub("^BRAAK$", "BRAAK_AD", colnames(metadata))
  colnames(metadata) = gsub("^apoe_genotype$", "ApoE_gt", colnames(metadata))
  write.csv(metadata[,prioritized], file = file.path(outDir, "clinical_metadata.csv"), row.names=F)
  write.csv(metadata, file = "~/Desktop/clinical_metadata_full.csv", row.names=F)
  write.csv(metadata[,prioritized], file = FINAL_METADATA_PATH, row.names=F)
  
  # Save metadata to Synapse
  file = synStore(File(path="~/Desktop/clinical_metadata.csv", parent=SYNAPSE$root), used="https://github.com/DiseaseNeuroGenomics/psychad_metadata/blob/main/psychad_clinical_metadata_merge.R")
  file = synStore(File(path="~/Desktop/clinical_metadata_full.csv", parent=SYNAPSE$individual_files), used="https://github.com/DiseaseNeuroGenomics/psychad_metadata/blob/main/psychad_clinical_metadata_merge.R")
  
  # Copy merge & querying code for metadata to public github
  cmd_copy_metadata = paste0("cp ", unlist(PSYCHAD_METADATA_CODES), " ", PSYCHAD_METADATA_GITHUB_MINERVA)
  #system(paste0(cmd_copy_metadata, collapse=";"))
  #print(paste0("cd ", PSYCHAD_METADATA_GITHUB_MINERVA, "; git commit -a; git push"))

  # Paths to the code from Jaro's private bitbucket. We will copy code from merging & querying clinical metadata from there to his clone of psychad-metadata public github  
  PSYCHAD_METADATA_CODES = list(
    metadata_create = "/sc/arion/projects/roussp01a/jaro/atacseq_ad/NYGC_AD_R01/other_projects/psychad_clinical_metadata_merge.R",
    metadata_querying = "/sc/arion/projects/roussp01a/jaro/atacseq_ad/NYGC_AD_R01/other_projects/psychad_clinical_metadata_querying.R"
  )
  
  # Save genotype convertor
  write.csv(metadata[,c("SubID", "SNParray_HBBC", "SNParray_CommonMind", "WGS_CommonMind", "WGS_RUSH", "WGS_Ampad", "SNParray_Microglia", "SNParray_PsychAD", "ADSP_SampleId")], file=file.path(outDir, "NPSAD_gt_converter_Jaro.csv"), quote=F, row.names=F)
}

########################################################################################################
# Save sample-level metadata to Synapse & Minerva
########################################################################################################

{
  allInfo$ID_very_short = sapply(allInfo$ID_short, function(x) {
    idx = strsplit(x, split = "_")[[1]]
    paste0(idx[1:length(idx)-1], collapse="_")
  })
    
  demuxStats = read.csv(DEMUX_STATS, sep="\t")
  demuxStats = demuxStats[nchar(demuxStats$donors) > 0,]
  demuxStatsDf = do.call("rbind.data.frame", lapply(1:nrow(demuxStats), function(i) {
    donors = strsplit(demuxStats[i, "SubIDs"], split=",")[[1]]
    cells = strsplit(demuxStats[i, "cell_counts"], split=",")[[1]]
    df = cbind.data.frame(paste0(donors, "_", gsub("\\-cDNA", "", demuxStats[i,"sample"])), cells)
    colnames(df) = c("sample_id", "cells")
    df
  }))
  rownames(demuxStatsDf) = demuxStatsDf$sample_id
  
  prioritized_allInfo = c("ID_short", "SubID", "poolName", "chr1reads")
  allInfoWcells = cbind.data.frame(allInfo[,prioritized_allInfo], as.numeric(demuxStatsDf[allInfo$ID_very_short, "cells"]))
  colnames(allInfoWcells) = c("snRNAseq_ID", "SubID", "poolName", "chr1_reads", "cell_count")
  write.csv(allInfoWcells, file=FINAL_ALLINFO_PATH, quote=F, row.names=F)
}

########################################################################################################
# Final plots
########################################################################################################

{
  all_dx = setdiff(unlist(DIAGNOSIS_CLASSIFICATION), DIAGNOSIS_BENIGN)
  all_dx = all_dx[sapply(all_dx, function(x) sum(na.omit(metadata[,x])==1) > 20)] # Keep only categories with non-zero counts
  deleterious_dx = setdiff(all_dx, DIAGNOSIS_BENIGN)
  metadata$Control = sapply(1:nrow(metadata), function(i) 
    ifelse((sum(na.omit(as.integer(metadata[i,all_dx]))) > 0) | (metadata[i,"CERAD"] %in% c(2,3,4)) | (metadata[i,"BRAAK"] %in% c(1,2,3,4,5,6)), 0, 1))
  
  listInput = lapply(c(setdiff(all_dx, c(DIAGNOSIS_BENIGN, "MCI", "Dementia")), "Control"), function(colName) {
    metadata[which(as.integer(metadata[,colName])==1), "SubID"]
  })
  names(listInput) = c(setdiff(all_dx, c(DIAGNOSIS_BENIGN, "MCI", "Dementia")), "Control")
  mpdf("diagnosis_stratification", width=15, height=6); print(upset(fromList(listInput), order.by = "freq", nsets = length(all_dx))); dev.off()
}

{
  cellCountDf = do.call("rbind.data.frame", sapply(1:nrow(demuxStats), function(i) { 
    x = as.numeric(strsplit(demuxStats$cell_counts[i], split=",")[[1]])
    x = x[order(x, decreasing=T)]
  }))
  rownames(cellCountDf) = demuxStats$sample
  colnames(cellCountDf) = paste0("libcells_", 1:6)
  
  cellCountDf$suma = sapply(1:nrow(cellCountDf), function(i) { sum(cellCountDf[i,1:6]) })
  for(lib in paste0("libcells_", 1:6)) {
    cellCountDf[lib] = cellCountDf[,lib] / cellCountDf$suma
  }
  cellCountDf = cellCountDf[,!colnames(cellCountDf) %in% c("suma")]
  
  cellCountDf$poolName = demuxStats$sample
  cellCountDf$poolName = ordered(cellCountDf$poolName, levels=cellCountDf$poolName[order(cellCountDf$libcells_1, decreasing=T)])
  
  cellCountMeltDf = reshape2::melt(cellCountDf)
  cellCountMeltDf$variable = ordered(cellCountMeltDf$variable, levels=paste0("libcells_", 1:6))
  
  myPlot = ggplot(cellCountMeltDf, aes(fill=variable, y=value, x=poolName)) + geom_bar(position="fill", stat="identity") + xlab("Pool") + ylab("Proportion of cells") +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
  mpdf(paste0("sampleCount_vireo_distribution"), width=8, height=4); print(myPlot); dev.off()
  
  myPlot = ggplot(cellCountMeltDf, aes(fill=variable, y=value, x=variable)) + geom_boxplot() + xlab("Sample in pool") + ylab("Proportion of cells") +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 0.5,
                       axis.text.y=element_text(colour = "black"), axis.ticks.x=element_blank(), legend.position = "none") + ylim(c(0,1))
  mpdf(paste0("sampleCount_vireo_boxplot"), width=8, height=4); print(myPlot); dev.off()
}

{
  DEMUX_STATS_HTO = "/sc/arion/projects/roussp01a/jaro/project_psychAD/metadata/snrnaseq_qc_from_prashant_old.csv"
  demuxHtoStats = read.csv(DEMUX_STATS_HTO, skip=2)[,c("SAMPLE", "N_CELLS_AFTER_DEMUX_CS")]
  demuxHtoStats = demuxHtoStats[!duplicated(demuxHtoStats$SAMPLE),]
  
  cellCountDf = do.call("rbind.data.frame", sapply(1:nrow(demuxHtoStats), function(i) { 
    x = sapply(strsplit(demuxHtoStats$N_CELLS_AFTER_DEMUX_CS[i], split=",")[[1]], function(x) as.numeric(strsplit(x, ":")[[1]][2]) )
    x
  }))
  rownames(cellCountDf) = demuxHtoStats$sample
  colnames(cellCountDf) = paste0("libcells_", 1:6)
  cellCountDf = cellCountDf[,1:6]
  
  cellCountDf$suma = sapply(1:nrow(cellCountDf), function(i) { sum(cellCountDf[i,1:6]) })
  for(lib in paste0("libcells_", 1:6)) {
    cellCountDf[lib] = cellCountDf[,lib] / cellCountDf$suma
  }
  cellCountDf = cellCountDf[,!colnames(cellCountDf) %in% c("suma")]
  
  cellCountDf$poolName = demuxHtoStats$SAMPLE
  cellCountDf$poolName = ordered(cellCountDf$poolName, levels=cellCountDf$poolName[order(cellCountDf$libcells_1, decreasing=T)])
  
  cellCountMeltDf = reshape2::melt(cellCountDf)
  cellCountMeltDf$variable = ordered(cellCountMeltDf$variable, levels=paste0("libcells_", 1:6))
  
  myPlot = ggplot(cellCountMeltDf, aes(fill=variable, y=value, x=poolName)) + geom_bar(position="fill", stat="identity") + xlab("Pool") + ylab("Proportion of cells") +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
  mpdf(paste0("sampleCount_hto_distribution"), width=8, height=4); print(myPlot); dev.off()
  
  myPlot = ggplot(cellCountMeltDf, aes(fill=variable, y=value, x=variable)) + geom_boxplot() + xlab("Sample in pool") + ylab("Proportion of cells") +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 0.5,
                       axis.text.y=element_text(colour = "black"), axis.ticks.x=element_blank(), legend.position = "none") + ylim(c(0,1))
  mpdf(paste0("sampleCount_hto_boxplot"), width=8, height=4); print(myPlot); dev.off()
}







