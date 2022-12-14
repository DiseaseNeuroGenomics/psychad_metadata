library(UpSetR)
library(reshape)
library(ggplot2)
library(ggsci)
library(synapser)

# Load data from Synapse
synLogin('','')
metadataAll = read.csv(synGet(entity='syn44716176')$path, header = T)
metadata = metadataAll[metadataAll$snRNAseq_ID_count > 0,]
rownames(metadata) = metadata$SubID

PSYCHAD_METADATA_CODES = list(
  metadata_create = "/sc/arion/projects/roussp01a/jaro/atacseq_ad/NYGC_AD_R01/other_projects/psychad_clinical_metadata_merge.R",
  metadata_querying = "/sc/arion/projects/roussp01a/jaro/atacseq_ad/NYGC_AD_R01/other_projects/psychad_clinical_metadata_querying.R"
)
PSYCHAD_METADATA_GITHUB_MINERVA = "/sc/arion/projects/roussp01a/jaro/repositories/psychad_metadata/"     # Minerva clone of github code for creating (this file) & querying clinical metadata (we will copy scripts there)

##########################
# Sets of diagnosis levels
DIAGNOSIS_CLASSIFICATION = list(
  "neurodegenerative" = c("AD", "MCI", "Dementia", "PD", "PD_uncertain_plus_encephalitic", "DLBD", "FTD", "ALS", "Others_Neurodegenerative"),
  "neurological" = c("MS", "PSP", "Epilepsy", "Seizures", "Tumor", "Migraine_headaches", "Head_Injury", "Vascular", "Others_Neurological"), 
  "neuropsychiatric" = c("SCZ", "MDD", "BD_unspecific", "BD_I", "BD_II", "PTSD", "ADHD", "OCD", "Tardive_Dyskinesia_Neuroleptic_induced", "Schizoaffective_bipolar", "Schizoaffective_depressive", 
                         "Anorexia", "Bulimia", "Anxiety", "Binge_Purge", "Eating_disorder", "Others_Neuropsychiatric"),
  "metabolic" = c("Diabetes_mellitus_unspecified", "TD_I", "TD_II"))
DIAGNOSIS_BENIGN = c("Anxiety", "Migraine_headaches", DIAGNOSIS_CLASSIFICATION$metabolic)

##########################
# UpsetR showing overlaps across diagnosis groups (here it is for *clinical* diagnosis; pathology-level metrics not taken into account)
#metadata=metadata[metadata$PD==1,]
all_dx = setdiff(unlist(DIAGNOSIS_CLASSIFICATION), DIAGNOSIS_BENIGN)
#all_dx = all_dx[sapply(all_dx, function(x) sum(na.omit(metadata[,x])==1) > 0)] # Keep only categories with non-zero counts
deleterious_dx = setdiff(all_dx, DIAGNOSIS_BENIGN)
#metadata$Control = sapply(1:nrow(metadata), function(i) 
#  ifelse((sum(na.omit(as.integer(metadata[i,all_dx]))) > 0) | (metadata[i,"CERAD"] %in% c(2,3,4)) | (metadata[i,"BRAAK_AD"] %in% c(1,2,3,4,5,6)), 0, 1))
#metadata$Control = sapply(1:nrow(metadata), function(i) 
#  ifelse((sum(na.omit(as.integer(metadata[i,all_dx]))) > 0) | (metadata[i,"CERAD"] %in% c(2,3,4)) | (metadata[i,"BRAAK_AD"] %in% c(1,2,3,4,5,6)), 0, 1))
#metadata$Control = sapply(1:nrow(metadata), function(i) sum(na.omit(as.integer(metadata[i, all_dx]))) == 0 )

listInput = lapply(c(setdiff(all_dx, c(DIAGNOSIS_BENIGN, "MCI", "Dementia"))), function(colName) {
  print(colName)
  metadata[which(as.integer(metadata[,colName]) == 1), "SubID"]
})
names(listInput) = c(setdiff(all_dx, c(DIAGNOSIS_BENIGN, "MCI", "Dementia")))
upset(fromList(listInput), order.by = "freq", nsets = length(all_dx))

#hist(metadata[metadata$Control==1, "Age"], breaks = 20)
#hist(metadata[(metadata$Control==1) & (metadata$Brain_bank=="HBCC"), "Age"], breaks = 20)

##########################
# List of sample names per diagnosis groups

metadata[,unlist(DIAGNOSIS_CLASSIFICATION)][is.na(metadata[,unlist(DIAGNOSIS_CLASSIFICATION)])] = 0
metadata$SCZ_ALL = sapply(1:nrow(metadata), function(i) ifelse(sum(na.omit(metadata[i, c("SCZ", "Schizoaffective_bipolar", "Schizoaffective_depressive")])) > 0, 1, 0))
metadata$BD_ALL = sapply(1:nrow(metadata), function(i) ifelse(sum(na.omit(metadata[i, c("BD_unspecific", "BD_I", "BD_II")])) > 0, 1, 0))
metadata$PD_ALL = sapply(1:nrow(metadata), function(i) ifelse(sum(na.omit(metadata[i, c("PD", "PD_uncertain_plus_encephalitic")])) > 0, 1, 0))
metadata$TD_ALL = sapply(1:nrow(metadata), function(i) ifelse(sum(na.omit(metadata[i, c("Diabetes_mellitus_unspecified", "TD_I", "TD_II")])) > 0, 1, 0))

SubID_groups = list(
  "controls_neuropathological" = na.omit(sapply(1:nrow(metadata), function(i) {
    ifelse(((metadata[i,"CERAD"] %in% c(1, 2)) & (metadata[i,"BRAAK_AD"] %in% c(0, 1, 2))) | ((metadata[i,"Brain_bank"] == "HBCC")), 
           metadata[i,"SubID"], NA)
  })),
  "controls_neuropathological_clinical" = na.omit(sapply(1:nrow(metadata), function(i) {
    ifelse((sum(na.omit(as.integer(metadata[i,setdiff(deleterious_dx, c("MCI", "Dementia", "AD"))]))) == 0) & 
             (((metadata[i,"CERAD"] %in% c(1, 2)) & (metadata[i,"BRAAK_AD"] %in% c(0, 1, 2))) | ((metadata[i,"Brain_bank"] == "HBCC"))), 
           metadata[i,"SubID"], NA)
  })),
  "controls_supercontrols" = na.omit(sapply(1:nrow(metadata), function(i) {
    ifelse((sum(na.omit(as.integer(metadata[i,setdiff(deleterious_dx, c("AD"))]))) == 0) & 
             (((metadata[i,"CERAD"] %in% c(1)) & (metadata[i,"BRAAK_AD"] %in% c(0, 1, 2))) | ((metadata[i,"Brain_bank"] == "HBCC"))), 
           metadata[i,"SubID"], NA)
  })),
  "controls_supercontrols_for_aging" = na.omit(sapply(1:nrow(metadata), function(i) {
    ifelse((sum(na.omit(as.integer(metadata[i,setdiff(deleterious_dx, c("AD"))]))) == 0) & 
             (((metadata[i,"CERAD"] %in% c(1)) & (metadata[i,"BRAAK_AD"] %in% c(0, 1, 2))) | ((metadata[i,"Brain_bank"] == "HBCC"))) &
             (metadata[i,"Brain_bank"] != "RUSH"), 
           metadata[i,"SubID"], NA)
  })),
  "AD" = na.omit(sapply(1:nrow(metadata), function(i) {
    ifelse((sum(na.omit(as.integer(metadata[i,setdiff(deleterious_dx, c("MCI", "Dementia", "AD"))]))) == 0) & 
             (metadata[i,"AD"] == 1), 
           metadata[i,"SubID"], NA)
  })),
  "AD_strict" = na.omit(sapply(1:nrow(metadata), function(i) {
    ifelse((sum(na.omit(as.integer(metadata[i,setdiff(deleterious_dx, c("MCI", "Dementia", "AD"))]))) == 0) & 
             ((metadata[i,"CERAD"] %in% c(4)) & (metadata[i,"BRAAK_AD"] %in% c(3, 4, 5, 6)) & (metadata[i,"Dementia"]) & (metadata[i,"Brain_bank"] != "HBCC")), 
           metadata[i,"SubID"], NA)
  })),
  "CERAD_1" = na.omit(sapply(1:nrow(metadata), function(i) {
    ifelse((metadata[i,"CERAD"] == 1) & (sum(na.omit(as.integer(metadata[i,setdiff(deleterious_dx,c("AD", "MCI", "Dementia"))]))) == 0),
           metadata[i,"SubID"], NA)
  })),
  "CERAD_2" = na.omit(sapply(1:nrow(metadata), function(i) {
    ifelse((metadata[i,"CERAD"] == 2) & (sum(na.omit(as.integer(metadata[i,setdiff(deleterious_dx,c("AD", "MCI", "Dementia"))]))) == 0),
           metadata[i,"SubID"], NA)
  })),
  "CERAD_3" = na.omit(sapply(1:nrow(metadata), function(i) {
    ifelse((metadata[i,"CERAD"] == 3) & (sum(na.omit(as.integer(metadata[i,setdiff(deleterious_dx,c("AD", "MCI", "Dementia"))]))) == 0),
           metadata[i,"SubID"], NA)
  })),
  "CERAD_4" = na.omit(sapply(1:nrow(metadata), function(i) {
    ifelse((metadata[i,"CERAD"] == 4) & (sum(na.omit(as.integer(metadata[i,setdiff(deleterious_dx,c("AD", "MCI", "Dementia"))]))) == 0),
           metadata[i,"SubID"], NA)
  })),
  "Braak_0" = na.omit(sapply(1:nrow(metadata), function(i) {
    ifelse((metadata[i,"BRAAK_AD"] == 0) & (sum(na.omit(as.integer(metadata[i,setdiff(deleterious_dx,c("AD", "MCI", "Dementia"))]))) == 0),
           metadata[i,"SubID"], NA)
  })),
  "Braak_1" = na.omit(sapply(1:nrow(metadata), function(i) {
    ifelse((metadata[i,"BRAAK_AD"] == 1) & (sum(na.omit(as.integer(metadata[i,setdiff(deleterious_dx,c("AD", "MCI", "Dementia"))]))) == 0),
           metadata[i,"SubID"], NA)
  })),
  "Braak_2" = na.omit(sapply(1:nrow(metadata), function(i) {
    ifelse((metadata[i,"BRAAK_AD"] == 2) & (sum(na.omit(as.integer(metadata[i,setdiff(deleterious_dx,c("AD", "MCI", "Dementia"))]))) == 0),
           metadata[i,"SubID"], NA)
  })),
  "Braak_3" = na.omit(sapply(1:nrow(metadata), function(i) {
    ifelse((metadata[i,"BRAAK_AD"] == 3) & (sum(na.omit(as.integer(metadata[i,setdiff(deleterious_dx,c("AD", "MCI", "Dementia"))]))) == 0),
           metadata[i,"SubID"], NA)
  })),
  "Braak_4" = na.omit(sapply(1:nrow(metadata), function(i) {
    ifelse((metadata[i,"BRAAK_AD"] == 4) & (sum(na.omit(as.integer(metadata[i,setdiff(deleterious_dx,c("AD", "MCI", "Dementia"))]))) == 0),
           metadata[i,"SubID"], NA)
  })),
  "Braak_5" = na.omit(sapply(1:nrow(metadata), function(i) {
    ifelse((metadata[i,"BRAAK_AD"] == 5) & (sum(na.omit(as.integer(metadata[i,setdiff(deleterious_dx,c("AD", "MCI", "Dementia"))]))) == 0),
           metadata[i,"SubID"], NA)
  })),
  "Braak_6" = na.omit(sapply(1:nrow(metadata), function(i) {
    ifelse((metadata[i,"BRAAK_AD"] == 6) & (sum(na.omit(as.integer(metadata[i,setdiff(deleterious_dx,c("AD", "MCI", "Dementia"))]))) == 0),
           metadata[i,"SubID"], NA)
  })),
  "AD_DLBD" = na.omit(sapply(1:nrow(metadata), function(i) {
    ifelse((sum(na.omit(metadata[i,setdiff(deleterious_dx,c("AD", "MCI", "Dementia", "DLBD"))])) == 0) &
             (metadata[i,"AD"] == 1) & (metadata[i,"DLBD"] == 1), metadata[i,"SubID"], NA)
  })),
  "AD_Vascular" = na.omit(sapply(1:nrow(metadata), function(i) {
    ifelse((sum(na.omit(metadata[i,setdiff(deleterious_dx,c("AD", "MCI", "Dementia", "DLBD", "Vascular"))])) == 0) &
             (metadata[i,"AD"] == 1) & (metadata[i,"Vascular"] == 1), metadata[i,"SubID"], NA)
  })),
  "AD_SCZ" = na.omit(sapply(1:nrow(metadata), function(i) {
    ifelse((sum(na.omit(metadata[i,setdiff(deleterious_dx,c("AD", "MCI", "Dementia", "SCZ_ALL", "SCZ", "Schizoaffective_bipolar", "Schizoaffective_depressive"))])) == 0) &
             (metadata[i,"AD"] == 1) & (metadata[i,"SCZ_ALL"] == 1), metadata[i,"SubID"], NA)
  })),
  "SCZ" = na.omit(sapply(1:nrow(metadata), function(i) {
    ifelse((sum(na.omit(metadata[i,setdiff(deleterious_dx,c("SCZ_ALL", "SCZ", "Schizoaffective_bipolar", "Schizoaffective_depressive", "Dementia", "MCI"))])) == 0) & (metadata[i,"SCZ_ALL"] == 1), metadata[i,"SubID"], NA)
  })),
  "BD" = na.omit(sapply(1:nrow(metadata), function(i) {
    ifelse((sum(na.omit(metadata[i,setdiff(deleterious_dx,c("BD_ALL", "BD_unspecific", "BD_I", "BD_II", "Dementia", "MCI"))])) == 0) & (metadata[i,"BD_ALL"] == 1), metadata[i,"SubID"], NA)
  })),
  "PD" = na.omit(sapply(1:nrow(metadata), function(i) {
    ifelse((sum(na.omit(metadata[i,setdiff(deleterious_dx,c("PD", "PD_uncertain_plus_encephalitic", "Dementia", "MCI"))])) == 0) & (metadata[i,"PD"] == 1 | metadata[i,"PD_uncertain_plus_encephalitic"] == 1), metadata[i,"SubID"], NA)
  })),
  "DLBD" = na.omit(sapply(1:nrow(metadata), function(i) {
    ifelse((sum(na.omit(metadata[i,setdiff(deleterious_dx,c("DLBD", "Dementia", "MCI"))])) == 0) & (metadata[i,"DLBD"] == 1), metadata[i,"SubID"], NA)
  }))
)

individual_dx = deleterious_dx
individual_dx_counts = lapply(individual_dx, function(dx) {
  na.omit(sapply(1:nrow(metadata), function(i) {
    ifelse((sum(na.omit(metadata[i,setdiff(deleterious_dx,c(dx, "Dementia", "MCI"))])) == 0) & (metadata[i,dx] == 1), metadata[i,"SubID"], NA)
  }))
})
names(individual_dx_counts) = individual_dx
sapply(individual_dx_counts, function(dxx) length(dxx))

individual_dx = deleterious_dx
individual_dx_counts = lapply(individual_dx, function(dx) {
  na.omit(sapply(1:nrow(metadata), function(i) {
    ifelse((sum(na.omit(metadata[i,setdiff(deleterious_dx,c(dx, "AD", "Dementia", "MCI"))])) == 0) & (metadata[i,dx] == 1) & (metadata[i,"AD"] == 1), metadata[i,"SubID"], NA)
  }))
})
names(individual_dx_counts) = individual_dx
sapply(individual_dx_counts, function(dxx) length(dxx))

individual_dx2_counts = lapply(individual_dx, function(dx) {
  na.omit(sapply(1:nrow(metadata), function(i) {
    ifelse((metadata[i,dx] == 1), metadata[i,"SubID"], NA)
  }))
})
names(individual_dx2_counts) = individual_dx
sapply(individual_dx2_counts, function(dxx) length(dxx))

##########################
# Summary table of diagnostic groups

SubID_group_desc = list(
  "controls_neuropathological" = "CERAD=(1,2) and Braak=(0,1,2) and (secondary diagnosis allowed)",
  "controls_neuropathological_clinical" = "CERAD=(1,2) and Braak=(0,1,2) and (secondary diagnosis not allowed; MCI/dementia status not checked)",
  "controls_supercontrols" = "CERAD=(1), Braak=(0,1,2) and (secondary diagnosis not allowed;  MCI/dementia is checked too)",
  "controls_supercontrols_for_aging" = "CERAD=(1), Braak=(0,1,2) and (secondary diagnosis not allowed;  MCI/dementia is checked too); RUSH samples excluded",
  "AD" = "CERAD=(2,3,4), Braak=(3,4,5,6) and (must be clinically MCI or Dementia) and (secondary diagnosis not allowed); One of variables (CERAD, Braak, MSSM:CDR/RUSH:cogdx) can be NA",
  "AD_strict" = "CERAD=(4) and Braak=(3,4,5,6) and (must be clinically Dementia) and (secondary diagnosis not allowed)",
  "CERAD_1" = "(CERAD=1) and (secondary diagnosis not allowed;AD|MCI|Dementia status not checked)",
  "CERAD_2" = "(CERAD=2) and (secondary diagnosis not allowed;AD|MCI|Dementia status not checked)",
  "CERAD_3" = "(CERAD=3) and (secondary diagnosis not allowed;AD|MCI|Dementia status not checked)",
  "CERAD_4" = "(CERAD=4) and (secondary diagnosis not allowed;AD|MCI|Dementia status not checked)",
  "Braak_0" = "(Braak=0) and (secondary diagnosis not allowed;AD|MCI|Dementia status not checked)",
  "Braak_1" = "(Braak=1) and (secondary diagnosis not allowed;AD|MCI|Dementia status not checked)",
  "Braak_2" = "(Braak=2) and (secondary diagnosis not allowed;AD|MCI|Dementia status not checked)",
  "Braak_3" = "(Braak=3) and (secondary diagnosis not allowed;AD|MCI|Dementia status not checked)",
  "Braak_4" = "(Braak=4) and (secondary diagnosis not allowed;AD|MCI|Dementia status not checked)",
  "Braak_5" = "(Braak=5) and (secondary diagnosis not allowed;AD|MCI|Dementia status not checked)",
  "Braak_6" = "(Braak=6) and (secondary diagnosis not allowed;AD|MCI|Dementia status not checked)",
  "AD_DLBD" = "(AD) and (DLBD) and (secondary diagnosis not allowed;AD|MCI|Dementia status not checked)",
  "AD_Vascular" = "(AD) and (Vascular) and (secondary diagnosis not allowed;AD|MCI|Dementia status not checked)",
  "AD_SCZ" = "(AD) and (SCZ|Schizoaffective_bipolar|Schizoaffective_depressive) and (secondary diagnosis not allowed)",
  "SCZ" = "(SCZ|Schizoaffective_bipolar|Schizoaffective_depressive) and (secondary diagnosis not allowed)",
  "BD" = "(BD_unspecific|BD_I|BD_II) and (secondary diagnosis not allowed)",
  "PD" = "(PD|PD_uncertain_plus_encephalitic) and (secondary diagnosis not allowed)",
  "DLBD" = "(DLBD) and (secondary diagnosis not allowed)"
)

SubID_group_issues = list(
  "controls_neuropathological" = "If CERAD, Braak are not defined, sample is ignored",
  "controls_neuropathological_clinical" = "If CERAD, Braak, or CDR(MSSM)/cogdx(RUSH) are not defined, sample is ignored",
  "controls_supercontrols" = "If CERAD, Braak, or CDR(MSSM)/cogdx(RUSH) are not defined, sample is ignored",
  "controls_supercontrols_for_aging" = "If CERAD, Braak, or CDR(MSSM)/cogdx(RUSH) are not defined, sample is ignored",
  "AD" = "",
  "AD_strict" = "",
  "CERAD_1" = "",
  "CERAD_2" = "",
  "CERAD_3" = "",
  "CERAD_4" = "",
  "Braak_0" = "",
  "Braak_1" = "",
  "Braak_2" = "",
  "Braak_3" = "",
  "Braak_4" = "",
  "Braak_5" = "",
  "Braak_6" = "",
  "AD_DLBD" = "",
  "AD_Vascular" = "",
  "AD_SCZ" = "",
  "SCZ" = "Many SCZ+somethings (eating disorders) makes this 'pure' set smaller",
  "BD" = "Many BD+somethings (eating disorders) makes this 'pure' set smaller",
  "PD" = "",
  "DLBD" = ""
)

summary_df = do.call("cbind.data.frame", list(
  "name" = names(SubID_group_desc), 
  "size" = sapply(names(SubID_group_desc), function(cat) length(SubID_groups[[cat]])),
  "size_mssm" = sapply(names(SubID_group_desc), function(cat) sum(SubID_groups[[cat]] %in% metadata[metadata$Brain_bank=="MSSM","SubID"])),
  "size_rush" = sapply(names(SubID_group_desc), function(cat) sum(SubID_groups[[cat]] %in% metadata[metadata$Brain_bank=="RUSH","SubID"])),
  "size_hbcc" = sapply(names(SubID_group_desc), function(cat) sum(SubID_groups[[cat]] %in% metadata[metadata$Brain_bank=="HBCC","SubID"])),
  "size_mssm" = sapply(names(SubID_group_desc), function(cat) mean(metadata[SubID_groups[[cat]][SubID_groups[[cat]] %in% metadata[metadata$Brain_bank=="MSSM","SubID"]],"Age"] )),
  "size_rush" = sapply(names(SubID_group_desc), function(cat) mean(metadata[SubID_groups[[cat]][SubID_groups[[cat]] %in% metadata[metadata$Brain_bank=="RUSH","SubID"]],"Age"] )),
  "size_hbcc" = sapply(names(SubID_group_desc), function(cat) mean(metadata[SubID_groups[[cat]][SubID_groups[[cat]] %in% metadata[metadata$Brain_bank=="HBCC","SubID"]],"Age"] )),
  "desc" = unlist(SubID_group_desc),
  "notes" = unlist(SubID_group_issues)))

#write.csv(summary_df, file="~/Desktop/psychad_distribution.csv")

##
metadata$sum = 0
metadata$sum = sapply(1:nrow(metadata), function(i) sum(na.omit(as.integer(metadata[i,all_dx]))) )
table(metadata$sum)
#write.csv(metadata, file="~/Desktop/psychad_distribution_2.csv")
##

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

corrMatrix = metadata[c(SubID_groups$controls_neuropathological_clinical, SubID_groups$AD), c("CERAD", "Age", "CDRScore", "BRAAK_AD", "AD", "MCI", "Brain_bank")] # "PLAQUE"
#corrMatrix = corrMatrix[corrMatrix$Brain_bank=="MSSM",]
corrMatrix[corrMatrix$Brain_bank=="RUSH","CDRScore"] = sapply(which(corrMatrix$Brain_bank=="RUSH"), function(i) { if(corrMatrix[i,"AD"]==0 & corrMatrix[i,"MCI"] == 0) { 0 } else if(corrMatrix[i,"MCI"]) { 0.75 } else { 3 } })
#corrMatrix$CERAD = sapply(corrMatrix$CERAD, function(cerad) { if(is.na(cerad)) { NA } else if(cerad==2) { 4 } else if(cerad==4) { 2 } else { cerad } }) # skip this mess
corrMatrix = corrMatrix[,c("Age", "CERAD","CDRScore", "BRAAK_AD")] # "PLAQUE"
#corrMatrix = corrMatrix[complete.cases(corrMatrix),]
corrMatrixOut = data.frame(sapply(1:ncol(corrMatrix), function(x) sapply(1:ncol(corrMatrix), function(y) { cor.test(as.numeric(corrMatrix[,x]), as.numeric(corrMatrix[,y]), method="spearman")$estimate } )))
colnames(corrMatrixOut) = c("Age of death", "CERAD", "CDRScore", "BRAAK_AD") # "PLAQUE"
rownames(corrMatrixOut) = colnames(corrMatrixOut)
corrMatrixOut = get_upper_tri(corrMatrixOut)
corrMatrixOut$ID = rownames(corrMatrixOut)
corrMatrixOut = melt(corrMatrixOut, id.vars="ID", factorsAsStrings=T)
corrMatrixOut$variable = ordered(corrMatrixOut$variable,levels=c("Age of death", "CERAD", "PLAQUE", "CDRScore", "BRAAK_AD"))
corrMatrixOut$ID = ordered(corrMatrixOut$ID,levels=c("Age of death", "CERAD", "PLAQUE", "CDRScore", "BRAAK_AD"))

corr = ggplot(corrMatrixOut, aes(variable, ID, fill = value)) +
  geom_tile() + scale_fill_material("red") +
  coord_equal() + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x=element_text(colour = "black"), axis.text.y=element_text(colour = "black"), axis.ticks.x=element_blank(),
                                     axis.title.y=element_blank(),  axis.ticks.y=element_blank())
corr = corr + geom_text(aes(variable, ID, label = round(value, 2)), color = "black", size = 3) 

metadata = metadata[!((metadata$Brain_bank == "RUSH") & is.na(metadata$Sex)),]
sum(!is.na(metadata[metadata$Brain_bank=="MSSM","PLAQUE"])) / length(metadata[metadata$Brain_bank=="MSSM","PLAQUE"])
sum(!is.na(metadata[metadata$Brain_bank=="RUSH","PLAQUE"])) / length(metadata[metadata$Brain_bank=="RUSH","PLAQUE"])
sum(!is.na(metadata[metadata$Brain_bank=="MSSM","CERAD"])) / length(metadata[metadata$Brain_bank=="MSSM","CERAD"])
sum(!is.na(metadata[metadata$Brain_bank=="RUSH","CERAD"])) / length(metadata[metadata$Brain_bank=="RUSH","CERAD"])
any(sapply(which(metadata$Brain_bank=="MSSM"), function(i) { is.na(metadata[i,"CERAD.1"]) & is.na(metadata[i,"CERAD_2"]) }))
sum(!is.na(metadata[metadata$Brain_bank=="RUSH","BRAAK_AD"])) / length(metadata[metadata$Brain_bank=="RUSH","BRAAK_AD"])

#
print("Samples by Brain_bank")
table(metadata$Brain_bank)

length(metadataAll$SubID)
nrow(metadataAll[!is.na(metadataAll$SNParray_PsychAD),])
nrow(metadataAll[!is.na(metadataAll$Imaging_XENum),])
#table(unique(MSSM_imaging_meta$SubID) %in% MSSM_imaging_meta$SubID)
#table(unique(MSSM_imaging_meta$SubID) %in% metadata$SubID)

#
print("Number of clinically diagnosed diseases per donor")
metadata$sum = sapply(1:nrow(metadata), function(i) sum(na.omit(as.integer(metadata[i, unlist(DIAGNOSIS_CLASSIFICATION)]))))
table(metadata$sum)

saveRDS(SubID_groups, file="/sc/arion/projects/roussp01a/jaro/clinical_metadata_sampleSets.RDS")

# Copy merge & querying code for metadata to public github
cmd_copy_metadata = paste0("cp ", unlist(PSYCHAD_METADATA_CODES), " ", PSYCHAD_METADATA_GITHUB_MINERVA)
sapply(cmd_copy_metadata, system)

## upload on synapse
file <- synStore(File(path="/sc/arion/projects/roussp01a/jaro/clinical_metadata_sampleSets.RDS", parent="syn22399913"), used="https://github.com/DiseaseNeuroGenomics/psychad_metadata/blob/main/psychad_clinical_metadata_querying.R")
