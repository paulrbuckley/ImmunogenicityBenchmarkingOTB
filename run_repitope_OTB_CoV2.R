# Title     : TODO
# Objective : TODO
# Created by: paulbuckley

options(java.parameters="-Xmx60G")  ## allow JAVA to use large memory space
library(tidyverse)
library(data.table)
library(Repitope)
# Read in epitopes for testing
FullDataset = Epitope_Import(OtherFileNames = list("CoV2_TestData_ForRepitope.csv"))
# Read in Repitope MHC-I model features, for training the Repitope model as published by Ogishi et al.
featureDF_MHCI <- fst::read_fst("FeatureDF_MHCI_Weighted.10000_RepitopeV3.fst", as.data.table=T)
# Read in features of epitopes for testing
featureDF_COV2_Test_Peps = fst::read_fst("FeatureDF_COV2_EPsWeighted.10000.fst", as.data.table=T)
# 764 unique peptides reported below as several peptides appear multiple times in the context of different HLA.
if(FullDataset %>% filter(Peptide %in% featureDF_COV2_Test_Peps$Peptide) %>% nrow == FullDataset %>% nrow){
  print("Full peptides are in feature DF")
}

# Train the Repitope MHC-I Out the box Model
MHCI_HUMAN_MODEL=Immunogenicity_TrainModels(
  featureDF=featureDF_MHCI[Peptide%in%MHCI_Human$Peptide,],
  metadataDF=MHCI_Human[,.(Peptide, Immunogenicity)],
  featureSet=MHCI_Human_MinimumFeatureSet,
  seedSet = 1:5,
  coreN = parallel::detectCores(logical = T)-4
)
# Extrapolate model to unseen peptides
ALL_PREDICTIONS = Immunogenicity_Predict(list(featureDF_COV2_Test_Peps[Peptide%in%FullDataset$Peptide,]),MHCI_HUMAN_MODEL)
# Write out the predictions
write.csv(ALL_PREDICTIONS$ScoreDT_1,file= "SARS_COV_2_DATASET/REPITOPE_OTB/ALL_PREDICTIONS.csv",quote=F,row.names = F)
