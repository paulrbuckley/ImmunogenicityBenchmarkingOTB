# Title     : TODO
# Objective : TODO
# Created by: paulbuckley
# Created on: 22/03/2021

options(java.parameters="-Xmx60G")  ## allow JAVA to use large memory space
library(tidyverse)
library(data.table)
library(Repitope)

FullDataset = Epitope_Import(OtherFileNames = list("TESLA_TestData_ForRepitope.csv"))

# Read in feature set for training MHCI Human Repitope model.
featureDF_MHCI <- fst::read_fst("FeatureDF_MHCI_Weighted.10000_RepitopeV3.fst", as.data.table=T)
# Read in feature set for the test GBM data.
featureDF_Test_Peps = fst::read_fst("FeatureDF_NEOANTIGEN_TESLA_EPsWeighted.10000.fst", as.data.table=T)

# Generate Out the box MHCI Human model
MHCI_HUMAN_MODEL=Immunogenicity_TrainModels(
  featureDF=featureDF_MHCI[Peptide%in%MHCI_Human$Peptide,],
  metadataDF=MHCI_Human[,.(Peptide, Immunogenicity)],
  featureSet=MHCI_Human_MinimumFeatureSet,
  seedSet = 1:5,
  coreN = parallel::detectCores(logical = T)-4
)

ALL_PREDICTIONS = Immunogenicity_Predict(list(featureDF_Test_Peps[Peptide%in%FullDataset$Peptide,]),MHCI_HUMAN_MODEL)
write.csv(ALL_PREDICTIONS$ScoreDT_1,file= "NEOANTIGEN_TESLA/REPITOPE_OTB/ALL_PREDICTIONS.csv",quote=F,row.names = F)
