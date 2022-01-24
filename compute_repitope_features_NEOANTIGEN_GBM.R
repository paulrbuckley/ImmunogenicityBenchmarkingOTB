# Title     : TODO
# Objective : TODO
# Created by: paulbuckley
# Created on: 22/03/2021

options(java.parameters="-Xmx60G")  ## allow JAVA to use large memory space
library(tidyverse)
library(data.table)
library(Repitope)
# Create fragment library using default settings.
fragLibDT <- CPP_FragmentLibrary(TCRSet_Public, fragLenSet=3:11, maxFragDepth=100000, seedSet=1:5)
fst::write_fst(fragLibDT, "FragLibrary_2021.fst", compress=0)
FullDataset = Epitope_Import(OtherFileNames = list("GBM_TestData_ForRepitope.csv"))

# Compute the features of our test peptides using default settings
# # Features [MHC-I]
 featureDFList_MHCI_EPS <- Features(
   peptideSet=unique(c(FullDataset$Peptide)),
   fragLib="FragLibrary_2021.fst",
   aaIndexIDSet="all",
   fragLenSet=3:8,
   fragDepth=10000,
   fragLibType="Weighted",
   seedSet=1:5,                                   ## must be the same random seeds used for preparing the fragment library
   coreN=parallel::detectCores()-2,        ## parallelization
   tmpDir="NEOANTIGEN_GBM/REPITOPE_OTB/temp_v3"   ## where intermediate files are stored
 )


saveFeatureDFList(featureDFList_MHCI_EPS, "FeatureDF_NEOANTIGEN_GBM_EPs")