# Title     : TODO
# Objective : TODO
# Created by: paulbuckley

options(java.parameters="-Xmx60G")  ## allow JAVA to use large memory space
library(tidyverse)
library(data.table)
library(Repitope)
# Create fragment library using default settings
fragLibDT <- CPP_FragmentLibrary(TCRSet_Public, fragLenSet=3:11, maxFragDepth=100000, seedSet=1:5)
fst::write_fst(fragLibDT, "FragLibrary_2021.fst", compress=0)
# Read in the epitope dataset
FullDataset = Epitope_Import(OtherFileNames = list("CoV2_TestData_ForRepitope.csv"))

# Compute the features of our test peptides using default settings
# # Features [MHC-I]
 featureDFList_MHCI_COV2_EPS <- Features(
   peptideSet=unique(c(FullDataset$Peptide)),
   fragLib="FragLibrary_2021.fst",
   aaIndexIDSet="all",
   fragLenSet=3:8,
   fragDepth=10000,
   fragLibType="Weighted",
   seedSet=1:5,                                   ## must be the same random seeds used for preparing the fragment library
   coreN=parallel::detectCores()-2,        ## parallelization
   tmpDir="REPITOPE_OTB/temp_v5"   ## where intermediate files are stored
 )


saveFeatureDFList(featureDFList_MHCI_COV2_EPS, "FeatureDF_COV2_EPs")