## ---- eval=FALSE---------------------------------------------------------
#  source("http://bioconductor.org/biocLite.R")
#  biocLite(MethylMix)

## ---- eval=FALSE---------------------------------------------------------
#  cancerSite <- "OV"
#  targetDirectory <- paste0(getwd(), "/")
#  GetData(cancerSite, targetDirectory)

## ---- eval=FALSE---------------------------------------------------------
#  cancerSite <- "OV"
#  targetDirectory <- paste0(getwd(), "/")
#  
#  library(doParallel)
#  cl <- makeCluster(5)
#  registerDoParallel(cl)
#  GetData(cancerSite, targetDirectory)
#  stopCluster(cl)

## ---- eval=FALSE, tidy=TRUE----------------------------------------------
#  cancerSite <- "OV"
#  targetDirectory <- paste0(getwd(), "/")
#  
#  cl <- makeCluster(5)
#  registerDoParallel(cl)
#  
#  # Downloading methylation data
#  METdirectories <- Download_DNAmethylation(cancerSite, targetDirectory)
#  # Processing methylation data
#  METProcessedData <- Preprocess_DNAmethylation(cancerSite, METdirectories)
#  # Saving methylation processed data
#  saveRDS(METProcessedData, file = paste0(targetDirectory, "MET_", cancerSite, "_Processed.rds"))
#  
#  # Downloading gene expression data
#  GEdirectories <- Download_GeneExpression(cancerSite, targetDirectory)
#  # Processing gene expression data
#  GEProcessedData <- Preprocess_GeneExpression(cancerSite, GEdirectories)
#  # Saving gene expression processed data
#  saveRDS(GEProcessedData, file = paste0(targetDirectory, "GE_", cancerSite, "_Processed.rds"))
#  
#  # Clustering probes to genes methylation data
#  METProcessedData <- readRDS(paste0(targetDirectory, "MET_", cancerSite, "_Processed.rds"))
#  res <- ClusterProbes(METProcessedData[[1]], METProcessedData[[2]])
#  
#  # Putting everything together in one file
#  toSave <- list(METcancer = res[[1]], METnormal = res[[2]], GEcancer = GEProcessedData[[1]], GEnormal = GEProcessedData[[2]], ProbeMapping = res$ProbeMapping)
#  saveRDS(toSave, file = paste0(targetDirectory, "data_", cancerSite, ".rds"))
#  
#  stopCluster(cl)

## ---- tidy=TRUE----------------------------------------------------------
library(MethylMix)
library(doParallel)
data(METcancer)
data(METnormal)
data(GEcancer)
head(METcancer[, 1:4])
head(METnormal)
head(GEcancer[, 1:4])

## ---- tidy=TRUE, warning=F-----------------------------------------------
MethylMixResults <- MethylMix(METcancer, GEcancer, METnormal)

## ---- tidy=TRUE, eval=FALSE----------------------------------------------
#  library(doParallel)
#  cl <- makeCluster(5)
#  registerDoParallel(cl)
#  MethylMixResults <- MethylMix(METcancer, GEcancer, METnormal)
#  stopCluster(cl)

## ---- tidy=TRUE----------------------------------------------------------
MethylMixResults$MethylationDrivers
MethylMixResults$NrComponents
MethylMixResults$MixtureStates
MethylMixResults$MethylationStates[, 1:5]
MethylMixResults$Classifications[, 1:5]
# MethylMixResults$Models

## ---- tidy=TRUE, eval=F--------------------------------------------------
#  # Plot the most famous methylated gene for glioblastoma
#  plots <- MethylMix_PlotModel("MGMT", MethylMixResults, METcancer)
#  plots$MixtureModelPlot

## ---- tidy=TRUE, eval=F--------------------------------------------------
#  # Plot MGMT also with its normal methylation variation
#  plots <- MethylMix_PlotModel("MGMT", MethylMixResults, METcancer, METnormal = METnormal)
#  plots$MixtureModelPlot

## ---- tidy=TRUE, eval=F--------------------------------------------------
#  # Plot a MethylMix model for another gene
#  plots <- MethylMix_PlotModel("ZNF217", MethylMixResults, METcancer, METnormal = METnormal)
#  plots$MixtureModelPlot

## ---- tidy=TRUE, eval=F--------------------------------------------------
#  # Also plot the inverse correlation with gene expression (creates two separate plots)
#  plots <- MethylMix_PlotModel("MGMT", MethylMixResults, METcancer, GEcancer, METnormal)
#  plots$MixtureModelPlot
#  plots$CorrelationPlot

## ---- eval = FALSE, tidy=TRUE--------------------------------------------
#  # Plot all functional and differential genes
#  for (gene in MethylMixResults$MethylationDrivers) {
#      MethylMix_PlotModel(gene, MethylMixResults, METcancer, METnormal = METnormal)
#  }

## ---- tidy=TRUE, echo = FALSE--------------------------------------------
sessionInfo()

