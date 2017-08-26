### R code from vignette source 'DriverNet-Overview.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(width=60)


###################################################
### code chunk number 2: preliminaries
###################################################
library(DriverNet)

data(samplePatientMutationMatrix)
data(samplePatientOutlierMatrix)
data(sampleInfluenceGraph)
data(sampleGeneNames)



###################################################
### code chunk number 3: options
###################################################
options(width=60)


###################################################
### code chunk number 4: preliminaries
###################################################
# The main function to compute drivers
driversList = computeDrivers(samplePatientMutationMatrix, 
              samplePatientOutlierMatrix,sampleInfluenceGraph,
              outputFolder=NULL, printToConsole=FALSE)

drivers(driversList)[1:10]


###################################################
### code chunk number 5: options
###################################################
options(width=60)


###################################################
### code chunk number 6: preliminaries
###################################################
# random permute the gene labels to compute p-values
randomDriversResult = computeRandomizedResult(
  patMutMatrix=samplePatientMutationMatrix,
  patOutMatrix=samplePatientOutlierMatrix, 
  influenceGraph=sampleInfluenceGraph,
  geneNameList= sampleGeneNames, outputFolder=NULL, 
  printToConsole=FALSE,numberOfRandomTests=20, weight=FALSE, 
  purturbGraph=FALSE, purturbData=TRUE)


###################################################
### code chunk number 7: options
###################################################
options(width=60)


###################################################
### code chunk number 8: preliminaries
###################################################
# Summarize the results
res = resultSummary(driversList, randomDriversResult, 
      samplePatientMutationMatrix,sampleInfluenceGraph,
      outputFolder=NULL, printToConsole=FALSE)
res[1:2,]


