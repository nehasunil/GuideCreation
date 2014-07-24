#Input Variables

#dataSet: Users can import other data sets if they are interested in other data.
dataSet ="Roadmap + Encode"

#focus: Name of cell type or individual
#The users add the name of the specific cell type or individual or import if necessary 
focus = NULL

#chromState: Users choose whether they are interested in looking into enhancers, promoters, 
#genes, or a combination by choosing chromatin states of interest.
chromState = NULL

#geneOntology: Based on the testing method, users can input what kinds of genes they are 
#interested in. Users input biological process ID(s).
geneOntology = "GO:0040007"

#Sampling Method
#Users can choose how to sample guide sequence with multiple inputs

#guideNumber: integer value for the number of guides to choose.
guideNumber=100000

#guidesPerRegion: integer value of the number of guides for each sampled region
#recalculates if input value is given for numberOfRegions 
#guidesPerRegion = guideNumber/numberOfRegions
guidesPerRegion=5

#numberOfRegions: integer value for the number of regions to be sampled
numberOfRegions=25000
#recalculates if input value is given for guidesPerRegion
#numberOfRegions = guideNumber/guidesPerRegion

#fractionCellTypeSpecific: decimal value from 0 to 1 for what percent of the 
#sampled region should be specific to the cell type
fractionCellTypeSpecific = 1

#weightingForSampling: vector with weights for each type of region that user can input
#Default: calculated based on fractionCellTypeSpecific
weightingForSampling=NULL

#guidesPerRegionVector: vector with number of guides for each region
#sum must equal guideNumber
#Default: equal weighting (guidesPerRegion) for all values
guidesPerRegionVector=NULL

#minFractionFocusInCluster: decimal value for the minimum percentage of a cluster that
#has to be unique to the focus in order for cluster to represent focus-specific regions
minFractionFocusInCluster=.80


#Other Variables
#cluster: the cluster that contains the highest percentage of cell type specific enhancers
cluster=NULL

#sampledRegions: the regions that have been sampled to make guides from
sampledRegions=NULL

#Methods
#accessors
#getDataSet() - returns name of data set (String)
getDataSet <- function()
{
  return(dataSet)
}
#getFocus() - returns name of focus (String)
getFocus <- function()
{
  return(focus)
}
#getChromState() - returns vector with all the chromatin states chosen (vector of Strings)
getChromState <- function()
{
  return(chromState)
}
#getGeneOntology() - returns vector of strings with gene ontology selections
getGeneOntology <- function()
{
  return(geneOntology)
}
#getGuideNumber() - returns number of guides to create (int)
getGuideNumber <- function()
{
  return(guideNumber)
}
#getGuidesPerRegion() - returns number of guides per sampled region
getGuidesPerRegion <- function()
{
  return(guidesPerRegion)
}
#getNumberOfRegions() - returns number of regions to be chosen
getNumberOfRegions <- function()
{
  return(numberOfRegions)
}
#getFractionCellTypeSpecific() - returns fraction of regions that should be cell type specific
getFractionCellTypeSpecific <- function()
{
  return(fractionCellTypeSpecific)
}
#getWeightingForSampling() - return weighting for sampling vector
getWeightingForSampling <- function()
{
  return(weightingForSampling)
}
#getGuidesPerRegionVector() - returns guides per region vector
getGuidesPerRegionVector <- function()
{
  return(getGuidesPerRegionVector)
}

#modifiers:
#setDataSet(String data) - sets dataSet to data
setDataSet <- function(data)
{
  dataSet=data
}
#setFocus(String focusName) - sets focus to focusName
setFocus <- function(focusName)
{
  focus=focusName
}
#setChromState(vector chromatinState) adds chromatinState to chromState
setChromState <- function(chromatinState)
{
  chromState=chromatinState
}
#addChromState(String chromatinState) adds chromatinState to chromState
addChromState <- function(chromatinState)
{
  chromState=c(chromState,chromatinState)
}
#removeChromState(String chromatinState) removes chromatinState from chromState
removeChromState <- function(chromatinState)
{
  index=match(chromatinState,chromState)
  if(!is.na(index))
  {
    chromState=chromState[-index]
  }
  else
  {
    print(paste(chromatinState,"cannot be removed"))
  }
}
#setGeneOntology(vector go) - sets geneOntology to go
setGeneOntology <- function(go)
{
  geneOntology=go
}
#addGeneOntology(String go) - adds go to geneOntology
addGeneOntology <- function(go)
{
  geneOntology=c(chromeState,focusName)
}
#removeGeneOntology(String go) - removes go from geneOntology
removeGeneOntology <- function(go)
{
  index=match(go,geneOntology)
  if(!is.na(index))
  {
    geneOntology=go[-index]
  }
  else
  {
    print(paste(go,"cannot be removed"))
  }
}
#setGuideNumber(int guideNum) - sets guideNumber to guideNum
setGuideNumber <- function(guideNum)
{
  guideNumber=guideNum
}
#setGuidesPerRegion(int numGuidesPerRegion) sets guidesPerRegion to numGuidesPerRegion and changes numberOfRegions accordingly
setGuidesPerRegion <- function(numGuidesPerRegion)
{
  guidesPerRegion = numGuidesPerRegion
  numberOfRegions = floor(guideNumber/guidesPerRegion)
}
#setNumberOfRegions(int numRegions) - sets numberOfRegions to numRegions and changes guidesPerRegion accordingly
setNumberOfRegions <- function(numRegions)
{
  numberOfRegions = numRegions
  guidesPerRegion= floor(guideNumber/numberOfRegions)
}
#setFractionCellTypeSpecific(Double fractionCTSpecific) - sets fractionCellTypeSpecific to fractionCTSpecific
setFractionCellTypeSpecific <- function(fractionCTSpecific)
{
  fractionCellTypeSpecific = fractionCTSpecific
}
#setWeightingForSampling(vector weightForSampling)
setWeightingForSampling <- function(weightForSampling)
{
  fractionCellTypeSpecific = fractionCTSpecific
}
#setGuidesPerRegionVector(vector GuidesPerRegionVect)
setGuidesPerRegionVector <- function(guidesPerRegionVect)
{
  guidesPerRegionVector = guidesPerRegionVect
}
#createInfoMatrix(vector cols) - returns a sparse matrix with each region grouped by cluster as the rows, 
#and the attributes in the parameter cols as the columns. All the valid column options include cluster, 
#number of transcription factor binding sites, and gene ontology.
createInfoMatrix <- function(cols)
{
}
#IntersectionTFSites(dataframe(BED files) regions) - returns a sparse matrix with each transcription factor 
#and its binding sites in regions 
intersectionTFSites <- function(regions)
{
  #bedtools for R
  #find intersections
}
#selectCluster() - returns cluster (can only be done with Roadmap data) with highest percentage of its regions 
#from the cell type of focus and sets cluster instance variable to the cluster (set min percentage)
selectCluster <- function()
{
  #find percent in each cluster
  #check if any percent is greater than minFractionFocusInCluster
  #set cluster instance variable
  return(cluster)
}
#calcWeightingForSampling() - calculates weighting of regions based on fraction of cell type specific
calcWeightingForSampling <- function()
{
  focusedCount= round(fractionCellTypeSpecific *guideNumber)
 
}
#sampleRegions(dataframe(BED files) remRegion) - samples specified number of regions using input values 
#ignoring regions in remRegion using weightingForSampling, sets sampledRegions to the regions sampled, 
#and returns sampledRegions
sampleRegions <- function(remRegions)
{
}
#createGuides(dataframe(BED files) regions) - creates guides for the given regions,
#the number of guides per region equals the input value guidesPerRegion
createGuides <- function(regions)
{
}




