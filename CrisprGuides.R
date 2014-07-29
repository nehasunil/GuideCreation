#Input Variables

#dataSet: Users can import other data sets if they are interested in other data. 
dataSet ="Roadmap + Encode"

#dataDirectory: Users can import directory of data here if they are interested in other data
dataDirectory = NULL

#clusterDirectory: Users can import directory of clusters here if they are interested in other data
clusterDirectory = NULL

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
  #bedtools for R?
  #find intersections
  #
}
#selectCluster() - returns cluster (can only be done with Roadmap data) with highest percentage of its regions 
#from the cell type of focus and sets cluster instance variable to the cluster (set min percentage)
selectCluster <- function()
{
  files <- list.files(path=clusterDirectory)
  
  clust=read.table(paste(clusterDirectory,files[1],sep=''))
  count=0
  for(k in 1:nrow(clust))
  {
    if(!is.na(match(paste(as.character(clust[k,1]),":",clust[k,2],"-",clust[k,3],sep=''), paste(as.character(focus[,4])))))
    { count=count+1 }
  }
  percentage=count/nrow(clust)
  k562_dist=cbind(count,percentage)
  
  for(i in 2:length(files))
  {
    clust=read.table(paste(clusterDirectory,files[i],sep=''))
    count=0
    for(k in 1:nrow(clust))
    {
      if(!is.na(match(paste(as.character(clust[k,1]),":",clust[k,2],"-",clust[k,3],sep=''), paste(as.character(focus[,4])))))
      { count=count+1 }
    }
    percentage=count/nrow(clust)
    k562_dist=rbind(k562_dist,c(count,percentage))
  }
  
  rownames(k562_dist)=c(1:length(files))
  
  if (max(k562_dist[,2])>minFractionFocusInCluster)
  {
    maxInd=which.max(k562_dist[,2])
    cluster=read.table(paste(clusterDirectory,files[maxInd],sep=''))
    return(cluster)
  }
}
#calcWeightingForSampling() - calculates weighting of regions based on fraction of cell type specific
#parameter = data frame with gene ontology in first column and weighting for each GO in second column
calcWeightingForSampling <- function(go)
{
  if(nrow(go)>0)
  {
    weighting = rep(0,times=nrow(focus))
    for(i in 1:nrow(focus))
    {
      goIndex=match(focus$GeneOntology[i],go[,1])
      if(!is.na(goIndex))
      {
          if(length(goIndex)==1)
          {
            weighting[i]=go[goIndex,2]
          }
          else
          {
            weights=go[goIndex,2]
            weighting[i]=max(weights)
          }
      }
    }
    return(weighting)
  }
  else
  {
    return(rep(1,times=nrow(focus)))
  }
}
#sampleRegions(dataframe(BED files) remRegion) - samples specified number of regions using input values 
#ignoring regions in remRegion given indices of guides of interest (weighted), 
#
#sets sampledRegions to the regions sampled
#and returns sampledRegions
sampleRegions <- function(weighted)
{
  numCTSpecific=fractionCellTypeSpecific*guideNumber
  numNonSpecific=guideNumber-numCTSpecific
  
  focusWithWeighting=cbind(focus,weighted)
  
  ctSpecific=selectCluster()
  #for each cluster value, look for weighting value on focus and append
  ctSpecific$Weighting=0
  for(i in 1:nrow(ctSpecific))
  {
    index=match(ctSpecific[i,4],focusWithWeighting[,4])
    if(!is.na(index))
    {
      ctSpecific[i,ncol(ctSpecific)]=as.numeric(focusWithWeighting[index,ncol(focusWithWeighting)])
    }
  }
    
  if (nrow(ctSpecific) < numCTSpecific)
  {
    numCTspecific=nrow(ctSpecific)
    numNonSpecific = guideNumber-numCTSpecific
    print("There are not enough cell type specific regions")
    print(paste(numCTspecific,"cell type specific regions will be sampled"))
    print(paste(numNonSpecific,"nonspecific regions will be sampled"))
  }
  specificSample=ctSpecific[sample(1:nrow(ctSpecific),numCTSpecific),prob=ctSpecific[,nrow(ctSpecific)]]
  
  rest=focusWithWeighting
  for(i in 1:nrow(ctSpecific))
  {
    index=match(paste(ctSpecific[i,1],":",ctSpecific[i,2],"-",ctSpecific[i,3],sep=''),rest[,4])
    if(!is.na(index))
    {
      rest[-index,]
    }
  }
  nonSpecificSample=rest[sample(1:nrow(rest),numNonSpecific),prob=rest[,ncol(rest)]]
  return (cbind(specificSample,nonSpecificSample))
}
#createGuides(dataframe(BED files) regions) - creates guides for the given regions,
#the number of guides per region equals the input value guidesPerRegion
createGuides <- function(regions)
{
}




