## Most recent node for emdata4 77309
## Node for version 0.9 34d7a
newNode <- "20631"
oldNode <- "34d7a"

getNeuronMesh <- function(bodyid,cloud=FALSE,node=newNode){
  raw <- getNeuronMesh_raw(bodyid,cloud=cloud,node=node)
  ngMeshToMesh(raw)
}

getNeuronMesh_raw <- function(bodyid,cloud=FALSE,node=newNode){
  if (cloud)
    serv <- "https://hemibrain-dvid2.janelia.org/"
  else
    serv <- "https://emdata4.int.janelia.org:8900/"
  
  neuronQuery <-  paste0(serv,"api/node/",node,"/segmentation_meshes/key/",bodyid,".ngmesh")
  
  r <- httr::GET(neuronQuery)
  httr::content(r)
}

ngMeshToMesh <- function(rawMesh){
  fullLength <- as.integer(length(rawMesh)/4)
  nVertices <- readBin(rawMesh,integer())
  xyzVertices <- matrix(readBin(rawMesh[5:length(rawMesh)],numeric(),size=4,n=nVertices*3),nrow=3)/8
  faces <- matrix(readBin(rawMesh[(5+nVertices*3*4):length(rawMesh)],integer(),n=fullLength-nVertices*3-1),nrow=3)+1
  rgl::tmesh3d(xyzVertices,faces,homogeneous = FALSE)
}
