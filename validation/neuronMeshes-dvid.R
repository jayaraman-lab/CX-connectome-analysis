


#' Gets a neuron mesh from a dvid server and returns it as a mesh3d object
#' @param bodyid the bodyid of the neuron to return
#' @param server the dvid server to query (see "https://dvid.io/")
#' @param node the uuid of the node to query
getNeuronMeshDvid <- function(bodyid,server="https://hemibrain-dvid.janelia.org/",node=newNode){
  raw <- getNeuronMesh_raw(bodyid,server=server,node=node)
  ngMeshToMesh(raw)
}

getNeuronMesh_raw <- function(bodyid,server=NULL,node=NULL){
  
  neuronQuery <-  paste0(server,"api/node/",node,"/segmentation_meshes/key/",bodyid,".ngmesh")
  
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
