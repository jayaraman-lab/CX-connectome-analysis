# This file contain functions for:
# 0) A number of modular functions from 1-4 for performing many PCA/projection related manipulations of synapse and roi locations.
# 1) Performing PCA on meshes or synapse locations to get primary axes
# 2) Projecting a set of synapses or mesh points onto new axes
# 3) Plotting the distribution of synapses over ROI outlins
# 4) Getting synapse locations in a loop so the query doesnt time out


MeshOutline <- function(Mesh, Plane){
  
  if (Plane=="XY"){
    meshOutline <- ahull(x=Mesh$X, y=-Mesh$Y,alpha=100)
  } else if (Plane=="XZ"){
    meshOutline <- ahull(x=Mesh$X, y=Mesh$Z,alpha=100)
  } else if (Plane=="YZ"){
    meshOutline <- ahull(x=Mesh$Y, y=Mesh$Z,alpha=100)
  }
  
  outline = data.frame(meshOutline$arcs)
  outline = bind_rows(outline, outline[1,]) # close the shape
  
  return(outline)
}


MeshPoints <- function(roiMesh){
  roiMeshPts = data.frame(dotprops(roiMesh)$points)
  names(roiMeshPts) <- c("x","y","z")
  return(roiMeshPts)
}


getCOM <- function(pointsXYZ){
  return(c(mean(pointsXYZ$x),mean(pointsXYZ$y),mean(pointsXYZ$z)))
}

resetOrigin <- function(pointsXYZ, origin){
  pointsXYZ = transform(pointsXYZ, x = x-origin[1], y = y-origin[2], z = z-origin[3])
  return(pointsXYZ)
}


makeRotMatXY <- function(angle){
  # Rotate the first two PCs to align ROI in direction we want 
  angRad=angle/360*2*pi
  RotationMatrix=matrix( c(cos(angRad), sin(angRad), 0, 
                           -sin(angRad),cos(angRad),0,
                           0,0,1) , nrow = 3, ncol = 3)
  return(RotationMatrix)
}

makeRotMatYZ <- function(angle){
  # Rotate the first two PCs to align ROI in direction we want 
  angRad=angle/360*2*pi
  RotationMatrix=matrix( c(1,0,0,
                           0,cos(angRad), sin(angRad), 
                           0,-sin(angRad),cos(angRad)), 
                         nrow = 3, ncol = 3)
  return(RotationMatrix)
}

makeRotMatXZ <- function(angle){
  # Rotate the first two PCs to align ROI in direction we want 
  angRad=angle/360*2*pi
  RotationMatrix=matrix( c(cos(angRad),0,sin(angRad),
                           0,1,0,
                           -sin(angRad),0,cos(angRad)), 
                         nrow = 3, ncol = 3)
  return(RotationMatrix)
}

covPCA <- function(pointsXYZ){
  # calculate covariance matrix
  cov = data.frame(cov(pointsXYZ))
  # calculate eigenvectors and values (columns contain PCs, ranked from most variance accounted for to least)
  covEigen = eigen(cov) 
  return(covEigen)
}


changeBasis <- function(pointsXYZ, covEigen){
  # points should be centered at origin
  
  #Convert all points to the coordinate system defined by covariance matrix eigen vectors
  pts = matrix(0, nrow = length(pointsXYZ$x), ncol = 3)
  for (i in seq(1,length(pointsXYZ$x))) {
    pt = as.numeric(pointsXYZ[i,])  %*% covEigen$vectors
    pts[i,] = pt
  }
  pointsNewBase = data.frame(x = pointsXYZ$x, y = pointsXYZ$y, z = pointsXYZ$z,
                             X = pts[,1], Y = pts[,2], Z = pts[,3])
  
  return(pointsNewBase)
}



changeBasis_df <- function(Input_DF, NewAxes) {
  
  
  # Project just the x,y,z into new space
  pointsXYZ=data.frame(x=Input_DF$x, y=Input_DF$y, z=Input_DF$z)
  NewProjection=changeBasis(pointsXYZ, NewAxes)
  
  # Populate new data frame
  NewProjection_DF = Input_DF
  NewProjection_DF$x=NewProjection$x
  NewProjection_DF$y=NewProjection$y
  NewProjection_DF$z=NewProjection$z
  NewProjection_DF$X=NewProjection$X
  NewProjection_DF$Y=NewProjection$Y
  NewProjection_DF$Z=NewProjection$Z
  
  
  return(NewProjection_DF)
  
}


# Loop over bodyids and get synapse locations to avoid queries that time out
GetSynapseLocs <- function(BodyIDs, ROI, GROUP) {
  
  
  if (GROUP==TRUE){
    
    if (length(BodyIDs)>=25){Chunks=25}else{Chunks=1} # chunk of 25 works for FB
    
    Starts=seq(from = 1, to = floor(length(BodyIDs)/Chunks)*Chunks+1, by = Chunks)
    Stops=c(Starts[2:length(Starts)]-1, length(BodyIDs))
    
    
    for (bbb in 1:length(Starts)){
      print(bbb)
      if (bbb==1){
        
        if (ROI=="all"){SynLocs =  neuprint_get_synapses(BodyIDs[ Starts[bbb]:Stops[bbb] ])
        } else {SynLocs =  neuprint_get_synapses(BodyIDs[ Starts[bbb]:Stops[bbb] ], roi = ROI)}
        
        SynLocs =  mutate(SynLocs, type=neuprint_get_meta(bodyid)$type, partner_type=neuprint_get_meta(partner)$type,
                          x=as.numeric(x),y=as.numeric(y),z=as.numeric(z))
      } else {
        
        if (ROI=="all"){TempLocs =  neuprint_get_synapses(BodyIDs[ Starts[bbb]:Stops[bbb] ])
        } else {TempLocs =  neuprint_get_synapses(BodyIDs[ Starts[bbb]:Stops[bbb] ], roi = ROI)} 
        
        TempLocs =  mutate(TempLocs, type=neuprint_get_meta(bodyid)$type, partner_type=neuprint_get_meta(partner)$type,
                           x=as.numeric(x),y=as.numeric(y),z=as.numeric(z))
        SynLocs=rbind(SynLocs,TempLocs)
        remove(TempLocs)
      }
      
    }
    
  } else {
    
    if (ROI=="all"){SynLocs =  neuprint_get_synapses(BodyIDs)
    } else {SynLocs =  neuprint_get_synapses(BodyIDs, roi = ROI)}
    
    SynLocs =  mutate(SynLocs, type=neuprint_get_meta(bodyid)$type, partner_type=neuprint_get_meta(partner)$type,
                      x=as.numeric(x),y=as.numeric(y),z=as.numeric(z))
    
  }
  
  SynLocs$prepost[SynLocs$prepost==0]="Output"
  SynLocs$prepost[SynLocs$prepost==1]="Input"
  
  return(SynLocs)
  
}














  