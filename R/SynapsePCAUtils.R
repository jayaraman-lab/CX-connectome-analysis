# This file contain functions for:
# 0) A number of modular functions from 1-4 for performing many PCA/projection related manipulations of synapse and roi locations.
# 1) Performing PCA on meshes or synapse locations to get primary axes
# 2) Projecting a set of synapses or mesh points onto new axes
# 3) Plotting the distribution of synapses over ROI outlins


MeshOutline <- function(Mesh, Plane, alpha_rad=100){
  
  if (Plane=="XY"){
    meshOutline <- ahull(x=Mesh$X, y=-Mesh$Y,alpha=alpha_rad)
  } else if (Plane=="XZ"){
    meshOutline <- ahull(x=Mesh$X, y=Mesh$Z,alpha=alpha_rad)
  } else if (Plane=="YZ"){
    meshOutline <- ahull(x=Mesh$Y, y=Mesh$Z,alpha=alpha_rad)
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


rotatePoints <- function( x, y, z, rot = c(0,0,0), flip = c(1,-1,1) ){
  # rotate points
  points = data.matrix(data.frame(x=x,y=y,z=z))
  points = points %*% makeRotMatXY(rot[1])
  points = points %*% makeRotMatYZ(rot[2])
  points = points %*% makeRotMatXZ(rot[3])
  
  #make df, flip if necessary
  pointsDf = data.frame(x=flip[1]*points[,1],y=flip[2]*points[,2],z=flip[3]*points[,3])
  
  return(pointsDf)
}














  