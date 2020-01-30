# This file contain functions for:
# 1) Performing PCA on meshes or synapse locations to get primary axes
# 2) Projecting a set of synapses or mesh points onto new axes


PCA_SynapseRoi <- function(Input_DF, Rotate) {
  
  
  # Get synapse locations from input data frame, which must have x,y,z columns
  syn_locs=data.frame(x=Input_DF$x, y=Input_DF$y, z=Input_DF$z)
  
  
  # calculate covariance matrix
  synCov = cov(syn_locs)
  
  
  # calculate eigenvectors and values
  synCov_eigen = eigen(synCov) # vectors: columns contain PCs, ranked from most variance accounted for to least 
  
  
  # Rotate the first two PCs to align ROI in direction we want 
  Rotate_Rad=Rotate/360*2*pi
  RotationMatrix=matrix( c(cos(Rotate_Rad), sin(Rotate_Rad), 0, -sin(Rotate_Rad),cos(Rotate_Rad),0,0,0,1) , nrow = 3, ncol = 3)
  synCov_eigen$vectors =  synCov_eigen$vectors %*% RotationMatrix
  
  
  # Use the last eigen vector as the plane to bisect the ROI
  n_plane =c(synCov_eigen$vectors[1,3], synCov_eigen$vectors[2, 3], synCov_eigen$vectors[3, 3])

  
  # Zero mean the data 
  syn_locs_centered = transform(syn_locs, x = x-mean(x), y = y-mean(y), z = z-mean(z))

  
  # Make a matrix from the centered points
  StartMat = matrix(c(syn_locs_centered$x, syn_locs_centered$y, syn_locs_centered$z), nrow = length(syn_locs_centered$x), ncol = 3)
  colnames(StartMat) <- c("x","y","z")
  
  
  # Project the original data onto PCs and create a new data frame with the synapse locations (for plotting)
  NewProjection = StartMat %*% synCov_eigen$vectors
  syn_locs_centered_PCspace=data.frame(x=NewProjection[,1], y=NewProjection[,2], z=NewProjection[,3])
  
  
  # Plot in 3D to check that it captures the axes we want
  #nclear3d()
  plot3d(syn_locs_centered[c('x','y','z')], col='seagreen3',alpha=0.1)
  vectors3d(sqrt(synCov_eigen$values[1])*synCov_eigen$vectors[1:3, 1], col="black", lwd=2, radius=1/25) # Largest PC
  vectors3d(sqrt(synCov_eigen$values[2])*synCov_eigen$vectors[1:3, 2], col="gray", lwd=2, radius=1/25) # Second largest PC
  vectors3d(sqrt(synCov_eigen$values[3])*synCov_eigen$vectors[1:3, 3], col="red", lwd=2, radius=1/25)   # Smallest PC
  planes3d(n_plane[1],n_plane[2],n_plane[3], 0, col="gray", alpha=0.2)
  decorate3d(box=FALSE)

  # Make scatter plot of two projections
  p1<-ggplot(syn_locs_centered_PCspace, aes(x=x, y=y)) + geom_point(size=1, alpha = 0.05) +  xlab("x=PC1") + ylab("y=PC2")
  print(p1)
  
  p2<-ggplot(syn_locs_centered_PCspace, aes(x=x, y=z)) + geom_point(size=1, alpha = 0.05) +  xlab("x=PC1") + ylab("z=PC3")
  print(p2)
  
  p3<-ggplot(syn_locs_centered_PCspace, aes(x=y, y=z)) + geom_point(size=1, alpha = 0.05) +  xlab("y=PC2") + ylab("z=PC3")
  print(p3)
  
  return(synCov_eigen)
}


Project_NewCoordinates <- function(Input_DF, NewAxes) {
  
  
  # Make a matrix from the synapse or mesh points
  syn_locs=data.frame(x=Input_DF$x, y=Input_DF$y, z=Input_DF$z)
  StartMat = matrix(c(syn_locs$x, syn_locs$y, syn_locs$z), nrow = length(syn_locs$x), ncol = 3)
  colnames(StartMat) <- c("x","y","z")
  
  
  # Project the original data onto PCs and create a new data frame with the synapse locations (for plotting)
  NewProjection = StartMat %*% NewAxes$vectors
  NewProjection_DF = Input_DF
  NewProjection_DF$x=NewProjection[,1]
  NewProjection_DF$y=NewProjection[,2]
  NewProjection_DF$z=NewProjection[,3]


  return(NewProjection_DF)
  
}
  
  
  
  
  
  
  
  
  
  
  
  
  

