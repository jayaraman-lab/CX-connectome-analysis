# This file contain functions for:
# 0) A number of modular functions from 1-4 for performing many PCA/projection related manipulations of synapse and roi locations.
# 1) Performing PCA on meshes or synapse locations to get primary axes
# 2) Projecting a set of synapses or mesh points onto new axes
# 3) Plotting the distribution of synapses over ROI outlins
# 4) Getting synapse locations in a loop so the query doesnt time out
# Note, eventually, we should remove the large functions and work off the modular ones near the top


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
    
    # Make sure last block is not length 1
    if (tail(Starts,n=1)== tail(Stops,n=1)){
      Starts=Starts[1:length(Starts)-1]
      Stops=Stops[1:length(Stops)-1]
      Stops[length(Stops)]<-length(BodyIDs)
    }
    
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







  
  
  
  
PlotSynapseDistributions <- function(Synapses, ROI_Outline1, ROI_Outline2, ROI_Outline3, ROI, Cmax, PlotDir) { 
  
  
  # Get max and min synapse positions
  Extra=500
  ALLDATA=rbind(Synapses[,3:5],ROI_Outline1[,1:3],ROI_Outline2[,1:3])
  XLIM=c( min(c(ALLDATA$x))-Extra , max(c(ALLDATA$x))+Extra )
  YLIM=c( min(c(-ALLDATA$y))-Extra , max(c(-ALLDATA$y))+Extra ) # take negative here
  ZLIM=c( min(c(ALLDATA$z))-Extra , max(c(-ALLDATA$z))+Extra )
  
  
  # Plot all data to start
  ppp=ggplot() + geom_polygon(data=ROI_Outline1, aes(x = x, y = -y), colour='black', fill=NA, size = 0.75 ) +
    geom_polygon(data=ROI_Outline2, aes(x = x, y = -y), colour='black', fill=NA, size = 0.75 ) +
    geom_hex(data=Synapses, aes(x = x, y = -y), bins=50 ) + 
    scale_fill_gradientn(colors = brewer.pal(3,"Blues"), breaks=c(0,Cmax), limits=c(0, Cmax), oob = scales::squish) +
    xlim(XLIM[1], XLIM[2]) + ylim(YLIM[1], YLIM[2]) + ggtitle("ALL") +
    coord_fixed()   + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                            panel.background = element_blank(), axis.line = element_line(colour = "black"))
  ggsave(paste(PlotDir, ROI, "_RingNeuronSynDistribution_","ALL",".pdf",sep=""),
         plot = ppp, device='pdf', scale = 1, width = 8, height = 7, units ="in", dpi = 600, limitsize = TRUE)
  remove(ppp)
  
  
  # Make same plot for each neuron type now
  PlotHandles1 <- vector("list", length(unique(Synapses$type)))
  PlotHandles2 <- vector("list", length(unique(Synapses$type)))
  # Get 95% contour and plot all on one plot
  for (nnn in 1:length(unique(Synapses$type)) ){
    
    # Get the synapses to plot
    ToPlotSyns=subset(Synapses, Synapses$type == sort(unique(Synapses$type))[nnn] )
    
    # Plot the 2D hisogram without contour lines
    ppp=ggplot()  + geom_polygon(data=ROI_Outline1, aes(x = x, y = -y), colour='black', fill=NA, size = 0.75 ) +
      geom_polygon(data=ROI_Outline2, aes(x = x, y = -y), colour='black', fill=NA, size = 0.75 ) +
      geom_polygon(data=ROI_Outline3, aes(x = x, y = -y), colour='black', fill=NA, size = 0.75 ) +
      geom_hex(data=ToPlotSyns, aes(x = x, y = -y), bins=50 ) + 
      scale_fill_gradientn(colors = brewer.pal(3,"Blues"), breaks=c(0,Cmax), limits=c(0, Cmax), oob = scales::squish) +
      xlim(XLIM[1], XLIM[2]) + ylim(YLIM[1], YLIM[2]) + ggtitle(sort(unique(Synapses$type))[nnn]) +
      coord_fixed() + theme_void()  + theme(legend.position="none")
      
    
    #ggsave(paste(PlotDir, ROI, "_RingNeuronSynDistribution_",sort(unique(Synapses$type))[nnn],".png",sep=""),
           #plot = ppp, device='png', scale = 1, width = 8, height = 7, units ="in", dpi = 600, limitsize = TRUE)
    PlotHandles1[[nnn]] <- ppp
    remove(ppp)
    
    
    ContourType=1
    if (ContourType==1){
      # Get synapse 95% contour, but this seems to only work well for gaussian-like distributions
      d <- data.frame(x=ToPlotSyns$x,y=ToPlotSyns$y)
      kd <- ks::kde(d, compute.cont=TRUE, bgridsize=c(7,7))
      contour_95 <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]], z=estimate, levels=cont["2%"])[[1]])
      contour_95 <- data.frame(contour_95)
      contour_95$name=sort(unique(Synapses$type))[nnn]
      
      # Smooth contour 
      Contour=data.frame(x = rollmean(contour_95$x,4),
                         y = rollmean(contour_95$y,4),
                         name = sort(unique(Synapses$type))[nnn])
    } else if (ContourType==2){
      
      plot_data=ToPlotSyns[,c(3,4,9)]
      colnames(plot_data) <- c("X","Y","Label")
      plot_data$Label=as.factor(plot_data$Label)
      
      #Genrate the kernel density for each group
      newplot_data <- plot_data %>% group_by(Label) %>% do(Dens=kde2d(.$X, .$Y, h=1 , n=100, lims=c(XLIM,-YLIM)))
      
      #Transform the density in  data.frame
      newplot_data  %<>%  do(Label=.$Label, V=expand.grid(.$Dens$x,.$Dens$y), Value=c(.$Dens$z)) %>% do(data.frame(Label=.$Label,x=.$V$Var1, y=.$V$Var2, Value=.$Value))
      
      #Thresh=quantile(newplot_data$Value, c(.10)) 
      #newplot_data$Value[newplot_data$Value>Thresh]<-1
      #newplot_data$Value[newplot_data$Value<=Thresh]<-0
      
      ggplot()  + geom_polygon(data=ROI_Outline1, aes(x = x, y = -y), colour='black', fill=NA, size = 0.75 ) +
        geom_polygon(data=ROI_Outline2, aes(x = x, y = -y), colour='black', fill=NA, size = 0.75 ) +
        stat_contour(data=newplot_data, aes(x=x, y=-y, z=Value, fill=Label), geom="polygon", alpha=1) 
      
    }
    
    
    # Plot the 2D hisogram with contour lines
    ppp = ggplot() + geom_polygon(data=ROI_Outline1, aes(x = x, y = -y), colour='black', fill=NA, size = 0.75 ) +
      geom_polygon(data=ROI_Outline2, aes(x = x, y = -y), colour='black', fill=NA, size = 0.75 ) +
      geom_polygon(data=ROI_Outline3, aes(x = x, y = -y), colour='black', fill=NA, size = 0.75 ) +
      geom_hex(data=ToPlotSyns, aes(x = x, y = -y), bins=50 ) + 
      geom_polygon(data=Contour, aes(x=x, y=-y), colour='red', fill=NA, size = 0.5) +
      scale_fill_gradientn(colors = brewer.pal(3,"Blues"), breaks=c(0,Cmax), limits=c(0, Cmax), oob = scales::squish) +
      xlim(XLIM[1], XLIM[2]) + ylim(YLIM[1], YLIM[2]) + ggtitle(sort(unique(Synapses$type))[nnn]) +
      coord_fixed() + theme_void()  + theme(legend.position="none")
    
    
    #ggsave(paste(PlotDir, ROI, "_RingNeuronSynDistribution_Contour_",sort(unique(Synapses$type))[nnn],".png",sep=""),
           #plot = ppp, device='png', scale = 1, width = 8, height = 7, units ="in", dpi = 600, limitsize = TRUE)
    PlotHandles2[[nnn]] <- ppp
    remove(ppp)
    
    if (nnn==1){
      ContourAll=Contour
    } else {
      ContourAll=rbind(ContourAll,Contour)
    }
    
  }
  
  ppp=ggplot() + geom_polygon(data=ROI_Outline1, aes(x = x, y = -y), colour='black', fill=NA, size = 0.75 ) +
    geom_polygon(data=ROI_Outline2, aes(x = x, y = -y), colour='black', fill=NA, size = 0.75 ) +
    geom_polygon(data=ContourAll, aes(x=x, y=-y,color=name, fill = name), alpha=0.05) + theme_bw()
  ggsave(paste(PlotDir, ROI, "_RingNeuronSynDistribution",".png",sep=""),
         plot = ppp, device='png', scale = 1, width = 8, height = 7, units ="in", dpi = 600, limitsize = TRUE)
  
  
  
  Output= list("PlotHandles1"=PlotHandles1,"PlotHandles2"=PlotHandles2)
  return(Output)
  
  
}

  
  



  
  

























