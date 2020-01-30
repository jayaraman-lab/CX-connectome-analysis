# This file contain a functions getting the outline of an mesh

Get_MeshRoiOutline <- function(Mesh, Spacing, Smooth, ROI) {
  

  # Get 3D data points of mesh nodes
  DataToUse=data.frame(x=Mesh$vb[1,], y=Mesh$vb[2,],z=Mesh$vb[3,],name=ROI)
  
  
  # Get X and Y limits for image grid
  XLIM=c( floor(min(DataToUse$x/Spacing))*Spacing-Spacing*10 , ceiling(max(DataToUse$x/Spacing))*Spacing+Spacing*10  )
  YLIM=c( floor(min(DataToUse$z)/Spacing)*Spacing-Spacing*10  , ceiling(max(DataToUse$z)/Spacing)*Spacing+Spacing*10  ) 
  
  
  # Initialize matrix
  XSpacing = seq(from = XLIM[1], to = XLIM[2], by = Spacing)
  YSpacing = seq(from = YLIM[1], to = YLIM[2], by = Spacing)
  ImageMat = matrix(0, nrow = length(YSpacing) , ncol = length(XSpacing))
  rownames(ImageMat) = as.character(YSpacing)
  colnames(ImageMat) = as.character(XSpacing)
  
  
  # Populate image matrix
  DataToUse$x=round(DataToUse$x/Spacing)*Spacing
  DataToUse$y=round(DataToUse$y/Spacing)*Spacing
  DataToUse$z=round(DataToUse$z/Spacing)*Spacing
  for (val in 1:length(DataToUse$x)){ImageMat[ as.character(DataToUse$z[val]),
                                               as.character(DataToUse$x[val]) ] = 1}
  
  # Plot matrix 
  longData<-melt(ImageMat)
  p0<-ggplot(longData, aes(x = Var2, y = -Var1)) +  geom_raster(aes(fill=value))
  print(p0)
  
  # blur image and detect edges
  Cimage=as.cimg(t(ImageMat))
  p1<-plot(Cimage )
  print(p1)
  Cimage_blurr=isoblur(Cimage,Smooth) #3 is good
  p2<-plot(Cimage_blurr)
  print(p2)
  Cimage_Edges=cannyEdges(Cimage_blurr, alpha = 1)
  p3<-plot(Cimage_Edges)
  print(p3)
  
  # Convert back to dataframe and get coordinates
  Cimage_Edges_df <- as.data.frame(Cimage_Edges)
  Edges=data.frame(x=XSpacing[Cimage_Edges_df$x], y=0,z=YSpacing[Cimage_Edges_df$y], name="Edge")
  
  # Loop over data and get ordered points by finding nearest neighbor 
  EdgesTemp=Edges
  for (val in 1:length(EdgesTemp$x)){
    if (val==1){
      EdgesOrdered=data.frame(x=EdgesTemp$x[val], y=EdgesTemp$y[val], z=EdgesTemp$z[val], name='ROI')
      EdgesTemp[val,]=NA
    } else {
      Distances = ((EdgesTemp$x-EdgesOrdered$x[val-1])^2 + (EdgesTemp$z-EdgesOrdered$z[val-1])^2)^0.5
      Ind=which.min(Distances)
      EdgesOrdered=rbind(EdgesOrdered, data.frame(x=EdgesTemp$x[Ind], y=EdgesTemp$y[Ind], z=EdgesTemp$z[Ind], name='ROI') )
      EdgesTemp[Ind,]=NA
    }
  }
  EdgesOrdered[(length(EdgesTemp$x)+1),]=EdgesOrdered[1,]
  
  # Get distances to see if there are two edges or one
  DistanceAll= ((diff(EdgesOrdered$x))^2 + (diff(EdgesOrdered$z))^2)^0.5
  DistanceAll[length(DistanceAll)]=0
  
  # If there are two edges, divide into two. Either way, plot the data. 
  DistanceThresh=200
  if (sum(DistanceAll>DistanceThresh)==0){
    EdgesOrdered1=NA
    EdgesOrdered2=NA
  } else {
    # Divide data into two edges
    IND=which(DistanceAll>DistanceThresh)
    EdgesOrdered1=EdgesOrdered[1:IND ,]
    EdgesOrdered2=EdgesOrdered[((IND+1):(length(EdgesOrdered$x)-1)) ,]
  }
  
  Output= list("EdgesOrdered"=EdgesOrdered,"EdgesOrdered1"=EdgesOrdered1,"EdgesOrdered2"=EdgesOrdered2)
  #return(Output)
  
}

