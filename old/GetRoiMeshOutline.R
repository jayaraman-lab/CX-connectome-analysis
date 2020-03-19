# This file contain a functions getting the outline of an mesh

Get_MeshRoiOutline <- function(DataToUse, Spacing, Smooth, ROI, DIM) {
  

  # Get X and Y limits for image grid
  XLIM=c( floor(min(DataToUse[,DIM[1]]/Spacing))*Spacing-Spacing*10 , ceiling(max(DataToUse[,DIM[1]]/Spacing))*Spacing+Spacing*10  )
  YLIM=c( floor(min(DataToUse[,DIM[2]])/Spacing)*Spacing-Spacing*10  , ceiling(max(DataToUse[,DIM[2]])/Spacing)*Spacing+Spacing*10  ) 
  
  
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
  for (val in 1:length(DataToUse$x)){ImageMat[ as.character(DataToUse[,DIM[2]][val]),
                                               as.character(DataToUse[,DIM[1]][val]) ] = 1}
  
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
  Edges=data.frame(x=XSpacing[Cimage_Edges_df$x], y=YSpacing[Cimage_Edges_df$y], z=0, name="Edge")
  
  # Loop over data and get ordered points by finding nearest neighbor 
  EdgesTemp=Edges
  for (val in 1:length(EdgesTemp$x)){
    if (val==1){
      EdgesOrdered=data.frame(x=EdgesTemp$x[val], y=EdgesTemp$y[val], z=EdgesTemp$z[val], name='ROI')
      EdgesTemp[val,]=NA
    } else {
      Distances = ((EdgesTemp$x-EdgesOrdered$x[val-1])^2 + (EdgesTemp$y-EdgesOrdered$y[val-1])^2)^0.5
      Ind=which.min(Distances)
      EdgesOrdered=rbind(EdgesOrdered, data.frame(x=EdgesTemp$x[Ind], y=EdgesTemp$y[Ind], z=EdgesTemp$z[Ind], name='ROI') )
      EdgesTemp[Ind,]=NA
    }
  }
  EdgesOrdered[(length(EdgesTemp$x)+1),]=EdgesOrdered[1,]
  
  # Get distances to see if there are two edges or one
  DistanceAll= ((diff(EdgesOrdered$x))^2 + (diff(EdgesOrdered$y))^2)^0.5
  DistanceAll[length(DistanceAll)]=0
  
  # Rearrange axes
  
  
  # If there are two edges, divide into two. Either way, plot the data. 
  DistanceThresh=200
  if (sum(DistanceAll>DistanceThresh)==0){
    EdgesOrdered1=EdgesOrdered
    EdgesOrdered2=NA
    
    # Smooth ROI outline
    EdgesOrdered1_OUT=data.frame(x = rollmean(EdgesOrdered1$x,8),
                                 y = rollmean(EdgesOrdered1$y,8),
                                 z = rollmean(EdgesOrdered1$z,8),
                                 name=ROI)
    EdgesOrdered2_OUT=NA
    
    
    
    
  } else {
    # Divide data into two edges
    IND=which(DistanceAll>DistanceThresh)
    EdgesOrdered1=EdgesOrdered[1:IND ,]
    EdgesOrdered2=EdgesOrdered[((IND+1):(length(EdgesOrdered$x)-1)) ,]

    # Smooth ROI outline
    EdgesOrdered1_OUT=data.frame(x = rollmean(EdgesOrdered1$x,8),
                                 y = rollmean(EdgesOrdered1$y,8),
                                 z = rollmean(EdgesOrdered1$z,8),
                                 name=ROI)
    EdgesOrdered2_OUT=data.frame(x = rollmean(EdgesOrdered2$x,8),
                                 y = rollmean(EdgesOrdered2$y,8),
                                 z = rollmean(EdgesOrdered2$z,8),
                                 name=ROI)
   
  }
  
  
  
  Output= list("EdgesOrdered"=EdgesOrdered,"EdgesOrdered1_OUT"=EdgesOrdered1_OUT,"EdgesOrdered2_OUT"=EdgesOrdered2_OUT)
  #return(Output)
  
}

