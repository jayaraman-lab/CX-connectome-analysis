Get_SynLayerDistribution <- function(All_Neurons, Thresh, PlotDir){
  #' Function for getting the width of each neuron's synapse distribution in a direction tangent to the FB layer
  
  # Initialize dataframe that will contain all information about each neuron's synapse distribution
  Syn_Distribution <- data.frame(bodyid=as.numeric(), type=as.character(), Layer=as.character(), 
                                 Side=as.character(), PBglom=as.character(), FBcol=as.character(),
                                 X_Mean=as.numeric(), Y_Mean=as.numeric(), Z_Mean=as.numeric(), 
                                 X_STD_Orth=as.numeric(), Z_STD_Orth=as.numeric(),
                                 X_HW_Orth=as.numeric(), Z_HW_Orth=as.numeric(),
                                 Cross_XL_X=as.numeric(),  Cross_XL_Z=as.numeric(),
                                 Cross_XR_X=as.numeric(),  Cross_XR_Z=as.numeric(),
                                 Cross_YU_X=as.numeric(),  Cross_YU_Z=as.numeric(),
                                 Cross_YD_X=as.numeric(),  Cross_YD_Z=as.numeric(),
                                 Box_UL_X=as.numeric(),    Box_UL_Z=as.numeric(),
                                 Box_UR_X=as.numeric(),    Box_UR_Z=as.numeric(),
                                 Box_LL_X=as.numeric(),    Box_LL_Z=as.numeric(),
                                 Box_LR_X=as.numeric(),    Box_LR_Z=as.numeric(),
                                 LayerLength=as.numeric(), NumberOfSyns=as.numeric())
  
  # Get unique neuron types and layers
  Neuron_Types=unique(All_Neurons$type)
  Layers_All=sort(unique(All_Neurons$Layer))
  
  # Loop over neuron types
  for (ttt in 1:length(Neuron_Types)){
    SubData=subset(All_Neurons, type == Neuron_Types[ttt])
  
    # Loop over layers
    for (lll in 1:8){
      LayerData=subset(SubData, Layer==Layers_All[lll])
      
      # Get the outline and bisecting line for this layer
      eval(parse( text= paste("MID=Mid_L", as.character(lll),sep="") ))
      eval(parse( text= paste("OUTLINE=Outline_L", as.character(lll),"_XZ",sep="") ))
      
      # Loop over individual neurons
      NeuronIDs=unique(LayerData$bodyid)
      for (nnn in 1:length(NeuronIDs)){
        NeuronData=subset(LayerData, bodyid==NeuronIDs[nnn])
        
        # Rename "Z" as "Y" to match convention that y-axis is vertical in plots. This corresponds the the anterior-posterior axis.
        Mean_YY=mean(NeuronData$Y)                               # Get average dorsal-ventral position. 
        NeuronData=NeuronData[colnames(NeuronData) != "Y"]       # Get rid of dorsal-ventral information (i.e. "Y")
        colnames(NeuronData)[colnames(NeuronData) == "Z"] <-"Y"  # Rename Z as Y, the new anterior-posterior variable
        
        # Get X range (width along the medial-lateral axis)
        X_RangeMin = quantile(NeuronData$X, probs = 0.1) 
        X_RangeMax = quantile(NeuronData$X, probs = 0.9) 
        
        # Check if there are enough synapses
        if (length(NeuronData$X)>Thresh){
          
        #####################################################################################################################  
        ############# Project data into new coordinate frame whose main axis is locally tangent to FB the layer ############# 
          
          # Compute the length of this layer
          Layer_Length_Temp=sum((diff(MID$X)^2 + diff(MID$Y)^2)^0.5)
          
          # Get a segment of the bisecting line (i.e. MID) that overlaps with the synapse cloud.
          # If there are no parts of bisecting line within the synapse cloud, get the closest segment.
          Mid_subset=subset(MID, X>X_RangeMin & X<X_RangeMax)
          if (length(Mid_subset$X)<2){
            Center=mean(c(X_RangeMin,X_RangeMax))
            INDS= order(abs(MID$X - Center))[c(1,2)]
            Mid_subset=MID[INDS,]
          }
          
          # Compute the median slope of this segment
          TempSlope=median(diff(Mid_subset$Y)/diff(Mid_subset$X))
          
          # Define new axis (basis vectors should be length 1 and orthogonal to keep the shape of the synapse cloud)
          # Also, keep Ax_1 pointing towards right and Ax_2 pointing up, so that axes are consistent across clouds
          Ax_1=c(1,TempSlope)
          Ax_1=Ax_1/sum(Ax_1^2)^0.5
          
          Ax_2=c(1, -Ax_1[1]/Ax_1[2])
          Ax_2=Ax_2/sum(Ax_2^2)^0.5
          if (Ax_2[2]<0){
            Ax_2[1]=-Ax_2[1]
            Ax_2[2]=-Ax_2[2]} 

          # Make sure vectors are orthogonal
          if (Ax_1 %*% Ax_2 > 10^-15){print(paste("Vectors not orthogonal ttt=", as.character(ttt), "   lll = ", as.character(lll), "   nnn = ", as.character(nnn)))}
          
          # Project synapse locations into new coordinate frame
          NeuronData$X_New= as.matrix(NeuronData[c("X","Y")]) %*% Ax_1  
          NeuronData$Y_New= as.matrix(NeuronData[c("X","Y")]) %*% Ax_2  
          
          # Calculate vectors to get from new coordinate frame back to original coordinate frame
          Ax1_NewToOld= t(Ax_1 * as.matrix(c(1,-1)))
          Ax2_NewToOld= t(Ax_2 * as.matrix(c(-1,1)))
          
        #####################################################################################################################   
        ############# Compute summary statistics about the shape of the synapse cloud  ###################################### 
          
          # Calculate mean, STD, and half widths along new axes
          TempMean_New=as.data.frame(t(colMedians(as.matrix(NeuronData[c("X_New","Y_New")]))))
          TempSTD_New=as.data.frame(t(sapply(NeuronData[c("X_New","Y_New")], sd)))
          TempRange_New=data.frame(X = abs(quantile(NeuronData$X_New, probs = 0.05) - quantile(NeuronData$X_New, probs = 0.95))/2,
                                   Y = abs(quantile(NeuronData$Y_New, probs = 0.05) - quantile(NeuronData$Y_New, probs = 0.95))/2)
          
          # Calculate mean location in original coordinate frame 
          TempMean=data.frame(X = as.matrix(TempMean_New[c("X_New","Y_New")]) %*% t(Ax1_NewToOld), Y = as.matrix(TempMean_New[c("X_New","Y_New")]) %*% t(Ax2_NewToOld))
          
          # Compute Cross points in original coordinate frame (4 points making up a cross around mean synapse position)
          Cross_XL= data.frame(X = as.numeric(as.matrix((TempMean_New -  c(TempRange_New$X,0))) %*% t(Ax1_NewToOld)) , 
                               Y = as.numeric(as.matrix((TempMean_New -  c(TempRange_New$X,0))) %*% t(Ax2_NewToOld)) )
          Cross_XR= data.frame(X = as.numeric(as.matrix((TempMean_New +  c(TempRange_New$X,0))) %*% t(Ax1_NewToOld)) , 
                               Y = as.numeric(as.matrix((TempMean_New +  c(TempRange_New$X,0))) %*% t(Ax2_NewToOld)) )
          Cross_YU= data.frame(X = as.numeric(as.matrix((TempMean_New +  c(0,TempRange_New$Y))) %*% t(Ax1_NewToOld)) , 
                               Y = as.numeric(as.matrix((TempMean_New +  c(0,TempRange_New$Y))) %*% t(Ax2_NewToOld)) )
          Cross_YD= data.frame(X = as.numeric(as.matrix((TempMean_New -  c(0,TempRange_New$Y))) %*% t(Ax1_NewToOld)) , 
                               Y = as.numeric(as.matrix((TempMean_New -  c(0,TempRange_New$Y))) %*% t(Ax2_NewToOld)) )
          
          # Compute box points in original coordinate frame (4 points making up a box around mean synapse position)
          Box_UL= data.frame(X = as.numeric(as.matrix((TempMean_New +  TempRange_New*c(-1,1) )) %*% t(Ax1_NewToOld)) , 
                             Y = as.numeric(as.matrix((TempMean_New +  TempRange_New*c(-1,1) )) %*% t(Ax2_NewToOld)) )
          Box_UR= data.frame(X = as.numeric(as.matrix((TempMean_New +  TempRange_New*c(1,1) )) %*% t(Ax1_NewToOld)) , 
                             Y = as.numeric(as.matrix((TempMean_New +  TempRange_New*c(1,1) )) %*% t(Ax2_NewToOld)) )
          Box_LL= data.frame(X = as.numeric(as.matrix((TempMean_New +  TempRange_New*c(-1,-1) )) %*% t(Ax1_NewToOld)) , 
                             Y = as.numeric(as.matrix((TempMean_New +  TempRange_New*c(-1,-1) )) %*% t(Ax2_NewToOld)) )
          Box_LR= data.frame(X = as.numeric(as.matrix((TempMean_New +  TempRange_New*c(1,-1) )) %*% t(Ax1_NewToOld)) , 
                             Y = as.numeric(as.matrix((TempMean_New +  TempRange_New*c(1,-1) )) %*% t(Ax2_NewToOld)) )
          
        ######################################################################################################################   
        ############# Plot example data to make sure the synapse clouds are being measured accurately ########################
        
          if ( ttt*lll*nnn %% 50) ==0 ){
            p1<-ggplot() +  geom_path(data=OUTLINE, aes(x=c1, y=c2), size = 1, color="red") +
              geom_point(data=NeuronData, aes(x=X,y=Y)) + 
              geom_path(data=MID, aes(x=X, y=Y), size = 1, color="orange") + 
              geom_point(data=TempMean, aes(x=X,y=Y,size=50),color="red") + 
              geom_point(data=Cross_XL, aes(x=X,y=Y,size=50),color="cornflowerblue") + 
              geom_point(data=Cross_XR, aes(x=X,y=Y,size=50),color="blueviolet") + 
              geom_point(data=Cross_YU, aes(x=X,y=Y,size=50),color="darkred") + 
              geom_point(data=Cross_YD, aes(x=X,y=Y,size=50),color="gold") + 
              geom_point(data=Box_UL, aes(x=X,y=Y,size=50),color="darkslateblue")  + 
              geom_point(data=Box_UR, aes(x=X,y=Y,size=50),color="darkorange") + 
              geom_point(data=Box_LL, aes(x=X,y=Y,size=50),color="deeppink") +
              geom_point(data=Box_LR, aes(x=X,y=Y,size=50),color="darkslategray") + coord_fixed(ratio = 1) +
              ggtitle(paste(NeuronData$type[1],"_Layer_",as.character(lll),sep=""))
            
            ggsave(paste(PlotDir, "Cloud_", NeuronData$type[1],"_Layer_",as.character(lll),"_", as.character(ttt*lll*nnn),".png",sep=""),
                   plot = p1, device='png', scale = 1, width = 8, height = 5, units ="in", dpi = 150, limitsize = TRUE)
            remove(p1)
          }
          
        ##################################################################################################################### 
        ############# Add data to dataframe ################################################################################# 
          
          TempDF <- data.frame(bodyid=NeuronData$bodyid[1], type=NeuronData$type[1], Layer=NeuronData$Layer[1], 
                               Side=NeuronData$Side[1], PBglom=NeuronData$PBglom[1], FBcol=NeuronData$FBcol[1],
                                         X_Mean=TempMean$X, Y_Mean=Mean_YY, Z_Mean=TempMean$Y, 
                                         X_STD_Orth=TempSTD_New$X_New, Z_STD_Orth=TempSTD_New$Y_New,
                                         X_HW_Orth=TempRange_New$X, Z_HW_Orth=TempRange_New$Y,
                                         Cross_XL_X=Cross_XL$X,    Cross_XL_Z=Cross_XL$Y,
                                         Cross_XR_X=Cross_XR$X,    Cross_XR_Z=Cross_XR$Y,
                                         Cross_YU_X=Cross_YU$X,    Cross_YU_Z=Cross_YU$Y,
                                         Cross_YD_X=Cross_YD$X,    Cross_YD_Z=Cross_YD$Y,
                                         Box_UL_X=Box_UL$X,        Box_UL_Z=Box_UL$Y,
                                         Box_UR_X=Box_UR$X,        Box_UR_Z=Box_UR$Y,
                                         Box_LL_X=Box_LL$X,        Box_LL_Z=Box_LL$Y,
                                         Box_LR_X=Box_LR$X,        Box_LR_Z=Box_LR$Y,
                                         LayerLength = sum((diff(MID$X)^2 + diff(MID$Y)^2)^0.5),
                                         NumberOfSyns = length(NeuronData$X) )

          Syn_Distribution=rbind(Syn_Distribution,TempDF)
          
        }
      }
    }
  }
  return(Syn_Distribution)
}

