# This script contains unitilies for parsing FB ROI outlines



Get_Mid <- function(Outline){
  
  # Cut the tips off the layer
  Outline_Cut=CutTips_FB(Outline, 0.85)
  
  
  # Get upper and lower half of layer 
  Outline_Cut_UpLow=UpperLower_FB(Outline_Cut)
  
  
  # Get orthogonal mid points
  Outline_MID=BisectShape_Ortho(Outline_Cut_UpLow)
  
  # Extend bisect line to edges of outline
  Outline_MID=Extend_Mid(Outline_MID, Outline)
  
  p1<-ggplot() +  geom_path(data=Outline, aes(x=c1, y=c2), size = 1, color="red") +
    geom_path(data=Outline_Cut, aes(x=c1, y=c2), size = 1, color="black") + 
    geom_path(data=subset(Outline_Cut_UpLow, UpLow=="Upper"), aes(x=c1, y=c2), size = 1, color="blue")  +
    geom_path(data=subset(Outline_Cut_UpLow, UpLow=="Lower"), aes(x=c1, y=c2), size = 1, color="green") +
    geom_path(data=Outline_MID, aes(x=X, y=Y), size = 1, color="orange") + coord_fixed(ratio = 1)
    print(p1)
  
    return(Outline_MID)
}


Extend_Mid <- function(Outline_MID, Outline){
  
  # Fill in left side
  Left_LocalSlope=median(diff(Outline_MID$Y[1:10])/diff(Outline_MID$X[1:10]))
  Left_Intercept=Outline_MID$Y[1] -  Left_LocalSlope*Outline_MID$X[1]
  Left_Outline=subset(Outline, c1< (-3500))
  Left_Outline=Left_Outline[order(Left_Outline$c1),]
  Left_Line = data.frame(X= (Left_Outline$c1), Y= (Left_Outline$c1) *  Left_LocalSlope + Left_Intercept )
  Left_Distance= ((Left_Line$X - Left_Outline$c1)^2 + (Left_Line$Y - Left_Outline$c2)^2)^0.5
  Left_Ind=which(Left_Distance == min(Left_Distance))
  Left_X=Left_Outline$c1[Left_Ind]
  Left_Line=subset(Left_Line, X>=Left_X & X<Outline_MID$X[1]-10)
  
  Outline_MID=rbind(Left_Line,Outline_MID)

  # Fill in Right side
  Right_LocalSlope=median(diff(Outline_MID$Y[ seq( (length(Outline_MID$Y)-10),length(Outline_MID$Y)) ])/diff(Outline_MID$X[seq( (length(Outline_MID$Y)-10),length(Outline_MID$Y))]))
  Right_Intercept=Outline_MID$Y[length(Outline_MID$Y)] -  Right_LocalSlope*Outline_MID$X[length(Outline_MID$Y)]
  Right_Outline=subset(Outline, c1>4000)
  Right_Outline=Right_Outline[order(Right_Outline$c1),]
  Right_Line = data.frame(X= (Right_Outline$c1), Y= (Right_Outline$c1) *  Right_LocalSlope + Right_Intercept )
  Right_Distance= ((Right_Line$X - Right_Outline$c1)^2 + (Right_Line$Y - Right_Outline$c2)^2)^0.5
  Right_Ind=which(Right_Distance == min(Right_Distance))
  Right_X=Right_Outline$c1[Right_Ind]
  Right_Line=subset(Right_Line, X<Right_X & X>Outline_MID$X[length(Outline_MID$X)]+10)
  
  Outline_MID=rbind(Outline_MID,Right_Line)
  
  
  ggplot() +  geom_path(data=Outline, aes(x=c1, y=c2), size = 1, color="red") +
    geom_path(data=Outline_MID, aes(x=X, y=Y), size = 1, color="orange") + coord_fixed(ratio = 1) +
    geom_path(data=Left_Line, aes(x=X,y=Y),color="green")

  
  return(Outline_MID)
}



CutTips_FB <- function(Outline, PRC){
  Distance=abs(Outline$c1) + Outline$c2
  Distance_Thresh=quantile(Distance, probs = PRC) 
  
  Outline_Cut=subset(Outline, (abs(c1)+c2) < Distance_Thresh)
  return(Outline_Cut)
}




UpperLower_FB <- function(Outline){
  Distance= (diff(Outline$c1)^2 +  diff(Outline$c2)^2 )^0.5
  Breaks=which(Distance>600)
  Upper=Outline[seq( (Breaks[1]+1),Breaks[2]), ]
  Upper=Upper[order(Upper$c1),]
  Upper$UpLow="Upper"
  Lower=Outline[c(seq(1,Breaks[1]-1), seq((Breaks[2]+1),length(Outline$c1))) , ]
  Lower=Lower[order(Lower$c1),]
  Lower$UpLow="Lower"
  UpperLower=rbind(Upper, Lower)
  return(UpperLower)
}



BisectShape_Ortho <- function(UpperLower){
  
  # Get upper and lower boundaries
  Upper=subset(UpperLower, UpLow=="Upper")
  Lower=subset(UpperLower, UpLow=="Lower")
  
  # Smooth upper and low lines
  Upper=as.data.frame(rollapply(Upper[c("c1","c2")], 50, mean, partial = TRUE))
  Lower=as.data.frame(rollapply(Lower[c("c1","c2")], 50, mean, partial = TRUE))
  
  Mid=data.frame(X=numeric(), Y=numeric())
  for (ppp in 3:(length(Upper$c1)-3) ){
    
    # compute local slope and intercept
    Local_Xs=Upper$c1[seq(ppp-3,ppp+3)]
    Local_Ys=Upper$c2[seq(ppp-3,ppp+3)]
    Local_Slope = median(diff(Local_Ys)/diff(Local_Xs))
    Local_OrthoSlope=-1/Local_Slope
    Intercept=Upper$c2[ppp]-(Upper$c1[ppp]*Local_OrthoSlope)
    
    # make line that goes through this point
    XXX_Temp=Lower$c1
    Line_Temp=data.frame(X = XXX_Temp , Y = XXX_Temp*Local_OrthoSlope + Intercept)
    
    # Find point of lower edge that intersects the line 
    Distance= ((Line_Temp$X-Lower$c1)^2 + (Line_Temp$Y-Lower$c2)^2)
    IND=which(Distance == min(Distance))
    
    # Get average of two points
    Temp_df=data.frame( X = mean(c(Upper$c1[ppp],Lower$c1[IND])), Y = mean(c(Upper$c2[ppp],Lower$c2[IND])) )
    Mid=rbind(Mid,Temp_df)
    
    # Optional plotting for debugging
    PLOT=0
    if (PLOT==1){
      plot_point=data.frame( X = Upper$c1[ppp], Y =Upper$c2[ppp] )
      plot_point2=data.frame( X = Lower$c1[IND], Y =Lower$c2[IND] )
      ggplot() +  geom_path(data=Upper, aes(x=c1, y=c2), size = 1, color="orange") +
        geom_path(data=Lower, aes(x=c1, y=c2), size = 1, color="red") +
        geom_path(data=Line_Temp, aes(x=X, y=Y), size = 1, color="green") +
        geom_point(data=plot_point, aes(x=X, y=Y)) + geom_point(data=plot_point2, aes(x=X, y=Y)) + ylim(c( min(Lower$c2), max(Upper$c2)))
    }
  }
  
  # Order and smooth
  Mid=Mid[order(Mid$X),]
  Mid=as.data.frame(rollapply(Mid, 5, median, partial = TRUE))
  Mid=as.data.frame(rollapply(Mid, 50, mean, partial = TRUE))
  
  return(Mid)
}
