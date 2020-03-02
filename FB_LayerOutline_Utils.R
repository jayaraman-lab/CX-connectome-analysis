# This script contains unitilies for parsing FB ROI outlines




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
  Upper=as.data.frame(rollapply(Upper[c("c1","c2")], 25, mean, partial = TRUE))
  Lower=as.data.frame(rollapply(Lower[c("c1","c2")], 25, mean, partial = TRUE))
  
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
