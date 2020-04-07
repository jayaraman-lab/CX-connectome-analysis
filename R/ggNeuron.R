
neuronForGG <- function(neur,axis=c("x","y")){UseMethod("neuronForGG")}

neuronForGG.neuron <- function(neur,axis=c("x","y")){
  refT <- mutate(neur$d,x=!!as.name(toupper(axis[1])),y=!!as.name(toupper(axis[2])))
  refT <- addOutlines(refT)
  ptsOrd <- treeOrder(neur$SegList[[1]],neur$SegList,branchpoints(neur),refTable = refT)
  refT <- left_join(ptsOrd,refT,by = c("idx" = "PointNo"))
  refT <- mutate(refT,x=ifelse(side=="upper",upperX,lowerX),y=ifelse(side=="upper",upperY,lowerY)) %>% mutate(bodyid=neur$bodyid)
  refT %>% select(x,y,bodyid)
}

neuronForGG.neuronlist <- function(neurL,axis=c("x","y")){
  outD <- bind_rows(lapply(neurL,neuronForGG,axis=axis))
}

angle <- function(x, y) {
  atan2(y[2] - y[1], x[2] - x[1])
}

## x and y are vectors of length 2
perpUp <- function(x, y, len, a) {
  dx <- len*cos(a + pi/2)
  dy <- len*sin(a + pi/2)
  upper <- c(x[1] + dx, y[1] + dy)
}

perpLow <-  function(x, y, len, a) {
  dx <- len*cos(a + pi/2)
  dy <- len*sin(a + pi/2)
  lower <- c(x[1] - dx, y[1] - dy)
}
    

## x and y are vectors of length 2
perpStart <- function(x, y, len){
  perp(x, y, len, angle(x, y), 1)
}

addOutlines <- function(neuronD){
  out <- neuronD %>%
    mutate(xPar = ifelse(Parent %in% PointNo,neuronD$x[match(Parent,neuronD$PointNo)],neuronD$x[match(PointNo,neuronD$Parent)]),
           yPar = ifelse(Parent %in% PointNo,neuronD$y[match(Parent,neuronD$PointNo)],neuronD$y[match(PointNo,neuronD$Parent)]),
          angle = atan2(y-yPar,x-xPar),
          dx = W/2 * (cos(angle + pi/2)),
          dy = W/2 * (sin(angle + pi/2)),
          lowerX = xPar - dx,
          lowerY = yPar - dy,
          upperX = xPar +dx,
          upperY =yPar + dy) 
}

treeOrder <- function(subTree,segmentList,branchPts,refTable){
  if (tail(subTree,1) %in% branchPts){
    bP <- tail(subTree,1)
    childrenTrees <- segmentList[sapply(segmentList,function(l) l[1]==bP)]
    childrenSec <- sapply(childrenTrees,function(s) s[2])
    childrenTab <- arrange(filter(refTable,PointNo %in% childrenSec),desc(upperY))$PointNo
    childrenTrees <- childrenTrees[match(childrenTab,childrenSec)]
    return(bind_rows(list(data.frame(idx=subTree,side="upper",stringsAsFactors = FALSE),
                          bind_rows(lapply(childrenTrees,treeOrder,segmentList,branchPts,refTable)),
                          data.frame(idx=rev(subTree),side="lower",stringsAsFactors = FALSE))))
    
  }else{
    return(rbind(data.frame(idx=subTree,side="upper",stringsAsFactors = FALSE),data.frame(idx=rev(subTree),side="lower",stringsAsFactors = FALSE)))
  }
}


