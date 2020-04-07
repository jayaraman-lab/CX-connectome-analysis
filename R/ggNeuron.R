
neuronForGG <- function(neur){UseMethod("neuronForGG")}

neuronForGG.neuron <- function(neur){
  dList <- lapply(neur$SegList,function(l){filter(neur$d,PointNo %in% l)})
  outD <- bind_rows(dList,.id="Segment") %>% rename(x=X,y=Y,z=Z) %>% mutate(bodyid=neur$bodyid)
  outD
}

neuronForGG.neuronlist <- function(neurL){
  outD <- bind_rows(lapply(neurL,neuronForGG)) %>% mutate(Segment = interaction(Segment,bodyid))
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
                                                    upperY =yPar + dy) %>% arrange(PointNo)
  data.frame(x=c(out$lowerX,rev(out$upperX)),y=c(out$lowerY,rev(out$upperY)))
}

StatNeuron <- ggproto("StatNeuron", Stat,
                    required_aes = c("x","y","Parent","PointNo","W","group"),
                    default_aes = aes(color=NA),
                    compute_group = function(data, scales) {
                       addOutlines(data)
                    }
)

stat_neuron <- function(data = NULL, mapping = aes(Parent=Parent,PointNo=PointNo,W=W,group=Segment), geom = "polygon",
                       position = "identity", na.rm = FALSE, show.legend = NA, 
                       inherit.aes = TRUE, ...) {
  
  layer(
    stat = StatNeuron, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}
