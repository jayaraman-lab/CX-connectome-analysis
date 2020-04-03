neuronForGG <- function(neur){UseMethod("neuronForGG")}

neuronForGG.neuron <- function(neur){
  neur$d <- mutate(neur$d,segment=0)
  for (i in 1:length(neur$SegList)){
    newD <- mutate(neur$d,segment=ifelse(PointNo %in% neur$SegList[[i]],i,segment))
    neur$d <- distinct(bind_rows(neur$d,newD))
  }
  neur$d <- filter(neur$d,segment!=0) %>% rename(x=X,y=Y,z=Z) %>% mutate(bodyid=neur$bodyid)
  neur$d
}

neuronForGG.neuronlist <- function(neurL){
  bind_rows(lapply(neurL,neuronForGG))
}


neuronForGGMem <- memoise::memoise(neuronForGG)

ggNeuron <- function(neurGG,proj=c("x","y"),...){UseMethod("ggNeuron")}

ggNeuron.neuronlist <- function(neurGG,proj=c("x","y"),...){
  neurGG <- neuronForGGMem(neurGG)
  ggNeuron(neurGG,proj=proj,...)
}

ggNeuron.neuron <- function(neurGG,proj=c("x","y"),...){
  neurGG <- neuronForGGMem(neurGG)
  ggNeuron(neurGG,proj=proj,...)
}

ggNeuron.data.frame <- function(neurGG,proj=c("x","y"),...){
  
 geom_path(data=neurGG,aes(x=!!as.name(proj[1]),y=!!as.name(proj[2]),group=interaction(segment,bodyid),color=bodyid),...)
}

ggSoma <- function(neurGG,proj=c("x","y"),...){UseMethod("ggSoma")}

ggSoma.neuronlist <- function(neurGG,proj=c("x","y"),...){
  neurGG <- neuronForGGMem(neurGG)
  ggSoma(neurGG,proj=proj,...)
}

ggSoma.neuron <- function(neurGG,proj=c("x","y"),...){
  neurGG <- neuronForGGMem(neurGG)
  ggSoma(neurGG,proj=proj,...)
}

ggSoma.data.frame <- function(neurGG,proj=c("x","y"),...){
    geom_circle(data=distinct(neurGG %>% filter(Label==1)),aes(x0=!!as.name(proj[1]),y0=!!as.name(proj[2]),r = W,fill=bodyid,color=bodyid),inherit.aes = FALSE,...)  # The soma
}