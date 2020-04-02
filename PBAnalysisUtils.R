#Generate a von Mises activity profile
vonMisesProf <- function(){
  # Make a bump with a von Mises profile
  k = 2.5
  m = pi
  rng = c(1:16)*(2*pi/15) - 2*pi/15
  actProf = exp(k*cos(rng-m))/(2*pi*besselI(k,0))
  actProf = (actProf - min(actProf))/(max(actProf)-min(actProf))
  rawProf = data.frame(nameid = rng,act = actProf)
  
  return(actProf)
  
  #  ggplot(rawProf,aes(x = nameid,y=act)) + geom_line() + 
  #    theme_classic() + coord_cartesian(expand=FALSE) + scale_x_continuous(breaks = seq(0,2*pi, by=pi/2), labels = c("0","pi/2","pi","3pi/2","2pi"))
  #  ggsave("EPGActProfile.pdf", plot = last_plot(), device='pdf', 
  #         path = "C:\\Users\\turnerevansd\\Dropbox (HHMI)\\FIBSEM CX Paper Jan 2020\\Figures\\PB\\",
  #        scale = 1.5, width = 20, height = 15, units ="cm", dpi = 600, limitsize = TRUE)
}

# Convert a connection table into a weight matrix
WMFromConTab <- function(conTab, preids, postids){
  
  justWeights = matrix(0L,nrow=length(preids),ncol=length(postids),)
  
  for (i in 1:length(preids)){
    for (j in 1:length(postids)){
      w = conTab[which(conTab$nameid == preids[i] & conTab$partnerid == postids[j]),]$ROIweight
      if (length(w) > 0)
        justWeights[i,j] = w
    }
  }
  return(justWeights)
}

# Order glomeruli according to their location in the PB
PBGlomSort <- function(gloms){
  sortedGloms = append(
    sort(
      as.character(gloms[which(grepl("L",gloms))]),
      decreasing = TRUE),
    sort(
      as.character(gloms[which(grepl("R",gloms))]),
      decreasing = FALSE))
  
  return(sortedGloms)
}

# Get the EPG to PB-X type transform matrix
EPGToD7Mat <- function(PBNronTypeStr){
  
  # Get the EPG to D7 connections
  EPG_2_Delta7 <- getConnectionTable_forSubset(
    unique(neuprint_search("EPG\\\\\\\\(.*")$bodyid),
    unique(neuprint_search("Delta7.*")$bodyid),
    "PB")
  EPG_2_Delta7$nameid = paste(as.character(EPG_2_Delta7$name.from), as.character(EPG_2_Delta7$from), sep = "_")
  EPG_2_Delta7$partnerid = paste(as.character(EPG_2_Delta7$name.to), as.character(EPG_2_Delta7$to), sep = "_")
  
  # Get the D7 to EPG connections weight matrix
  Delta7_2_PBTp <- getConnectionTable_forSubset(
    unique(neuprint_search("Delta7.*")$bodyid),
    unique(neuprint_search(PBNronTypeStr)$bodyid),
    "PB")
  Delta7_2_PBTp$nameid = paste(as.character(Delta7_2_PBTp$name.from), as.character(Delta7_2_PBTp$from), sep = "_")
  Delta7_2_PBTp$partnerid = paste(as.character(Delta7_2_PBTp$name.to), as.character(Delta7_2_PBTp$to), sep = "_")
  
  # Get the EPG and D7 names and glomeruli
  allEPGids <-EPG_2_Delta7$nameid %>% as.character() %>% unique() %>% PBGlomSort()
  allD7ids <- EPG_2_Delta7$partnerid %>% as.character() %>% unique() %>% sort()
  allPBTpids <- Delta7_2_PBTp$partnerid %>% as.character() %>% unique() %>% sort()
  
  # Create weight matrices for the connections
  justWeights_EPG_2_D7 <-  WMFromConTab(EPG_2_Delta7,allEPGids,allD7ids)
  justWeights_D7_2_PBTp <-  WMFromConTab(Delta7_2_PBTp,allD7ids,allPBTpids)
  
  # Multiply them together to get one weight matrix
  transMat <- (justWeights_EPG_2_D7 %*% justWeights_D7_2_PBTp)
  rownames(transMat) <- allEPGids
  colnames(transMat) <- allPBTpids
  
  return(transMat)
}

# Given an activity profile, calculate the output profile
EPGToD7Transform <- function(PBTpStr,actProf){
  
  # Order the activity profile according to the PB glomeruli
  actProf$glom <- factor(actProf$glom, levels = PBGlomSort(actProf$glom))
  actProf <- actProf[order(actProf$glom),]
  
  # Get the EPG to D7 to EPG transform matrix
  transMat <- EPGToD7Mat(PBTpStr)
  
  # Get the EPG and output type names
  allEPGids = rownames(transMat)
  allPBTpids = colnames(transMat)
  
  # Get the glomeruli ids for the EPGs
  allEPGGloms = strsplit(allEPGids,"_") %>% lapply(function(x){tail(x,2)[1]}) %>% unlist()
  allPBTpGloms = strsplit(allPBTpids,"_") %>% 
    lapply(function(x){head(x,3)[1+length(strsplit(PBTpStr[1],"_")[[1]])]}) %>% unlist()
  
  # Initialize data frames to hold all of the circular permutations
  EPGOutputs_all <- data.frame(glom = c(), act = c())
  EPGInputs_all <- data.frame(glom = c(), act = c())
  
  # Permute through the glomeruli
  for (pkGlom in 1:16){
    actProf$act <- c(tail(actProf$act, -1), head(actProf$act, 1))
    
    # Get the EPG activity given the assumed profile
    EPGact = numeric(length(allEPGids))
    for (EPG in 1:length(allEPGids)){
      EPGact[EPG] <- actProf[which(actProf$glom == allEPGGloms[EPG]),]$act
    }
    EPGProf <- data.frame(nameid = allEPGids, act = EPGact, glom = allEPGGloms)
    
    # Calculate the Delta 7 output activity onto the EPGs for the given EPG activity profile
    D7outAct = as.vector(EPGProf$act %*% transMat)
    D7outAct = (D7outAct - min(D7outAct))/(max(D7outAct)-min(D7outAct))
    PBTpInputs = data.frame(nameid = allPBTpids,
                            act= D7outAct,
                            glom = allPBTpGloms)
    
    # Average the activity over the given glomerulus
    EPGProf_byGlom <- EPGProf %>% group_by(glom) %>% summarize(act = mean(act))
    EPGProf_byGlom$glom <- factor(EPGProf_byGlom$glom,
                                  levels=PBGlomSort(unique(allEPGGloms)))
    EPGProf_byGlom <- EPGProf_byGlom[order(EPGProf_byGlom$glom),]
    EPGProf_byGlom$act <- c(tail(EPGProf_byGlom$act, pkGlom), 
                            head(EPGProf_byGlom$act, -pkGlom))
    colnames(EPGProf_byGlom) <- c("glom",paste0("act",pkGlom))
    
    PBTpInput_byGlom <- PBTpInputs %>% group_by(glom) %>% summarize(act = mean(act))
    PBTpInput_byGlom$glom <- factor(PBTpInput_byGlom$glom,
                                    levels=PBGlomSort(unique(allPBTpGloms)))
    PBTpInput_byGlom <- PBTpInput_byGlom[order(PBTpInput_byGlom$glom),]
    PBTpInput_byGlom$act <- c(tail(PBTpInput_byGlom$act, pkGlom), 
                              head(PBTpInput_byGlom$act, -pkGlom))
    colnames(PBTpInput_byGlom) <- c("glom",paste0("act",pkGlom))
    
    # Store the values across circular permutations
    if (pkGlom == 1){
      EPGOutputs_all <- EPGProf_byGlom
      PBTpInputs_all <- PBTpInput_byGlom
    } else {
      EPGOutputs_all <- merge(EPGOutputs_all, EPGProf_byGlom, by="glom")
      PBTpInputs_all <- merge(PBTpInputs_all, PBTpInput_byGlom, by="glom")
    }
  }  
  return(PBTpInputs_all)
}

#Calculate the mean and SD of D7 outputs
D7ShapeMeanSD <- function(inputAll){
  
  # Sort the PB gloms according to their location
  sortedGloms <- PBGlomSort(inputAll$glom)
  
  # Group by glom and calculate statistics
  D7OutputShape <- melt(inputAll,id="glom")
  D7OutputShape$glom <- factor(D7OutputShape$glom,levels = sortedGloms)
  D7OutputShape_stats <- D7OutputShape %>% 
    group_by(glom) %>% 
    summarize(mean = mean(value), sd = sd(value))
  
  return(D7OutputShape_stats)
}

# Fit a cosine to the D7 outpur profile
D7CosFit <- function(D7OutputShape_stats){
  df = data.frame(pts = c(1:nrow(D7OutputShape_stats)), weight = D7OutputShape_stats$mean)
  fit <- nls(weight ~ (C1 * cos(C2*(pts-theta))+C3), data=df, algorithm="port",
             start = list(C1=0.4, C2=0.8, theta=0, C3=0.4),
             lower = list(C1=0, C2=0, theta=0, C3=0),
             upper = list(C1=2, C2=2, theta=20, C3=1))
  return(fit)
}

# Plot the PB activity profile for a given set of activity values
PBActPlot <- function(D7OutputShape_stats,fit){
  
  dffit <- data.frame(glom = D7OutputShape_stats$glom)
  dffit$wMean <-  predict(fit, newdata=dffit)
  
  # Plot the data
  PBEPGToD7Prof <- ggplot() + 
    geom_line(data = D7OutputShape_stats, aes(x=glom,y=mean-sd,group=1),color="gray") +
    geom_line(data = D7OutputShape_stats, aes(x=glom,y=mean+sd,group=1),color="gray") +
    geom_line(data = D7OutputShape_stats, aes(x=glom,y=mean,group=1),color="black") +
    geom_line(data = dffit, aes(x=glom,y=wMean,group=1),color="red",linetype='dotted') +
    theme_classic() + coord_cartesian(ylim=c(0,1),expand=FALSE) +
    xlab('PB glom.') + ylab('activity')
  
  return(PBEPGToD7Prof)
}