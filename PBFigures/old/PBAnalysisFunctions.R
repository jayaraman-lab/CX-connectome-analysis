# Generate a von Mises activity profile
vonMisesProf <- function(){
  # Make a bump with a von Mises profile
  k <- 2.5
  m <- pi
  rng <- c(0:15)*(pi/8)
  actProf <- exp(k*cos(rng-m))/(2*pi*besselI(k,0))
  actProf <- (actProf - min(actProf))/(max(actProf)-min(actProf))
  rawProf <- data.frame(nameid = rng,act = actProf)
  
  return(actProf)
}

# Shift to align bumps
alignBumps <- function(DSnrons, Dat, shift){
  # List all the glomeruli 
  allGloms <- c("R9","R8","R7","R6","R5","R4","R3","R2","R1",
                "L1","L2","L3","L4","L5","L6","L7","L8","L9")
  
  # Sort the EPG glomeruli by their order in the PB
  gloms_PBSort <- allGloms[which(!grepl('9',allGloms))]
  
  allAlignedDat <- list()
  
  # Step through each downstream neuron
  for (ds in 1:length(DSnrons)){
    # Get the relevant data
    DatNow <- Dat[[DSnrons[ds]]]
    if (length(DatNow) == 0)
      next
    
    # Find the average activity per glomerulus
    DatNow$glom <- DatNow$name %>% lapply(function(x) sub('.*_(L[1-9]|R[1-9]).*','\\1',x)) %>% unlist()
    DatNow_byGlom <- DatNow %>% group_by(glom, peakGlom) %>% summarize(meanAct = mean(act))
    
    # Mark the R and L PB
    DatNow_byGlom$RL <- DatNow_byGlom$glom %>% lapply(function(x) substring(x,1,1)) %>% unlist()
    
    # Shift the activity
    for (g in 1:length(gloms_PBSort)){
      
      # Split right and left
      RActDat <- DatNow_byGlom %>% filter(peakGlom == gloms_PBSort[g], RL == "R")
      RActDat$alignAct <- 0
      LActDat <- DatNow_byGlom %>% filter(peakGlom == gloms_PBSort[g], RL == "L")
      LActDat$alignAct <- 0
      
      # Shift the activity
      RAct <- RActDat$meanAct
      LAct <- LActDat$meanAct
      if (g < 9){
        RAct <- c(tail(RAct, g), head(RAct, -g))
        LAct <- c(tail(LAct, -g), head(LAct, g))
      } else {
        RAct <- c(tail(RAct, (g-8)), head(RAct, -(g-8)))
        LAct <- c(tail(LAct, -(g-8)), head(LAct, (g-8)))
      }
      RAct <- c(tail(RAct, shift), head(RAct, -shift))
      LAct <- c(tail(LAct, -shift), head(LAct, shift))
      
      # Put the shifted values back in the array
      RActDat$alignAct <- RAct
      LActDat$alignAct <- LAct
      
      shiftDat <- rbind(RActDat,LActDat)
      
      if (g == 1)
        alignedDat <- shiftDat
      else
        alignedDat <- rbind(alignedDat, shiftDat)
    }
    
    alignedDat$glom <- factor(alignedDat$glom,levels=allGloms) 
    alignedDat$peakGlom <- factor(alignedDat$peakGlom,levels=allGloms)
    
    allAlignedDat[[DSnrons[ds]]] <- alignedDat
  }
  
  return (allAlignedDat)
}

# Plot the aligned bumps
alignedPlts <- function(DSnrons, alignedEPGToDS, alignedEPGToD7ToDS){
  for (ds in 1:length(DSnrons)){
    
    D7Dat <- melt(alignedEPGToD7ToDS[[DSnrons[ds]]],id.vars = c('glom','peakGlom','RL'))
    D7Dat$input <- "D7"
    D7Dat$value <- D7Dat$value/max(D7Dat$value)
    
    if (length(alignedEPGToPFNs[[PFNs[ds]]]) > 0){
      EPGDat <- melt(alignedEPGToDS[[DSnrons[ds]]],id.vars = c('glom','peakGlom','RL'))
      EPGDat$input <- "EPG"
      EPGDat$value <- EPGDat$value/max(EPGDat$value)
      
      allDat <- rbind(EPGDat,D7Dat)
    } else
      allDat <- D7Dat
    
    plt <- ggplot(allDat) + geom_tile(aes(x=glom,y=as.factor(peakGlom),fill=value)) + 
      facet_grid(rows = vars(input), cols = vars(variable)) +
      theme_paper() + coord_equal(ratio=1) + theme(axis.text.x = element_text(angle = 90, hjust=0, vjust =1)) + ggtitle(DSnrons[ds])
    
    if (ds == 1)
      plts <- plt
    else
      plts <- plts + plt
  }
  return(plts)
}

# Calculate the mean bump shape
alignedStats <- function(DSnrons, alignedEPGToDS, alignedEPGToD7ToDS){
  
  allBumpStats <- list()
  for (ds in 1:length(DSnrons)){
    
    # Split the data by right and left
    bump_D7_R <- alignedEPGToD7ToDS[[DSnrons[ds]]] %>% filter(grepl("R",peakGlom))
    bump_D7_L <- alignedEPGToD7ToDS[[DSnrons[ds]]] %>% filter(grepl("L",peakGlom))
    
    # Find the R and L mean shape
    bump_D7_R_stats <- bump_D7_R %>% group_by(glom) %>% summarize(mean = mean(alignAct), sd = mean(alignAct))
    bump_D7_R_stats$RL <- "R"
    bump_D7_R_stats$inpt <- "D7"
    bump_D7_L_stats <- bump_D7_L %>% group_by(glom) %>% summarize(mean = mean(alignAct), sd = mean(alignAct))
    bump_D7_L_stats$RL <- "L"
    bump_D7_L_stats$inpt <- "D7"
    
    if (length(alignedEPGToDS[[DSnrons[ds]]]) > 0){
      # Split the data by right and left
      bump_EPG_R <- alignedEPGToDS[[DSnrons[ds]]] %>% filter(grepl("R",peakGlom))
      bump_EPG_L <- alignedEPGToDS[[DSnrons[ds]]] %>% filter(grepl("L",peakGlom))
      
      # Find the R and L mean shape
      bump_EPG_R_stats <- bump_EPG_R %>% group_by(glom) %>% summarize(mean = mean(alignAct), sd = mean(alignAct))
      bump_EPG_R_stats$RL <- "R"
      bump_EPG_R_stats$inpt <- "EPG"
      bump_EPG_L_stats <- bump_EPG_L %>% group_by(glom) %>% summarize(mean = mean(alignAct), sd = mean(alignAct))
      bump_EPG_L_stats$RL <- "L"
      bump_EPG_L_stats$inpt <- "EPG"
      
      # Combine into one dataframe
      bump_stats <- rbind(rbind(bump_EPG_R_stats, bump_EPG_L_stats),rbind(bump_D7_R_stats, bump_D7_L_stats))
    } else
      bump_stats <- rbind(bump_D7_R_stats, bump_D7_L_stats)
    
    allBumpStats[[DSnrons[ds]]] <- bump_stats
  }
  return(allBumpStats)
}

# Compare the inverted D7 profile to the direct EPG input
normBumpProf <- function(data, RLglom, RLin){
  guideVals <- c('R profile - from EPG neurons' = 'royalblue',
                 'R profile - from D7 neurons (inverted)' = 'royalblue',
                 'L profile - from EPG neurons' = 'tan1',
                 'L profile - from D7 neurons (inverted)' = 'tan1')
  
  plt <- ggplot()
  if (nrow(filter(data, grepl(RLglom,glom) & inpt == "EPG")) > 0)
    plt <- plt + geom_line(data = filter(data, grepl(RLglom,glom) & inpt == "EPG" & RL == RLin), 
                           aes(x=glom,y=(mean-min(mean))/(max(mean)-min(mean)),group=1, color=paste(RLin,'profile - from EPG neurons', sep = ' ')))
  plt <- plt + geom_line(data = filter(data, grepl(RLglom,glom) & inpt == "D7"  & RL == RLin), 
                         aes(x=glom,y=1-(mean-min(mean))/(max(mean)-min(mean)),
                             group=1, color=paste(RLin,'profile - from D7 neurons (inverted)',sep = ' ')), linetype = 2) +
    theme_paper() + scale_color_manual(values = guideVals, name = 'profile legend') +
    scale_x_discrete(breaks = allGloms,
                     limits = allGloms[which(grepl(RLglom,allGloms))])
  
  return(plt)
}
inputComp <- function(bumpStats, title){
  
  bumpProf_PB_R_comp_fromR <- normBumpProf(bumpStats, "R", "R")
  bumpProf_PB_R_comp_fromL <- normBumpProf(bumpStats, "R", "L")
  bumpProf_PB_L_comp_fromR <- normBumpProf(bumpStats, "L", "R") + ggtitle(title)
  bumpProf_PB_L_comp_fromL <- normBumpProf(bumpStats, "L", "L")
  
  fromRplts <- (bumpProf_PB_L_comp_fromR / bumpProf_PB_R_comp_fromR)
  fromLplts <- (bumpProf_PB_L_comp_fromL / bumpProf_PB_R_comp_fromL)
  
  BumpProf_comp <- (fromRplts | fromLplts)
  
  return(BumpProf_comp)
}
normBumpProfAngs <- function(data, RLin, inptTp){
  guideVals <- c('R PB: R profile - from EPG neurons (inverted)' = 'royalblue',
                 'L PB: R profile - from EPG neurons (inverted)' = 'royalblue',
                 'R PB: L profile - from EPG neurons (inverted)' = 'tan1',
                 'L PB: L profile - from EPG neurons (inverted)' = 'tan1',
                 'R PB: R profile - from D7 neurons (inverted)' = 'royalblue',
                 'L PB: R profile - from D7 neurons (inverted)' = 'royalblue',
                 'R PB: L profile - from D7 neurons (inverted)' = 'tan1',
                 'L PB: L profile - from D7 neurons (inverted)' = 'tan1')
  
  plt <- ggplot()
  if (nrow(filter(data, grepl("R", glom) & RL == RLin & inpt == inptTp)) > 0)
    plt <- plt + geom_line(data = filter(data, grepl("R", glom) & RL == RLin & inpt == inptTp), 
                           aes(x=ang,y=1-(mean-min(mean))/(max(mean)-min(mean)),group=1, 
                               color=paste('R PB: ',RLin,' profile - from D7 neurons (inverted)', sep = '')), 
                           linetype = 2) +
    geom_line(data = filter(data, grepl("L", glom) & RL == RLin & inpt == inptTp), 
              aes(x=ang,y=1-(mean-min(mean))/(max(mean)-min(mean)),group=1, 
                  color=paste('L PB: ',RLin,' profile - from D7 neurons (inverted)', sep = '')),
              linetype = 3) +
    theme_paper() + scale_color_manual(values = guideVals) +
    scale_x_continuous(breaks = round(sort(unique(filter(data, grepl("R", glom) & RL == RLin & inpt == inptTp)$ang)), digits = 2)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))
  
  return(plt)
}
inputCompAngs <- function(bumpStats, title, inptTp){
  pltR <- normBumpProfAngs(bumpStats, "R", inptTp) + ggtitle(paste(title, '-',inptTp, ' input',sep=''))
  pltL <- normBumpProfAngs(bumpStats, "L", inptTp)
  
  plts <- pltR + pltL
  return(plts)
}

cosFit <- function(data){
  df = data.frame(angs = data$ang, weight = data$mean)
  fit <- nls(weight ~ (C1 * cos(angs-theta)+C2), data=df, algorithm="port",
             start = list(C1=1, theta=pi, C2=0),
             lower = list(C1=0, theta=0, C2=-1),
             upper = list(C1=2, theta=2*pi, C2=2))
  return(fit)
}
vMFit <- function(data){
  df = data.frame(angs = data$ang, weight = data$mean)
  fit <- nls(weight ~ (C1 * exp(k*cos(angs-theta))/(2*pi*besselI(k,0))), data=df, algorithm="port",
             control = c(warnOnly = TRUE),
             start = list(C1=1.7, k = 2.2, theta=0),
             lower = list(C1=0, k = 0.5, theta=-pi),
             upper = list(C1=5, k = 10, theta=pi))
  return(fit)
}

meanNorm <- function(mean){
  normedMean <- (mean-min(mean))/(max(mean)-min(mean))
  return(normedMean)
}


phaseDiffCalc <- function(data, inptTp){
  
  RDat_Rin <- filter(data, grepl("R", glom) & RL == "R" & inpt == inptTp)
  RDat_Rin$mean <- meanNorm(RDat_Rin$mean)
  
  
  LDat_Rin <- filter(data, grepl("L", glom) & RL == "R" & inpt == inptTp)
  LDat_Rin$mean <- meanNorm(LDat_Rin$mean)
  
  
  RDat_Lin <- filter(data, grepl("R", glom) & RL == "L" & inpt == inptTp)
  RDat_Lin$mean <- meanNorm(RDat_Lin$mean)
  
  LDat_Lin <- filter(data, grepl("L", glom) & RL == "L" & inpt == inptTp)
  LDat_Lin$mean <- meanNorm(LDat_Lin$mean)
  
  
  if (inptTp == "D7"){
    RRFit <- cosFit(RDat_Rin)
    RLFit <- cosFit(LDat_Rin)
    LRFit <- cosFit(RDat_Lin)
    LLFit <- cosFit(LDat_Lin)
  }
  if (inptTp == "EPG") {
    RRFit <- vMFit(RDat_Rin)
    RLFit <- vMFit(LDat_Rin)
    LRFit <- vMFit(RDat_Lin)
    LLFit <- vMFit(LDat_Lin)
  }
  
  RinDiff <- as.numeric(coef(RLFit)['theta']) - as.numeric(coef(RRFit)['theta']) 
  LinDiff <- as.numeric(coef(LLFit)['theta']) - as.numeric(coef(LRFit)['theta'])
  
  diffs <- data.frame(RL = c("R","L"),
                      diff = c(RinDiff,LinDiff))
  
  return(diffs)
}