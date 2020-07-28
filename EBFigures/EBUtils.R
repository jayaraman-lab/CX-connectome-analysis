# Order the columnar neurons by their location in the EB
EBGlomOrder <- function(name,id,PBEB){
  # Specify the order of glomeruli
  if (PBEB == "EB"){
    glomOrder <- list(
      EPG = c("R1","L8","R2","L7","R3","L6","R4","L5","R5","L4","R6","L3","R7","L2","R8","L1"),
      EPGt = c("R9","L9"),
      PEG = c("R1","L8","R2","L7","R3","L6","R4","L5","R5","L4","R6","L3","R7","L2","R8","L1"),
      PEN1 = c("R2","L9","R3","L8","R4","L7","R5","L6","R6","L5","R7","L4","R8","L3","R9","L2"),
      PEN2 = c("R2","L9","R3","L8","R4","L7","R5","L6","R6","L5","R7","L4","R8","L3","R9","L2"))  
  }
  if (PBEB == "PB"){
    glomOrder <- list(
      EPG = c("L8","L7","L6","L5","L4","L3","L2","L1",
              "R1","R2","R3","R4","R5","R6","R7","R8"),
      PEG = c("L8","L7","L6","L5","L4","L3","L2","L1",
              "R1","R2","R3","R4","R5","R6","R7","R8"),
      EPGt = c("R9","L9"),
      PEN1 = c("L9","L8","L7","L6","L5","L4","L3","L2",
               "R2","R3","R4","R5","R6","R7","R8","R9"),
      PEN2 = c("L9","L8","L7","L6","L5","L4","L3","L2",
               "R2","R3","R4","R5","R6","R7","R8","R9"))
  }


# Form a unique nameid from the PB type, the PB gloms where it arborizes, and the bodyid
name  <- gsub("\\s*\\([^\\)]+\\)","",as.character(name))
name  <- gsub("PEN_a","PEN1",as.character(name))
name  <- gsub("PEN_b","PEN2",as.character(name))
nameid <- paste(name, as.character(id), sep='-')

# Extract the unique names and types
allNames <- nameid %>% unique() %>% sort()
nameLabels = name %>% unique() %>% sort()
types = strsplit(allNames,'_') %>% lapply(function(x){x[[1]]}) %>% unique() %>% unlist()
types = gsub("([\\(\\)])", "\\\\\\1",types)

# Sort the names by the order of the glomeruli in the PB and exchange the bodyid for a number
newNames = c()
for (tp in 1:length(types)){
  nmsNow <- allNames[which(grepl(paste0('^',types[tp],"_"),allNames))]
  gloms <- sapply(nmsNow,function(n) strsplit(n,'_')[[1]][2]) %>% unlist() %>% as.character() 
  gloms <- sapply(gloms,function(n) strsplit(n,'-')[[1]][1]) %>% unlist() %>% as.character() 
  glomSort <- sapply(as.character(unlist(glomOrder[types[tp]])),
                     function(g) which(grepl(g,gloms))) %>% unlist() %>% as.numeric()
  glomSort <- append(glomSort, which(!(gloms %in% as.character(unlist(glomOrder[types[tp]])))))
  nmOrder <- nmsNow[glomSort]
  allNames[which(grepl(paste0('^',types[tp],"_"),allNames))] <- nmOrder
  
  nms <- nameLabels[which(grepl(paste0('^',types[tp],"_"),nameLabels))]
  for (n in 1:length(nms)){
    nmsNow <- nmOrder[which(grepl(nms[n],nmOrder))]
    nmsNow <- paste(nms[n],seq(1:length(nmsNow)),sep='-')
    nmOrder[which(grepl(nms[n],nmOrder))] <- nmsNow
  }
  newNames = append(newNames,nmOrder)
}
# Create a lookup table between the old and new names
nmSwap <- data.frame(oldNm = allNames, newNm = newNames)

# Swap in the new names
nameid <- lapply(nameid, function(x) nmSwap$newNm[match(x, nmSwap$oldNm)]) %>%
  factor(levels = nmSwap$newNm)

}