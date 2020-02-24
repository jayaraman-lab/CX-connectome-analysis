supertype <- function(types){
  supertype <- types
  supertype[grepl("FB.*",types)] <- "FB tangential"
  supertype[grepl("PFL.*",types)] <- "PFL"
  supertype[grepl("PFN.*",types)] <- "PFN"
  supertype[grepl("PFR.*",types)] <- "PFR"
  supertype[grepl("PEN.*",types)] <- "PEN"
  supertype[grepl("LN.*",types)] <- "LN"
  supertype[grepl("Delta0.*",types)] <- "Delta0"
  supertype[grepl("Delta6.*",types)] <- "Delta6"
  supertype[grepl("ExR.*",types)] <- "ExR"
  supertype[grepl("FC.*",types)] <- "FC"
  supertype[grepl("^FR.*",types)] <- "FR"
  supertype[grepl("FS.*",types)] <- "FS"
  supertype[grepl("LCN.*",types)] <- "LN"
  supertype[grepl("^R[1-6].*",types)] <- "Ring neuron"
  supertype[grepl("SA.*",types)] <- "SA"
  supertype[grepl("SpsP.*")] <- "SPS-PB"
  supertype
}

megatype <- function(supertypes){
  megatype <- supertypes
  megatype[supertypes %in% c("PFL","PFN","PFR","PFGs")] <- "FB Columnar"
  megatype[supertypes %in% c("EPG","PEG","PEN")] <- "EB Columnar"
  megatype[supertypes %in% c("Delta6","Delta0")] <- "FB interneurons"
  megatype[supertypes %in% c("FC","FR","FS")] <- "FB outputs"
  megatype
}