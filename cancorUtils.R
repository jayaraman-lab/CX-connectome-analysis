library(PMA)
cancorInputOutput <- function(inputs,outputs){
  inputMat <- pivot_wider(inputs,id_cols=c("type.from","roi"),names_from = type.to,values_from = knownWeightRelative,values_fill = list(knownWeightRelative=0))
  outputMat <- pivot_wider(outputs,id_cols=c("type.to","roi"),names_from = type.from,values_from = knownOutputContribution,values_fill = list(knownOutputContribution=0))
  #inputMat <- filter(inputMat,(rowSums(inputMat[,3:ncol(inputMat)])>0.2))
  #outputMat <- filter(outputMat,(rowSums(outputMat[,3:ncol(outputMat)])>0.2))
  inputMatRaw <- t(data.matrix(inputMat[,3:ncol(inputMat)]))
  colnames(inputMatRaw) <- inputMat$type.from
  outputMatRaw <- t(data.matrix(outputMat[,3:ncol(outputMat)]))
  colnames(outputMatRaw) <- outputMat$type.to
  #estim.regul(inputMatRaw, outputMatRaw, grid1 = seq(0.0001,0.0002,length.out = 2), grid2 = seq(0.0001,0.0002,length.out = 2), plt = TRUE)
  perm.out <- CCA.permute(inputMatRaw,outputMatRaw,typex = "standard",typez="standard",upos=TRUE,vpos = TRUE,standardize = TRUE)
  CCA(inputMatRaw,outputMatRaw,typex = "standard",typez="standard",upos=TRUE,vpos = TRUE,K=100,standardize = TRUE,penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz,
      v=perm.out$v.init,xnames=inputMat$type.from,znames=outputMat$type.to)
}