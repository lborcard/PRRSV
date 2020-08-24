es_score<-function(day){
  ## 
  library(Biobase)
  library(dplyr)
  library(limma)
  library(GSEABase)
  library(GSVA)
  library(sigPathway)
  library(readxl)
  library(tidyverse)
  library(ggpubr)
  #countfile<-read.delim(file = "p470/Day0/extended/p470.counts.txt",check.names = FALSE,row.names = 1)
  #the count file that is read by read.delim depends on the day argument above
  countfile<-read.delim(file = paste0("p470/Day",day,"/extended/p470.counts.txt"),check.names = FALSE,row.names = 1)
  
  #the countfile must be read as a matrix to be used in gsva function
  
  eset_matrix <- as.matrix(countfile)
  btmgeneset<-gmxToG("BTM gene sets.gmx")# from the sigpathway package
  
  #small loop to create a genesets list
  for (i in 1:346){
    
    genesets[i]<-btmgeneset[[i]]$title
  }
  
  results_gsva<-lapply(X = btmgeneset,FUN = function(x){gsva(expr=eset_matrix, 
                                                             gset.idx.list=x, 
                                                             annotation,
                                                             method="gsva",
                                                             kcdf="Poisson",#poisson for discrete numbers
                                                             abs.ranking=FALSE,
                                                             min.sz=1,#minimal size of each genesets
                                                             max.sz=Inf,#maximal size of each genesets
                                                             parallel.sz=0,
                                                             parallel.type="SOCK",
                                                             mx.diff=TRUE,
                                                             tau=1,#tau=switch(method,gsva=1, ssgsea=0.25, NA),
                                                             ssgsea.norm=TRUE,
                                                             verbose=TRUE)})
  
  results_gsva_df<-data.frame(matrix(unlist(results_gsva),#we transform the list in a dataframe
                                     nrow = length(results_gsva),
                                     byrow = T),
                              stringsAsFactors = FALSE)
  
  results_gsva_df<-t(results_gsva_df)#we transpose because we want the genesets as col and samples as rows
  
  sampleid<-colnames(countfile)#we take the names of the samples
  
  colnames(results_gsva_df)<- genesets
  
  #add the rownames of the samples (now as rows)
  results_gsva_df$sampleid<-sampleid
  
  #assign a new name to the data frame with Day + day of the count file + _ES (enrichment score)
  assign(paste("Day",day,"_ES",sep = ""),results_gsva_df)
}