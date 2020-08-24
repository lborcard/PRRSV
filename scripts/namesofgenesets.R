#extract the names of the genesets
genesets_ids<-data.frame()
btmfamilies$geneid<-(str_extract(btmfamilies$NAME,pattern = "(?<=\\()(.\\d.{0,4}?)(?=\\))"))
genesets_ids$name<-str_remove(genesetsnames,pattern =  "(?<=\\()(.\\d.{0,4}?)(?=\\))")
library(tidyverse)
