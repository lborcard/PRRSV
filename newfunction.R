newfunction <- function(data_for_corr= superdfadj,
                            Daystocompare=14,#2nd for the enrichment score
                            celltype="%prol. t-cells",
                            genesets,
                            treatmentofinterest=c("HP"),
                            corr_value=0){
  #
  library(tidyverse)
  library(corrplot)
  library(Hmisc)
  library(corrr)
  library(data.table)
  #we create the dataframe from the data we want to check
  myfavcorr<-data_for_corr %>%
    group_by()%>%
    filter(name==celltype)%>%
    filter(day.x==Daystocompare[1])%>%#we set the day we want to use for the FACS data
    filter(treatment%in%treatmentofinterest)%>%
    select(animalid,genesets,es,valueadj,day.y)%>%
    pivot_wider(names_from = c(genesets,day.y),values_from = es)%>%
    select(-animalid)
  res.cor <- correlate(myfavcorr)
  res.cor$id_col<-rownames(res.cor)
  str(res.cor)
  res.cor_long<-res.cor%>%pivot_longer(cols = c(-valueadj,-id_col),names_to = c("es","day"),names_sep = "_")
  return(res.cor)
}

