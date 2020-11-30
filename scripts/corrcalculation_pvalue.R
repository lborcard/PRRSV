#correlation with pvalues
corrcalculation_pvalue <- function(data_for_corr= superdf_clean,
                            Daystocompare=c(56,#1st number day of the facs data 
                                            3),#2nd for the enrichment score
                            celltype="%prol. t-cells",
                            genesets,
                            treatmentofinterest=c("HP"),
                            corr_value=0){
  #
  library(tidyverse)
  library(corrplot)
  #library(Hmisc)
  library(corrr)
  library(data.table)
  #we create the dataframe from the data we want to check
  myfavcorr<-data_for_corr %>%
    filter(name==celltype)%>%
    filter(day.x==Daystocompare[1])%>%#we set the day we want to use for the FACS data
    filter(day.y==Daystocompare[2])%>%#we set the day we want to use for the enrichment data
    filter(treatment%in%treatmentofinterest)%>%
    select(animalid,genesets,es,value)%>%
    pivot_wider(names_from = genesets,values_from = es)
  
  myfavcorr<- myfavcorr[-c(1)]
  corr_pvalues<-rcorr(as.matrix(myfavcorr))
  myvalues<-flattenCorrMatrix(corr_pvalues$r,corr_pvalues$P)%>%
    filter(p<0.05,row=="value" & abs(cor)>corr_value)%>%
    mutate(rsquared=cor^2*sign(cor),geneid=column)%>%
    select(geneid,rsquared,p)
  print(myvalues)

  #myvalues <- myvalues%>% filter(p<0.05,row=="value" & abs(cor)>corr_value)%>%
   # mutate(rsquared=cor^2*sign(cor))

  
  return(myvalues)
}

