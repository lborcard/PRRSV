corrcalculation <- function(data_for_corr= superdfadj,Daystocompare=c(28,3),
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
    filter(name==celltype)%>%
    filter(day.x==Daystocompare[1])%>%#we set the day we want to use for the FACS data
    filter(day.y==Daystocompare[2])%>%#we set the day we want to use for the enrichment data
    filter(treatment%in%treatmentofinterest)%>%
    select(animalid,genesets,es,valueadj)%>%
    pivot_wider(names_from = genesets,values_from = es)

  
  row.names(myfavcorr)<-myfavcorr$animalid
  myfavcorr<- myfavcorr[c(-1,-2)]
  print(myfavcorr)
  #myfavcorr_cor<-cor(myfavcorr)
  res.cor <- correlate(myfavcorr)
  print(head(res.cor))
  plotrescor<-res.cor %<>%
    focus("valueadj")%>%
    gather(-rowname, key = "colname", value = "cor") %>% 
    filter(abs(cor)>corr_value)%>%
    print()%>%
    #filter(colname=="value")%>%
    ggplot(aes(rowname, celltype,fill=cor)) +

    
    
    geom_tile()+
    ggtitle(paste("Day",Daystocompare[1],(celltype)," vs ","Day",Daystocompare[2],"(ES)",paste(treatmentofinterest),sep = ""),
            subtitle = paste("cutoff p<0.05"," correlation value >abs(",corr_value,")",sep = ""))+
    #coord_fixed(ratio = 1)+
    coord_equal()+
    coord_flip() +
    theme_classic()
  print((res.cor$data))
  return(plotrescor)
  myfavcorr_cor_df<-as.data.frame(myfavcorr_cor)
  zero_df <- myfavcorr_cor_df 
  #corrplot(as.matrix(zero_df),method = "color")
  corr_pvalues<-rcorr(as.matrix(zero_df))
  myvalues<-flattenCorrMatrix(corr_pvalues$r,corr_pvalues$P)
  myvalues <- myvalues%>% filter(p<0.005,row=="value" & abs(cor)>corr_value)%>%
    arrange(cor)
  myvalues$row<- celltype
  
  return(myvalues)
}

