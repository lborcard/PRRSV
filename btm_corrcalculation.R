btm_corrcalculation <- function(data_for_corr= superdfadj,
                                Daystocompare=c(28,3),
                            celltype="%prol. t-cells",
                            genesets,
                            treatmentofinterest=c("HP"),
                            corr_value=0.8){
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
  res.cor<-res.cor %<>%
    focus("valueadj")%>%
    gather(-rowname, key = "colname", value = "cor") %>% 
    filter(abs(cor)>corr_value)%>%
    print()
  res.cor_btm<-merge(btmfamilies_char,res.cor,by.x = "geneid",by.y = "rowname")
  print(head(res.cor_btm))
    #filter(colname=="value")%>%
  plotrescor<-res.cor_btm %>%
    filter(NAME!=str_detect(string = NAME,pattern = "TBA"))%>%
    filter(Family_name!="various")%>%
    ggplot(aes(y=NAME, x=Family_name,fill=cor)) +

    geom_tile()+
    ggtitle(paste("Day",Daystocompare[1],(celltype)," vs ","Day",Daystocompare[2],"(ES)",paste(treatmentofinterest),sep = ""),
            subtitle = paste("cutoff p<0.05"," correlation value >abs(",corr_value,")",sep = ""))+
    coord_fixed(ratio = 1)+
    
    
    theme_classic()+
    theme(axis.text.x = element_text(angle = 45,hjust = 1))
    #coord_equal()
  return(plotrescor)

}

