btm_corrcalculation <- function(data_for_corr= superdfadj,
                                Daystocompare=c(28,3),
                            celltype="%th.",
                            genesets,
                            myvalue="value",
                            treatmentofinterest="HP",
                            corr_value=0.8){
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
    select(animalid,genesets,es,myvalue)%>%
    pivot_wider(names_from = genesets,values_from = es)
  
  
  
  row.names(myfavcorr)<-myfavcorr$animalid
  myfavcorr<- myfavcorr[c(-1,-2)]
  #myfavcorr_cor<-cor(myfavcorr)
  res.cor <- correlate(myfavcorr,method = "kendall")
  res.cor<-res.cor %<>%
    focus(myvalue)%>%
    gather(-rowname, key = "colname", value = "cor") %>% 
    mutate(rsquared=sign(cor)*cor^2)%>%
    #filter(rsquared>corr_value)%>%
    print()
  res.cor_btm<-merge(btmfamilies_char,res.cor,by.x = "geneid",by.y = "rowname")
    #filter(colname=="value")%>%

  days<-paste0(Daystocompare[1],"vs",Daystocompare[2])
  namesettings<-paste(treatmentofinterest,celltype,sep = "_")
  res.cor_btm$daystocompare<-days
  res.cor_btm$namegraph<-namesettings
  print(head(res.cor_btm))
  return(res.cor_btm)
  #########
  #########
  #########
  plotrescor<-res.cor_btm %>%
    filter(NAME!=str_detect(string = NAME,
                            pattern = "TBA"))%>%
    filter(Family_name!="various")%>%
    filter(Family_name!="cell cycle")%>%
    ggplot(aes(y=NAME, x=Family_name,fill=cor)) +
    geom_tile()+
    ggtitle(paste("Day",
                  Daystocompare[1],
                  (celltype)
                  ," vs ","Day",
                  Daystocompare[2],
                  "(ES)",
                  paste(treatmentofinterest),
                  sep = ""),
            subtitle = paste("cutoff p<0.05 ",
                             expression(R^2),
                             ">",corr_value,sep = ""))+
    coord_fixed(ratio = 1)+
    
    
    theme_classic()+
    scale_fill_gradient2(
      #low = "grey2",
      #high = "blue",
      na.value = "grey50",
      guide = "colourbar",
      aesthetics = "fill",
      midpoint=0,
      limits=c(-1,1)
    )+
    theme(axis.text.x = element_text(angle = 45,hjust = 1))
    #coord_equal()
  return(plotrescor)

}

