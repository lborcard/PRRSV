ab_correlation <- function(data_for_corr= superdf_minus,
                                Daystocompare=c(28,"3min0"),
                                ab_day=7,
                                celltype="%th",
                                genesets,
                                myvalue=c("ab_inhib","valueadj"),
                                treatmentofinterest="HP",
                                corr_value=0){
  #
  library(tidyverse)
  library(corrplot)
  #library(Hmisc)
  library(corrr)
  library(data.table)
  
  #we create the dataframe from the data we want to check
  
  myfavcorr1<-data_for_corr %>%
    filter(name==celltype)%>%
    filter(day_facs==Daystocompare[1])%>%#we set the day we want to use for the FACS data
    filter(day_min==Daystocompare[2])%>%
    filter(day_ab==ab_day)%>%#we set the day we want to use for the enrichment data
    filter(treatment%in%treatmentofinterest)%>%
    select(animalid,genesets,es_min,myvalue[1])%>%
    pivot_wider(names_from = genesets,values_from = es_min)
  
  
  print(myfavcorr1)
  row.names(myfavcorr1)<-myfavcorr1$animalid
  myfavcorr1<- myfavcorr1[c(-1)]
  print(myfavcorr1)
  #myfavcorr_cor<-cor(myfavcorr)
  res.cor1 <- correlate(myfavcorr1,method = "pearson")
  res.cor1<-res.cor1 %<>%
    focus(myvalue[1])%>%
    gather(-rowname, key = "colname", value = "cor") %>% 
    mutate(rsquared=sign(cor)*cor^2)%>%
    #filter(rsquared>corr_value)%>%
    print()
  ######
  ### second correlation
  ######
  myfavcorr2<-data_for_corr %>%
    filter(name==celltype)%>%
    filter(day_facs==Daystocompare[1])%>%#we set the day we want to use for the FACS data
    filter(day_min==Daystocompare[2])%>%
    filter(day_ab==ab_day)%>%#we set the day we want to use for the enrichment data
    filter(treatment%in%treatmentofinterest)%>%
    select(animalid,genesets,es_min,myvalue[2])%>%
    pivot_wider(names_from = genesets,values_from = es_min)
  
  
  print(myfavcorr2)
  row.names(myfavcorr2)<-myfavcorr2$animalid
  myfavcorr2<- myfavcorr2[c(-1)]
  print(myfavcorr2)
  #myfavcorr_cor<-cor(myfavcorr)
  res.cor2 <- correlate(myfavcorr2,method = "pearson")
  res.cor2<-res.cor2 %<>%
    focus(myvalue[2])%>%
    gather(-rowname, key = "colname", value = "cor") %>% 
    mutate(rsquared=sign(cor)*cor^2)%>%
    #filter(rsquared>corr_value)%>%
    print()
  ##Bind the two correlations
  res.cor<-bind_cols(res.cor1,res.cor2)
  print(res.cor)
  res.cor_btm<-merge(btmfamilies_char,res.cor,by.x = "geneid",by.y = "rowname...1")
  #filter(colname=="value")%>%
  
  days<-paste0(Daystocompare[1],"vs",Daystocompare[2])
  namesettings<-paste(treatmentofinterest,celltype,sep = "_")
  res.cor_btm$daystocompare<-days
  res.cor_btm$namegraph<-namesettings
  print(head(res.cor_btm))
  #return(res.cor_btm)
  #########
  #########
  #########
  print(str(res.cor_btm))
  plotrescor<-res.cor_btm %>%
    filter(NAME!=str_detect(string = NAME,
                            pattern = "TBA"))%>%
    filter(Family_name=="ag presentation")%>%
    ggplot() +
    geom_tile(aes(x=ab_day,y=NAME),fill="rsquared...4")+
    geom_tile(aes(x=Daystocompare[2],y=NAME),fill="rsquared...8")+
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
    scale_fill_gradient2(
      #low = "grey2",
      #high = "blue",
      na.value = "grey50",
      guide = "colourbar",
      aesthetics = "fill",
      midpoint=0,
      limits=c(-1,1)
    )+
    coord_fixed(ratio = 1)+
    
    
    theme_classic()+
    
    theme(axis.text.x = element_text(angle = 45,hjust = 1))
  #coord_equal()
  return(plotrescor)
  
}

