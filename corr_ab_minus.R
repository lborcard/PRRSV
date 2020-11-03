ab_correlation <- function(data_for_corr= superdf_minus,
                                day_prol=c(14,28,56),
                                day_es="7min0",
                                ab_day=c(10,14,21,28,42,63),
                                celltype="%prol. t-cells",
                                genesets,
                                myvalue=c("ab_inhib","valueadj"),
                                treatmentofinterest="HP",
                                corr_value=0){
  
  library(tidyverse)
  library(corrplot)
  library(corrr)
  library(data.table)
  
  #we create the dataframe from the data we want to check
  res.cor1<-data.frame()
  for (x in ab_day) {
   myfavcorr1<-data_for_corr %>%
     filter(name==celltype)%>%
     filter(day_min==day_es)%>%
     filter(day_facs==28)%>%
     filter(day_ab==x)%>%#we set the day we want to use for the enrichment data
     filter(treatment%in%treatmentofinterest)%>%
     mutate(id_day=paste(animalid,day_ab,sep = "_"))%>%
     select(id_day,genesets,es_min,myvalue[1])%>%
     pivot_wider(names_from = genesets,values_from = es_min)
   
   
   print(myfavcorr1)
   row.names(myfavcorr1)<-myfavcorr1$animalid
   myfavcorr1<- myfavcorr1[-c(1)]
   print(myfavcorr1[3])
   
   res.cor <- correlate(myfavcorr1,method = "pearson")
   res.cor<-res.cor %<>%
     focus(myvalue[1])%>%
     gather(-rowname, key = "colname", value = "cor") %>% 
     mutate(ab_inhib_r2=sign(cor)*cor^2,geneid=rowname,day_ab=paste(x,"ab",sep = "_"))%>%
     select(geneid,ab_inhib_r2,day_ab)%>%
     print()
   
   res.cor1<-bind_rows(res.cor1,res.cor)
 }
  print(res.cor1)
  
  ######
  ### second correlation
  ######
  res.cor2<-data.frame()
  for (x in day_prol) {
    myfavcorr1<-data_for_corr %>%
      filter(name==celltype)%>%
      filter(day_min==day_es)%>%
      filter(day_facs==x)%>%
      filter(day_ab==3)%>%#we set the day we want to use for the enrichment data
      filter(treatment%in%treatmentofinterest)%>%
      mutate(id_day=paste(animalid,day_facs,sep = "_"))%>%
      select(id_day,genesets,es_min,myvalue[2])%>%
      pivot_wider(names_from = genesets,values_from = es_min)
    
    
    print(myfavcorr1)
    row.names(myfavcorr1)<-myfavcorr1$animalid
    myfavcorr1<- myfavcorr1[-c(1)]
    print(myfavcorr1[3])
    
    res.cor <- correlate(myfavcorr1,method = "pearson")
    res.cor<-res.cor %<>%
      focus(myvalue[2])%>%
      gather(-rowname, key = "colname", value = "cor") %>% 
      mutate(facs_r2=sign(cor)*cor^2,geneid=rowname,day_facs=paste(x,"facs",sep = "_"))%>%
      select(geneid,facs_r2,day_facs)%>%
      print()
    
    res.cor2<-bind_rows(res.cor2,res.cor)
  }
  print(res.cor2)
  #ab_inhib_r2 to 
  #day_ab day_facs in res.cor2
  #res.cor2 to res.cor
  
  ##Merge the two correlations
  
  res.cor<-merge(res.cor1,res.cor2)
  print(res.cor)
  res.cor_btm<-merge(btmfamilies_char,res.cor,by="geneid")

  days<-paste0(day_prol,"vs",day_es)
  namesettings<-paste(treatmentofinterest,celltype,sep = "_")
  res.cor_btm$daystocompare<-days
  res.cor_btm$namegraph<-namesettings
  print(head(res.cor_btm))
  print(ab_day)
  #return(res.cor_btm)
  #########
  #########
  #########
  print(str(res.cor_btm))
  plotrescor<-res.cor_btm %>%
    filter(NAME!=str_detect(string = NAME,
                            pattern = "TBA"))%>%
    filter(Family_name%in%c("ag presentation","B cells","IFN-I","inflammation","myeloid cells","NK/T cells"))%>%
    ggplot() +
    geom_tile(aes(x=as.factor(day_ab),y=NAME,fill=ab_inhib_r2))+
    geom_tile(aes(x=as.factor(day_facs),y=NAME,fill=facs_r2))+
    ggtitle(paste(str_replace(day_es,"min","/"),treatmentofinterest," "),
            subtitle = paste("cutoff p<0.05 ",
                             expression(R^2),
                             ">",corr_value,sep = ""))+
    scale_x_discrete()+
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
    facet_wrap(.~Family_name,ncol = 6)+ 
    
    
    theme_classic()+
    
    theme(axis.text.x = element_text(angle = 45,hjust = 1))
  #coord_equal()
  return(plotrescor)
  
}

