ab_facs_correlation <- function(data_for_corr= superdf_minus,
                                day_prol=c(14,28,42,56),
                                day_es=c("3min0","7min0","7min3"),
                                ab_day=c(10,14,21,28,42,63),
                                celltype="%prol. t-cells",
                                btm_family="ag presentation",
                                myvalue=c("ab_inhib","valueadj"),
                                treatmentofinterest="LP",
                                corr_value=0){
  
  library(tidyverse)
  library(corrplot)
  library(corrr)
  library(data.table)
  
 res<-data.frame()
  for (p in day_es) {
  #for loop containing the different days of the RNA seq data  
 
  res.cor1<-data.frame()
  for (x in ab_day) {
    #for loop for the antibody data
   myfavcorr1<-data_for_corr %>%
     filter(name==celltype)%>%
     filter(day_min==p)%>%
     
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
   
   res.cor <- correlate(myfavcorr1,method = "pearson") # calculation of the first correlation of AB vs Enrichment scores
   res.cor<-res.cor %<>%
     focus(myvalue[1])%>%
     gather(-rowname, key = "colname", value = "cor") %>% 
     mutate(ab_inhib_r2=sign(cor)*cor^2,geneid=rowname,day_ab=paste(x,"ab",sep = "_"))%>%
     select(geneid,ab_inhib_r2,day_ab)%>%
     print()
   
   res.cor1<-bind_rows(res.cor1,res.cor) # we bind the rows of the correlations (1 new one per loop through the Days)
 }
  print(res.cor1)
  
  ######
  ### second correlation for the facs data
  ######
  res.cor2<-data.frame()
  for (x in day_prol) {
    myfavcorr1<-data_for_corr %>%
      filter(name==celltype)%>%
      filter(day_min==p)%>%
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
 
  ##Merge the two correlations
  res.cor<-merge(res.cor1,res.cor2)
  # each loop we report the day of the enrichment score for further analysis
  res.cor<-mutate(res.cor,day_min=p) 
  #we add the newly created correlation data frame to the growing dataframe
  res<-bind_rows(res,res.cor)
  } 
 
  print(res)
  
  res.cor_btm<-merge(btmfamilies_char,res,by="geneid") 

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
  res.cor_btm$day_min<-factor(x=res.cor_btm$day_min,
                              levels = c("3min0","7min0","7min3"),
                              labels = c("3/0","7/0","7/3"))
  plotrescor<-res.cor_btm %>%
    filter(NAME!=str_detect(string = NAME,
                            pattern = "TBA"))%>%
    filter(Family_name==btm_family)%>%
    
    ggplot() +
    geom_tile(aes(x=as.factor(day_ab),y=NAME,fill=ab_inhib_r2))+
    geom_tile(aes(x=as.factor(day_facs),y=NAME,fill=facs_r2))+
    ggtitle(paste(celltype,treatmentofinterest," "),
            subtitle = paste("cutoff p<0.05 ",
                             expression(R^2),
                             ">",corr_value,sep = ""))+
    scale_x_discrete(limits=c("10_ab","14_ab","21_ab",
                              "28_ab","42_ab","63_ab","14_facs","28_facs","42_facs","56_facs"))+
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
    facet_wrap(.~day_min,ncol = 6)+ 
    
    ylab("Genesets")+
    xlab(element_blank())+
    labs(fill="R squared")+
    theme_classic()+
    
    theme(axis.text.x = element_text(angle = 45,hjust = 1),
          strip.text = element_text(size = 18))

  #coord_equal()
  ggsave(plot = plotrescor,
         filename = paste(str_replace(btm_family,"/","_"),
                                            treatmentofinterest,".png"),
         device = "png",dpi = 150,width = 10,height = 8)
  return(plotrescor)
  
}

