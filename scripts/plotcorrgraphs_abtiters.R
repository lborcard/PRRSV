##Antibody titers graphs
for(x in unique(superdf_abtiter$treatment)){
  for(y in cells_type_map){
    listofcorr<-lapply(list_days_abtiters,
                       FUN = function(p) {btm_corrcalculation_abtiters(treatmentofinterest = x,
                                                              celltype = y,
                                                              Daystocompare = p,
                                                              corr_value = 0)})
    mycorr<-bind_rows(listofcorr)
    actualpath<-paste0("figures/",x,"/",y,"/")
    print(actualpath)
    familynames<-unique(btmfamilies_char$Family_name)
    print(familynames)
    for(f in familynames){
      print(paste0(actualpath,"/",str_replace(f,pattern = "[\\/\\s\\-]",replacement = "_"),".png"))
      outfile<-paste0(actualpath,"/",str_replace(f,pattern = "[\\/\\s\\-]",replacement = "_"),".png")
      figname<-paste0(str_replace(f,pattern = "[\\/\\s\\-]",replacement = "_"),".png")
      plotrescor<-mycorr %>%
        filter(NAME!=str_detect(string = NAME,
                                pattern = "TBA"))%>%
        filter(Family_name==f)%>%
        
        ggplot(aes(y=NAME,
                   x=daystocompare,
                   fill=(cor/abs(cor))*cor^2)) +
        geom_tile()+
        coord_fixed(ratio = 1)+
        
        #facet_wrap(.~Family_name,nrow = 1)+
        ggtitle(label = paste(mycorr$namegraph
                              ,"_",
                              f))+
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
        theme(axis.text.x = element_text(angle = 45,
                                         hjust = 1))+
        ggsave(path = actualpath,filename = figname,device = "png")
      
    }
  }
}
