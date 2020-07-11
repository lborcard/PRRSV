for(x in treatment){
  for(y in cells_type_map){
    listofcorr<-lapply(list_days,
                       FUN = function(p) {btm_corrcalculation(treatmentofinterest = x,
                                                              celltype = y,
                                                              Daystocompare = p,
                                                              corr_value = 0)})
    mycorr<-bind_rows(listofcorr)
    familyname<-"ag presentation"
    plotrescor<-mycorr %>%
      filter(NAME!=str_detect(string = NAME,
                              pattern = "TBA"))%>%
      filter(Family_name==familyname)%>%
      
      ggplot(aes(y=NAME,
                 x=daystocompare,
                 fill=cor)) +
      geom_tile()+
      coord_fixed(ratio = 1)+
      
      #facet_wrap(.~Family_name,nrow = 1)+
      ggtitle(label = paste(mycorr$namegraph
                            ,"_",
                            familyname))+
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
      ggsave(path = "figures/",
             filename = paste(str_replace(mycorr$namegraph,
                                          pattern = "%",replacement = "percent_"),familyname,".png"),device = "png")
    
  }
}
