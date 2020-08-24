list_days<-list(c(14,0),c(14,3),c(14,7),c(28,0),c(28,3),c(28,7),c(42,0),c(42,3),c(42,7),c(56,0),c(56,3),c(56,7))
list_days<-list(c(14,0),c(14,3),c(14,7),c(28,0),c(28,3),c(28,7),c(42,0),c(42,3),c(42,7),c(56,0),c(56,3),c(56,7))
list_days_abtiters <- list(c(0,0),c(7,0),c(10,0),c(21,0),c(49,0),c(0,3),c(7,3),c(10,3),c(21,3),c(49,3),c(0,7),c(7,7),c(10,7),c(21,7),c(49,7))
cells_type_map<-cells_type[c(2,6,7,9,10,15,16)]
listofcorr<-lapply(list_days,
                   FUN = function(x) {btm_corrcalculation(Daystocompare = x,
                                                                    corr_value = 0)})
res.cor_btm<-bind_rows(listofcorr)
plotsofcorr<- function(x=family,corrtouse=lapply()){
map(x=unique(superdf$stimulus),map(y=cells_type_map,map()))

plotrescor<-res.cor_btm %>%
  filter(NAME!=str_detect(string = NAME,
                          pattern = "TBA"))%>%
  filter(Family_name=="ag presentation")%>%
  ggplot(aes(y=NAME, x=daystocompare,fill=cor)) +
  geom_tile()+
  coord_fixed(ratio = 1)+
  
  #facet_wrap(.~Family_name,nrow = 1)+
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
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  ggsave(path = "figures/",filename = paste(namegraph,"_",res.cor_btm$Family_name,".png"))

print(plotrescor)
dev.off()
}
return(plotrescor)