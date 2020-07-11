list_days<-list(c(14,0),c(14,3),c(14,7),c(28,0),c(28,3),c(28,7),c(42,0),c(42,3),c(42,7),c(56,0),c(56,3),c(56,7))
listofcorr<-lapply(list_days,
                   FUN = function(x) {btm_corrcalculation(Daystocompare = x,
                                                                    corr_value = 0)})
res.cor_btm<-bind_rows(listofcorr)
plotsofcorr<- function(x=family,corrtouse=lapply()){
  

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
  ggsave(paste(res.cor_btm$Family_name,".png"))

print(plotrescor)
dev.off()
}
return(plotrescor)