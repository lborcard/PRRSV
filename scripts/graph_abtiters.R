dpi_of_interest=49
day_of_es=3
genesetsofinterest=c("M139")
graph_abtiters<-superdf_abtiter%>%
  filter(dpi==dpi_of_interest)%>%
  filter(day.y==day_of_es)%>%
  filter(genesets==genesetsofinterest)%>%
  #cor(x = superdf$es,y = superdf$value,use="everything",method="pearson")%>%
  ggplot(aes(x=es,y=titers_log2,color=animalid,shape=treatment))+
  geom_point()+
  xlab(label = "enrichment score")+
  ylab("log2(titers)")+
  stat_smooth(method = "lm",col="black")+
  #ggtitle(geneofinterest)+
  theme_classic()+
  theme(aspect.ratio = 0.8)
  #facet_wrap(.~treatment)
print(graph_abtiters)
