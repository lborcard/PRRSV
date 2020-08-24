graph_cells<-superdf_clean%>%
  filter(name=="%prol. t-cells")%>%
  filter(day.x=="42")%>%
  #filter(day.y==c("7",""))%>%
  filter(genesets=="M67")%>%
  #cor(x = superdf$es,y = superdf$value,use="everything",method="pearson")%>%
  ggplot(aes(x=es,y=value,shape=stimulus,color=animalid))+
  geom_point()+
  xlab(label = "enrichment score")+
  ylab("%prol. t-cells")+
  stat_smooth(method = "lm",col="black")+
  ggtitle("M67")+
  theme_classic()+
  theme(aspect.ratio = 0.8)+
  facet_grid(day.y~treatment)
print(graph_cells)
