graph_cells<-function(genesets_graph="M123",day_facs_graph="28",treatment=c("HP","LP"),cells_graph="%prol. th")
  {superdf_minus%>%
  filter(name==cells_graph)%>%
  filter(day_facs==day_facs_graph)%>%
  #filter(day.y==c("7",""))%>%
  filter(genesets==genesets_graph)%>%
  filter(treatment%in%c("HP","LP"))%>%
  #cor(x = superdf$es,y = superdf$value,use="everything",method="pearson")%>%
  ggplot(aes(x=es_min,y=value,color=animalid))+
  geom_point()+
  xlab(label = "enrichment score")+
  ylab(cells_graph)+
  stat_smooth(method = "lm",col="black")+
  ggtitle(genesets_graph)+
  theme_classic()+
  theme(aspect.ratio = 0.8)+
  facet_grid(day_min~.)
  }
print(graph_cells)

graph_auc<-superdf_auc%>%
  filter(celltype=="%ctl")%>%
  filter(genesets=="M152.2")%>%
  #filter(treatment=="LP")%>%
  #cor(x = superdf$es,y = superdf$value,use="everything",method="pearson")%>%
  ggplot(aes(x=es_min,y=auc_cells,shape=day_min,color=animalid))+
  geom_point()+
  xlab(label = "enrichment score")+
  ylab("%prol. t-cells")+
  stat_smooth(method = "lm",col="black")+
  ggtitle("M152.2")+
  theme_classic()+
  theme(aspect.ratio = 0.8)+
  facet_grid(day_min~treatment)
print(graph_auc)
