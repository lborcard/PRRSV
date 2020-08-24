test<-superdfadj %>% group_by(day.y)%>% nest()
  #mutate(es3=filter(superdfadj,day.y=="3")$es/filter(superdfadj,day.y=="0")$es)

test_wider<-pivot_longer(completeES,
                              cols = -sampleid,
                              names_to = "genesets",
                              values_to = "es")
completeES_long_merge_wider<-pivot_wider(completeES_long_merge,
                                         id_cols = c(genesets,animalid),
                              names_from = day,
                              names_prefix = "es",
                              values_from = es)

g15<-(order(apply(filter(completeES_long_merge_wider,animalid=="G15")[,-c(1,2)],1,var),decreasing = TRUE)[1:100])

g14<-order(apply(filter(completeES_long_merge_wider,animalid=="G14")[,-c(1,2)],1,var),decreasing = TRUE)[1:100]

completeES_long_merge_wider$es3_0<-completeES_long_merge_wider$es3/completeES_long_merge_wider$es0

completeES_long_merge_wider$es7_0<-completeES_long_merge_wider$es7/completeES_long_merge_wider$es0

completeES_long_merge_wider$es7_3<-completeES_long_merge_wider$es7/completeES_long_merge_wider$es3

completeES_long_merge_wider_ratio<-completeES_long_merge_wider[,-c(3,4,5)]

completeES_long_merge_ratio<-pivot_longer(completeES_long_merge_wider_ratio,-c(genesets,animalid),names_prefix = "es",names_to = "day",values_to = "es")

##construction of the dataframes containing the antibodies inhibition
#Lp antibodies
lp_ab_titer <- read_excel("ab_titers.xlsx",
                          sheet = "LP_pigs", range = "h1:m13")

colnames(lp_ab_titer)<-lp_ab_titer[1,]
lp_ab_titer<-lp_ab_titer[-c(1),]
lp_ab_titer$day<-c(3,7,10,14,21,28,35,42,49,56,63)
lp_ab_titer<-pivot_longer(lp_ab_titer,cols = -c(day),names_to = "animalid",values_to = "ab_inhib")

## HP titers
hp_ab_titer <- read_excel("ab_titers.xlsx",
                          sheet = "HP_pigs", range = "N1:S13")

colnames(hp_ab_titer)<-hp_ab_titer[1,]
hp_ab_titer<-hp_ab_titer[-c(1),]
hp_ab_titer$day<-c(3,7,10,14,21,28,35,42,49,56,63)
hp_ab_titer<-pivot_longer(hp_ab_titer,cols = -c(day),names_to = "animalid",values_to = "ab_inhib")

##MLV titer
mlv_ab_titer <- read_excel("ab_titers.xlsx",
                          sheet = "MLV_pigs", range = "B1:G13")

colnames(mlv_ab_titer)<-mlv_ab_titer[1,]
mlv_ab_titer<-mlv_ab_titer[-c(1),]
mlv_ab_titer$day<-c(3,7,10,14,21,28,35,42,49,56,63)
mlv_ab_titer<-pivot_longer(mlv_ab_titer,cols = -c(day),names_to = "animalid",values_to = "ab_inhib")
