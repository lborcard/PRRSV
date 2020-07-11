
m3<-as.matrix(day7)
m3<-Day7_ES
colnames(m3)<-m3[1,]
rownames(m3)<-m3[,1]
m3<-m3[-1,-1]
mode(m3)<-"integer"
sel<-order(apply(Day7_ES,1,var),decreasing = TRUE)[1:346]
m3_sel<-m3[sel,]
pca<-prcomp(t(m3),scale. = TRUE)
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,3])
pca.data
pca.data$groups <-rep(c("mock","MLV","lp","hp"),each=6)
ggplot(data=pca.data, aes(x=X, y=Y, label=Sample,color=groups)) +
  #geom_text() +
  geom_point()+
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("My PCA Graph day7")
