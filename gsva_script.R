library(Biobase)
library(dplyr)
library(limma)
library(GSEABase)
library(GSVA)
library(sigPathway)
library(readxl)
library(tidyverse)
library(ggpubr)

## 
day <- 7
#countfile<-read.delim(file = "p470/Day0/extended/p470.counts.txt",check.names = FALSE,row.names = 1)
#the count file that is read by read.delim depends on the day argument above
countfile<-read.delim(file = paste0("p470/Day",day,"/extended/p470.counts.txt"),check.names = FALSE,row.names = 1)

#the countfile must be read as a matrix to be used in gsva function

eset_matrix <- as.matrix(countfile)

btmgeneset<-gmxToG("BTM gene sets.gmx")# from the sigpathway package

#small loop to create a genesets list to label them in the dataframe
for (i in 1:346){
        
        genesets[i]<-btmgeneset[[i]]$title
}

results_gsva<-lapply(X = btmgeneset,FUN = function(x){gsva(expr=eset_matrix, 
     gset.idx.list=x, 
     annotation,
     method="gsva",#we use the gsva method for the Enrichment scores
     kcdf="Poisson",#poisson for discrete numbers
     abs.ranking=FALSE,
     min.sz=1,#minimal size of each genesets
     max.sz=Inf,#maximal size of each genesets
     parallel.sz=0,
     parallel.type="SOCK",
     mx.diff=TRUE,
     tau=1,#tau=switch(method,gsva=1, ssgsea=0.25, NA),
     ssgsea.norm=TRUE,
     verbose=TRUE)})

results_gsva_df<-data.frame(matrix(unlist(results_gsva),#we transform the list in a dataframe
                                   nrow = length(results_gsva),
                                   byrow = T),
                            stringsAsFactors = FALSE)

results_gsva_df<-t(results_gsva_df)#we transpose because we want the genesets as col and samples as rows

sampleid<-colnames(countfile)#we take the names of the samples from the countfile

colnames(results_gsva_df)<- genesets #the genesets names are taken from the genesets vector created earlier

#add the rownames of the samples (now as rows)
results_gsva_df$sampleid<-sampleid

#assign a new name to the data frame with Day + day of the count file + _ES (enrichment score)
assign(paste("Day",day,"_ES",sep = ""),results_gsva_df)

######################################################################################################################################
######################################################################################################################################
###################At this point the 3 data frames  with the enrichment score should be ready#######################
######################################################################################################################################
######################################################################################################################################

#we bind the 3 data frames by rows
completeES<-rbind(Day0_ES_t,Day3_ES_t,Day7_ES_t)

#make a samplid column containing the infos about the  samples taken from the rownames
#completeES$sampleid <- row.names(completeES)

#import the data frame containing the info about the samples
samplesinfos <- read_excel("CH431.xlsx", 
                           col_names = FALSE, skip = 2,col_types = c("text"))
#we assign new names for the columns
colnames(samplesinfos)<- c("sampleid",
                           "virus",
                           "day",
                           "code",
                           "animalid",
                           "treatment",
                           "rqn")
#we do a pivot of the Enrichment score (ES) to get only one column with our ES
completeES_long<-pivot_longer(completeES,
                              cols = -sampleid,
                              names_to = "genesets",
                              values_to = "es")

# we clean the sample id column
completeES_long$sampleid<-gsub("_","",completeES_long$sampleid)# we replace the "_" by ""
## we replace "CH431" by "" in order to use the ID number for merging the data set
completeES_long$sampleid<-gsub("CH431","",completeES_long$sampleid)
samplesinfos$sampleid<-gsub("CH431-","",samplesinfos$sampleid)
## we set the sampleid col to numeric in the ES and the sampleinfos file
completeES_long$sampleid<-as.numeric(completeES_long$sampleid)

samplesinfos$sampleid<-as.numeric(samplesinfos$sampleid)

#######Major step we merge the ES data with the sample informations ####

completeES_long_merge<-merge(completeES_long,samplesinfos)

completeES_long_merge$day<-gsub("Day","",completeES_long_merge$day)

#we remove the trailing whitespace
completeES_long_merge$day<-str_trim(completeES_long_merge$day)

#we transform to factors
completeES_long_merge$genesets<-as.factor(as.character(completeES_long_merge$genesets))
str(completeES_long_merge)

#we take the proliferation data sets and do some cleaning

library(readxl)
proliferation <- read_excel("T cell data/Kick_PRRSV_proliferation.xlsx", 
                            skip = 1)
proliferation_longer<-pivot_longer(proliferation,cols=-1:-5)
colnames(proliferation_longer)<-tolower(colnames(proliferation_longer))
proliferation_longer$name<-as.factor(tolower(proliferation_longer$name))
names(proliferation_longer)<-c("well","day","treatment","stimulus","animalid","name","value" )

#we clean certain columns 

completeES_long_merge$treatment<-word(completeES_long_merge$treatment,1)#we extract the treatment value in a simpler form now ("MLV","MOCK"...)
completeES_long_merge$day<-factor(as.character(completeES_long_merge$day))
completeES_long_merge$animalid<-factor(completeES_long_merge$animalid)
proliferation_longer$animalid<-factor(proliferation_longer$animalid)
proliferation_longer$day<-factor(as.character(proliferation_longer$day))

#for the last step we merge together the data sets ES and proliferation data in 1 super data frame called superdf

superdf<-merge(proliferation_longer,completeES_long_merge,by=c("animalid","treatment"),all = TRUE)

##some graphs to explore the data 
graph2<-superdf%>%
        filter(name=="%prol. t-cells")%>%
        filter(day.x=="28")%>%
        filter(day.y=="3")%>%
        filter(genesets=="M0")%>%
        ggplot(aes(x=es,y=value,shape=stimulus,color=animalid))+
        geom_point()+
        ylab("%prol. t-cells")+
        xlab(label = "enrichment score")+
        theme_classic()+
        stat_smooth(method = "lm",col="black")+
        theme(aspect.ratio = 0.8)+
        facet_wrap(.~treatment)
graph3<-superdf_clean%>%
        filter(name=="%prol. t-cells")%>%
        filter(day.x=="42")%>%
        filter(day.y==c("7"))%>%
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
        facet_wrap(.~treatment)
print(graph3)
graph4<- superdf%>%
        filter(name=="%prol. treg" )%>%
        ggplot(aes(x=day.x,y=value),color=treatment)+
        geom_boxplot()+
        facet_wrap(~treatment)
print(graph4)
graphs<-ggarrange(graph2,graph3,graph4,common.legend = TRUE,labels = list("Day56vsDay7","Day28vsDay7"))
print(graphs)        
        
