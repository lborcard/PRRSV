PRRSV-GSVA
================

## Summary

PRRSV is a porcine virus. The idea of the project is to use single
sample gene set enrichment scores to correlate the early innate response
with antibody response at a later timepoint. 24 animal animals were
separated in 4 different groups:

Loading the necessary packages

``` r
## First specify the packages of interest
packages = c("tidyverse","Biobase","limma","GSEABase","sigPathway","readxl","ggpubr")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      BiocManager::install(x)
      library(x, character.only = TRUE)
    }
  }
)
```

## 

we can change the value for the day

day \<- 7

``` r
countfile<-read.delim(file = paste0("p470/Day",day,"/extended/p470.counts.txt"),check.names = FALSE,row.names = 1)
##the count file that is read by read.delim depends on the day argument above
```

small loop to create a genesets list

### Main function of the GSVA that will generate the Enrichment Scores

``` r
results_gsva<-lapply(X = btmgeneset,FUN = function(x){gsva(expr=eset_matrix, 
     gset.idx.list=x, 
     annotation,
     method="gsva",
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
```

``` r
results_gsva_df<-data.frame(matrix(unlist(results_gsva),#we transform the list in a dataframe
                                   nrow = length(results_gsva),
                                   byrow = T),
                            stringsAsFactors = FALSE)

results_gsva_df<-t(results_gsva_df)#we transpose because we want the genesets as col and samples as rows

sampleid<-colnames(countfile)#we take the names of the samples

colnames(results_gsva_df)<- genesets

#add the rownames of the samples (now as rows)


##assign a new name to the data frame with Day + day of the count file + _ES (enrichment score)
```

We bind the 3 data frames togther

``` r
#we bind the 3 data frames by rows
completeES<-rbind(Day0_ES_t,Day3_ES_t,Day7_ES_t)
head(completeES)
dim(completeES)
```

``` r
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
```

## Including Plots

You can also embed plots, for example:

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.
