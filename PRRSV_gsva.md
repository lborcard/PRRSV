PRRSV-GSVA
================

## Summary

PRRSV is a porcine virus. The idea of the project is to use single
sample gene set enrichment scores to correlate the early innate response
with antibody response at a later timepoint. 24 animal animals were
separated in 4 different groups:

## Analysis

Loading the necessary packages

``` r
## First specify the packages of interest
packages = c("tidyverse",
             "Biobase",
             "limma",
             "GSEABase",
             "sigPathway",
             "readxl",
             "ggpubr")

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
countfile<-read.delim(file = paste0("p470/Day",day,"/extended/p470.counts.txt"),
                      check.names = FALSE,row.names = 1)
##the count file that is read by read.delim depends on the day argument above
```

Small loop to create a genesets list from the GMX file

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


## assign a new name to the data frame with Day + day of the count file + _ES (enrichment score)
```

#### We bind the 3 data frames together

``` r
completeES<-rbind(Day0_ES_t,Day3_ES_t,Day7_ES_t)
head(completeES)
dim(completeES)
```

#### We do a pivot of the Enrichment score (ES) to get only one column with our ES and another with the names of the gensets

``` r
## 
completeES_long<-pivot_longer(completeES,
                              cols = -sampleid,
                              names_to = "genesets",
                              values_to = "es")
```

#### We clean the sample id column

``` r
completeES_long$sampleid<-gsub("_","",completeES_long$sampleid)# we replace the "_" by ""
## we replace "CH431" by "" in order to use the ID number for merging the data set
completeES_long$sampleid<-gsub("CH431","",completeES_long$sampleid)
samplesinfos$sampleid<-gsub("CH431-","",samplesinfos$sampleid)
## we set the sampleid col to numeric in the ES and the sampleinfos file
completeES_long$sampleid<-as.numeric(completeES_long$sampleid)
```

#### Import the data frame containing the info about the samples

``` r
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

samplesinfos$sampleid<-as.numeric(samplesinfos$sampleid)

#######Major step we merge the ES data with the sample informations ####

completeES_long_merge<-merge(completeES_long,samplesinfos)
```

We now have a dataframe with the infos of the samples with the
enrichmentscore values We proceed with cleaning of the data

``` r
completeES_long_merge$day<-gsub("Day","",completeES_long_merge$day)

#we remove the trailing whitespace
completeES_long_merge$day<-str_trim(completeES_long_merge$day)
```

``` r
## We transform to factors
completeES_long_merge$genesets<-as.factor(as.character(completeES_long_merge$genesets))
str(completeES_long_merge)
```

We take the proliferation data sets which contains the flow cytometry
data and do some cleaning

``` r
proliferation <- read_excel("T cell data/Kick_PRRSV_proliferation.xlsx", 
                            skip = 1)
proliferation_longer<-pivot_longer(proliferation,cols=-1:-5)

colnames(proliferation_longer)<-tolower(colnames(proliferation_longer))

proliferation_longer$name<-as.factor(tolower(proliferation_longer$name))

names(proliferation_longer)<-c("well","day","treatment","stimulus","animalid","name","value" )

proliferation_longer$animalid<-factor(proliferation_longer$animalid)
proliferation_longer$day<-factor(as.character(proliferation_longer$day))

str(proliferation_longer)
head(proliferation_longer)
##We change the values of the stimulus (1-3-4=LP,1-7-3 =HP,MED=MOCK,MLV= vr2234)
proliferation_longer$stimulus<-as.factor(proliferation_longer$stimulus)
levels(proliferation_longer$stimulus)<- c("LP","HP","MOCK","MLV")
```

We calculate the adjusted value for each value of proliferation, i.e we
subtract the value of the stimulation with “MED” (MOCK) to the
homologous stimulation

``` r
test_valueadj<-as.data.frame(proliferation_longer) 

test_valueadj$value<-replace_na(test_valueadj$value,replace = 0)

test_test_valueadj_MOCK <- filter(test_valueadj, stimulus == "MOCK")

test_test_valueadj_NONMOCK <- filter(test_valueadj, Group != "MOCK") 

test_valueadj <- left_join(
        test_valueadj,
        test_test_valueadj_MOCK,
        by = c("animal", "Day", "name","Group"),
        suffix = c("", "_mock")
) %>%
        mutate(adjusted = value - value_mock)

test_valueadj<-filter(test_valueadj,Group==stimulus)

test_valueadj$valueadj<-ifelse(sign(test_valueadj$adjusted)==-1,0,test_valueadj$adjusted)

proliferation_longer_adjusted<-test_valueadj[-c(3,7,8,9)]

colnames(proliferation_longer_adjusted)<- c("day","treatment","animalid","name","value","valueadj")
head(proliferation_longer_adjusted)
```

``` r
#we clean certain columns 

completeES_long_merge$treatment<-word(completeES_long_merge$treatment,1)#we extract the treatment value in a simpler form now ("MLV","MOCK"...)
completeES_long_merge$day<-factor(as.character(completeES_long_merge$day))
completeES_long_merge$animalid<-factor(completeES_long_merge$animalid)
```

    #for the last step we merge together the data sets ES and proliferation data in 1 super data frame called superdf
    
    superdf<-merge(proliferation_longer,
                   completeES_long_merge,
                   by=c("animalid","treatment"),all = TRUE)

head(superdf) str(superdf) \#\#create a superdf using the ratio between
the days7 and day3 /day0 \#\# we proliferation\_clean which contains
only the homologous treatment MOCK MOCK HP HP etc
superdf\_ratio\_clean\<-merge(proliferation\_longer\_clean,completeES\_long\_merge\_ratio,by=c(“animalid”),all
= TRUE) \`\`\`
