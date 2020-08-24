##script to open the titers file in excel and clean it
library(tidyverse)
abtiters<-readxl::read_xlsx("ab_titers.xlsx",col_types = "text")

##we remove other measurments from the treatment col because the ab_titers only contain 2 HP and LP

superdf_abtiter<-subset(superdfadj,treatment%in%c("HP","LP"))

unique(superdf_abtiter$treatment)

colnames(abtiters)<-c("treatment","animalid","dpi","titers")

unique(abtiters$titers)

##we change the titers to the format unique digit (i.e 4 8 16)
abtiters$titers_clean<-str_extract(abtiters$titers,"(?<=:)\\d+")

##we change titers to log2 format
abtiters$titers_log2<-log2(as.numeric(abtiters$titers_clean))

##we merge together the 2
superdf_abtiter<-merge(x=superdf_abtiter,y=abtiters)
