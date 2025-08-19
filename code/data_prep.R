rm(list=ls())
library(dplyr);library(tidyverse);library(stringr);library(lubridate);library(rlist)
library(here)
here()

source(here('code','reef_functions.R'))

###
fish<- read.csv(here('data','fish.csv'))
#colnames(fish)=c('id','formid','fish_memberid','region','site','speciesid','familyid','abundance','timecode','site4','site1','comment')
R<- subset(fish,region=='TWA')
surveys<- read.csv(here('data','TWAsurveys.csv'))

R[,13:31]<- surveys[match(R$formid,surveys$formid),2:20]

# let's get the year & month from the Date field (requires lubridate)
R$date<-ymd(R$date) #put into proper date format - note there are many with no dates (0000-00-00) so these will be dropped
R<-cbind(R,year=year(R$date))
R<-cbind(R,month=month(R$date))
R<-cbind(R,day=day(R$date))

# Thin  raw data to exclude dives shorter than 20min, longer than 120min
R<-R[as.numeric(R$btime)>20,]
R<-R[as.numeric(R$btime)<120,]

# Thin the data to remove night dives (start before 5am or after 8pm)
R<-R[as.numeric(R$start)>5,]
R<-R[as.numeric(R$start)<20,]

write.csv(R,here('data','R.csv'))

R<- read.csv(here('data','R.csv'))
rm(fish) #save space
gc()
hogfish<- reef_filter_sp(R=R,GZ=c(2101,2201,2301,3101,3201,3301,3403:3411),sp=0130,invert=0)

hogfish$habitat2<- NA
hogfish<- hogfish%>%
  mutate(
    habitat2=ifelse(habitat%in%c(0,1,11,12,7,8,9),'mixed',habitat2),
    habitat2=ifelse(habitat%in%c(4,5,6),'dropoff',habitat2),
    habitat2=ifelse(habitat%in%c(10),'artificial',habitat2),
    habitat2=ifelse(habitat%in%c(2),'highreef',habitat2),
    habitat2=ifelse(habitat%in%c(3),'lowreef',habitat2)
)

write.csv(hogfish,here('data','REEF_hogfish_2024.csv'))

#RVC data
devtools::install_github('jeremiaheb/rvc')

keys99_24 = rvc::getRvcData(years = 1999:2024, regions = c("FLA KEYS"))
x=keys99_24$sample_data
x$SSU_YEAR= paste(x$YEAR,x$SEC_SAMPLE_UNIT,x$STATION_NR,sep='_')
x1= x %>% select(SSU_YEAR,SPECIES_CD,everything())
x2= complete(x1,SSU_YEAR,nesting(SPECIES_CD),fill=list(NUM=0)) #ensures all non-sightings are recorded
zeros = anti_join(x2,x1)


sefls13_24 = rvc::getRvcData(years = 2013:2024, regions = c("SEFCRI"))

#species code = LAC MAXI (hogfish) 
hf_rvc_keys_99_24=subset(keys99_24$sample_data,SPECIES_CD=='LAC MAXI')
hf_rvc_keys_99_24$NUM2=base::ceiling(hf_rvc_keys_99_24$NUM) #round all fractions upwards to nearest whole number for negative binomial integers

hf_rvc_sefl_13_24=subset(sefls13_24$sample_data,SPECIES_CD=='LAC MAXI')
hf_rvc_sefl_13_24$NUM2=base::ceiling(hf_rvc_sefl_13_24$NUM) #round all fractions upwards to nearest whole number for negative binomial integers
hf_rvc_sefl_13_24=hf_rvc_sefl_13_24[,-12] #remove RUGOSITY_CD column
hf_rvc_sefl_13_24=hf_rvc_sefl_13_24[,-22] #remove sample_year and DEPTH_STRAT columns
hf_rvc_sefl_13_24=hf_rvc_sefl_13_24[,-22] #remove sample_year and DEPTH_STRAT columns

tort99_24 = rvc::getRvcData(years = 1999:2024, regions = c("DRY TORT"))

hf_rvc_drytort_99_24=subset(tort99_24$sample_data,SPECIES_CD=='LAC MAXI')
hf_rvc_drytort_99_24$NUM2=base::ceiling(hf_rvc_drytort_99_24$NUM) #round all fractions upwards to nearest whole number for negative binomial integers

write.csv(hf_rvc_keys_99_24,here('data','RVC_hogfish_flkeys_99_24.csv'))
write.csv(hf_rvc_sefl_13_24,here('data','RVC_hogfish_sefl_13_24.csv'))
write.csv(hf_rvc_drytort_99_24,here('data','RVC_hogfish_drytort_99_24.csv'))
