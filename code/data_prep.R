rm(list=ls())
library(dplyr)
library(here)
here()

source(here('code','reef_functions.R'))

#data filtering for REEF data not included - due to file size for repository - see REEF_hogfish_2024.csv for filtered data for hogfish0

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
