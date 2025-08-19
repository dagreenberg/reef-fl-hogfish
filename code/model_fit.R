library(cmdstanr);library(loo);library(bayesplot);library(dplyr);library(tidyverse);library(stringr);library(lubridate);library(rlist)
library(here)
here()
stanc_options = list("O1")
source(here('code','reef_functions.R'))
cmdstanr::set_cmdstan_path(path='C:/Users/greenbergda/Documents/.cmdstan/cmdstan-2.29.2')

#REEF model alone####
hogfish<- read.csv(here('data','REEF_hogfish_2024.csv'))
hogfish$habitat2<- NA
hogfish<- hogfish%>%
  mutate(
    habitat2=ifelse(habitat%in%c(0,1,11,12,7,8,9),'mixed',habitat2),
    habitat2=ifelse(habitat%in%c(4,5,6),'dropoff',habitat2),
    habitat2=ifelse(habitat%in%c(10),'artificial',habitat2),
    habitat2=ifelse(habitat%in%c(2),'highreef',habitat2),
    habitat2=ifelse(habitat%in%c(3),'lowreef',habitat2)
  )


hogfish$month_id<- paste(hogfish$year,hogfish$month,sep='_')
hogfish$m_year<- paste(hogfish$fish_memberid,hogfish$month,hogfish$year,sep='_')
hogfish=hogfish[complete.cases(hogfish$day),]

site_year=expand.grid(levels(factor(hogfish$geogr)),levels(factor(hogfish$year)))
site_year= site_year[order(site_year[,1],site_year[,2]),]
site_year[,3]=paste(site_year[,1],site_year[,2],sep="_")

dv_year=expand.grid(unique(hogfish$fish_memberid),unique(hogfish$year))
dv_year= dv_year[order(dv_year[,1]),]
dv_year[,3]=paste(dv_year[,1],dv_year[,2],sep="_")

hogfish$site_year_id<- match(paste(hogfish$geogr,hogfish$year,sep='_'),site_year[,3]) #this properly indexes each year-region combination
hogfish$dv_year_id<- match(paste(hogfish$fish_memberid,hogfish$year,sep='_'),dv_year[,3]) #this properly indexes each year-region combination

X<- matrix(data=c(scale(as.numeric(hogfish$btime)),scale(as.numeric(hogfish$averagedepth)),scale(as.numeric(hogfish$averagedepth)^2),scale(as.numeric(hogfish$visibility)),scale(as.numeric(hogfish$current)),hogfish$exp_binary),ncol=6,nrow=nrow(hogfish))

file_test_varsite<- file.path(cmdstan_path(), "gg models","site_var_mod_simple.stan")
mod_sv<- cmdstan_model(file_test_varsite)


data_2.1=list(y=hogfish$abundance2,
              N = nrow(hogfish),
              dmy=as.numeric(factor(hogfish$site_dmy)),
              N_dmy=length(unique(hogfish$site_dmy)),
              site=as.numeric(factor(hogfish$geogr)),
              N_site=length(unique(hogfish$geogr)),
              hab=as.numeric(factor(hogfish$habitat2)),
              N_hab=length(unique(hogfish$habitat2)),
              mth=as.numeric(factor(hogfish$month)),
              N_my=length(unique(hogfish$month_id)),
              my=as.numeric(factor(hogfish$month_id)),
              N_mth=length(unique(hogfish$month)),
              diver=as.numeric(factor(hogfish$fish_memberid)),
              N_dv=length(unique(hogfish$fish_memberid)),
              K=length(unique(hogfish$abundance2)),
              X=X,
              Z=ncol(X),
              TT=max(hogfish$year)-min(hogfish$year)+1,
              N_yr=length(unique(hogfish$year)),
              yr_index=sort(unique(as.numeric(factor(hogfish$year)))),
              year_id=as.numeric(factor(hogfish$year)),
              site_year_id=hogfish$site_year_id,
              N_R=1,
              dv_year_id=hogfish$dv_year_id)

fit2.1 <- mod_sv$sample(
  data = data_2.1,
  seed = 123, 
  chains = 4, 
  parallel_chains = 4,
  iter_warmup = 200,
  iter_sampling = 600,
  refresh = 100,
  adapt_delta = 0.99,
  max_treedepth = 20 # print update every 500 iters
)

params_gg_fl<- fit2.1$draws(format='draws_matrix',variables = c('cut','x','a_yr'))
gg_ts_fl<- ts_reef(hogfish)
State_space_timeseries_plot(sp='Hogfish',GZ=c('Florida'),params=params_gg_fl,TT=32,ts=gg_ts_fl,cols = c('turquoise4','royalblue4'),n.iter = 2400)

library(ggspatial);library(sf);library(rnaturalearth)


twa_sites<- read.csv(here('data','TWAgeog.csv'))
twa_sites$region.id<- substr(twa_sites$geogid, 1,4) #Get the reion id (first four digits)
twa_sites$lat_full<- twa_sites$lat #Copy of the full latitude - these are in degrees, minutes (with seconds as a fraction of minutes)
twa_sites$lon_full<- twa_sites$lon #Copy of the full longitude

twa_sites<-twa_sites%>% #Separate out degrees and minutes
  tidyr::separate(lat,into=c("lat_deg","lat_min"),sep=" ")%>%
  tidyr::separate(lon,into=c("lon_deg","lon_min"),sep=" ")
col.num<-c("lat_deg","lat_min","lon_deg","lon_min")
twa_sites[col.num]<-sapply(twa_sites[col.num],as.numeric) #will get warnings for aberrant entries
twa_sites$lat_dd<- twa_sites$lat_deg+twa_sites$lat_min/60
twa_sites$lon_dd<- twa_sites$lon_deg-twa_sites$lon_min/60
twa_sites<- twa_sites[complete.cases(twa_sites$lat_dd),] #remove sites/regions without spatial coordinates
twa_sites<- twa_sites[complete.cases(twa_sites$lon_dd),] #remove sites/regions without spatial coordinates

world <- ne_countries(scale = "medium", returnclass = "sf")

site_n = hogfish %>% group_by(geogr) %>% summarize(n=n())
site_n$lat=twa_sites$lat_dd[match(site_n$geogr,twa_sites$geogid)]
site_n$lon=twa_sites$lon_dd[match(site_n$geogr,twa_sites$geogid)]

main=ggplot(data = world) + 
  geom_sf(fill= 'antiquewhite') +  xlab("Longitude") + ylab("Latitude") +
  geom_point(data = site_n, mapping = aes(x = lon, y = lat), color = 'royalblue4', size = 2*log10(site_n$n), alpha = 0.7)+ 
  coord_sf(xlim = c(-85, -79), ylim = c(24, 28), expand = T) +
  annotate(geom = 'point', x = -84.5, y = 27, size = 2*log10(50)) +
  annotate(geom = 'point', x = -84.5, y = 26.5, size = 2*log10(100)) + 
  annotate(geom = 'point', x = -84.5, y = 26, size = 2*log10(500)) + 
  annotate(geom = 'point', x = -84.5, y = 25.5, size = 2*log10(2e3)) +
  annotate(geom = 'text', x = -84, y = 27, label='50') +
  annotate(geom = 'text', x = -84, y = 26.5, label='100') + 
  annotate(geom = 'text', x = -84, y = 26, label='500') + 
  annotate(geom = 'text', x = -84, y = 25.5, label='2000') +
  annotate(geom = 'text', x = -84.5, y = 27.5,label='n. surveys', size = 4)+
  theme(panel.background = element_rect(fill = 'white'),legend.title = element_blank())

main

#REEF-RVC combined models by region####

#seperate FL keys & SE and Dry Tortugas
fk_hogfish=subset(hogfish, geogr4 %notin% c(2301,3410)) #remove sites near Tampa & Dry Tortugas
fk_hogfish=subset(fk_hogfish, year>1998) #sync with start to RVC data (1999)

dt_hogfish=subset(hogfish, geogr4==3410)
dt_hogfish=subset(dt_hogfish, year>1998) #sync with start to RVC data (1999)

fk_hogfish$month_id<- paste(fk_hogfish$year,fk_hogfish$month,sep='_')
fk_hogfish$m_year<- paste(fk_hogfish$fish_memberid,fk_hogfish$month,fk_hogfish$year,sep='_')
fk_hogfish=fk_hogfish[complete.cases(fk_hogfish$day),]
fk_hogfish$R_id=rep(1,nrow(fk_hogfish))

dt_hogfish$month_id<- paste(dt_hogfish$year,dt_hogfish$month,sep='_')
dt_hogfish$m_year<- paste(dt_hogfish$fish_memberid,dt_hogfish$month,dt_hogfish$year,sep='_')
dt_hogfish=dt_hogfish[complete.cases(dt_hogfish$day),]
dt_hogfish$R_id=rep(2,nrow(dt_hogfish))
hogfish=rbind(fk_hogfish,dt_hogfish)

#RVC datasets
fk_hf_rvc=read.csv(here('data','RVC_hogfish_flkeys_99_24.csv'))
fk_hf_rvc$PSU_YEAR=paste(fk_hf_rvc$PRIMARY_SAMPLE_UNIT,fk_hf_rvc$YEAR)
sefl_hf_rvc=read.csv(here('data','RVC_hogfish_sefl_13_24.csv'))
sefl_hf_rvc$PSU_YEAR=paste(sefl_hf_rvc$PRIMARY_SAMPLE_UNIT,sefl_hf_rvc$YEAR)
dt_hf_rvc=read.csv(here('data','RVC_hogfish_drytort_99_24.csv'))
dt_hf_rvc$PSU_YEAR=paste(dt_hf_rvc$PRIMARY_SAMPLE_UNIT,dt_hf_rvc$YEAR)
fk_hf_rvc$REGION_ID=rep(1,nrow(fk_hf_rvc))
sefl_hf_rvc$REGION_ID=rep(1,nrow(sefl_hf_rvc))
dt_hf_rvc$REGION_ID=rep(2,nrow(dt_hf_rvc))
hf_rvc=rbind(fk_hf_rvc,sefl_hf_rvc,dt_hf_rvc)
hf_rvc$OCC=ifelse(hf_rvc$NUM2>0,1,0)
fk_hf_rvc=subset(hf_rvc,REGION_ID==1)
dt_hf_rvc=subset(hf_rvc,REGION_ID==2)


#Overview - RVC
sample_freq_fk_rvc=fk_hf_rvc %>%group_by(YEAR) %>% summarize(rvc.n.SSU=n(),rvc.n.PSU=n_distinct(PRIMARY_SAMPLE_UNIT),rvc.freq.occ=sum(OCC)/n())
sample_freq_dt_rvc=dt_hf_rvc %>%group_by(YEAR) %>% summarize(rvc.n.SSU=n(),rvc.n.PSU=n_distinct(PRIMARY_SAMPLE_UNIT),rvc.freq.occ=sum(OCC)/n())

#Overview - REEF
sample_freq_fk_reef=fk_hogfish %>%group_by(year) %>% summarize(reef.n.surveys=n(),reef.n.sites=n_distinct(geogr),reef.freq.occ=sum(occ)/n())
sample_freq_dt_reef=dt_hogfish %>%group_by(year) %>% summarize(reef.n.surveys=n(),reef.n.sites=n_distinct(geogr),reef.freq.occ=sum(occ)/n())

colnames(sample_freq_fk_rvc)[1]='year'
colnames(sample_freq_dt_rvc)[1]='year'

sample_freqs_fk=full_join(sample_freq_fk_rvc, sample_freq_fk_reef, by = "year")
sample_freqs_fk
sample_freqs_dt=full_join(sample_freq_dt_rvc, sample_freq_dt_reef, by = "year")

write.csv(sample_freqs_fk,'sampling_frequencies_hogfish_FlKeys.csv')
write.csv(sample_freqs_dt,'sampling_frequencies_hogfish_DryTort.csv')


#region by year
reg.year=expand.grid(levels(factor(hf_rvc$REGION_ID)),seq(min(hf_rvc$YEAR),max(hf_rvc$YEAR)))
reg.year[,3]=paste(reg.year[,1],reg.year[,2],sep='.')
reg.year=reg.year[order(reg.year[,1],reg.year[,2]),]

hf_rvc$REGION_YR_Id=match(paste(hf_rvc$REGION_ID,hf_rvc$YEAR,sep='.'),reg.year[,3])
hogfish$R_yr_id=match(paste(hogfish$R_id,hogfish$year,sep='.'),reg.year[,3])

#Combined model
X1<- matrix(data=c(scale(as.numeric(hf_rvc$DEPTH)),scale(as.numeric(hf_rvc$DEPTH^2))),ncol=2,nrow=nrow(hf_rvc))
X2<- matrix(data=c(scale(as.numeric(hogfish$btime)),scale(as.numeric(hogfish$averagedepth)),scale(as.numeric(hogfish$averagedepth)^2),scale(as.numeric(hogfish$visibility)),scale(as.numeric(hogfish$current)),hogfish$exp_binary),ncol=6,nrow=nrow(hogfish))

data_fk=list(y1 = hf_rvc$NUM2,
          y2 = hogfish$abundance2,
          N1 = nrow(hf_rvc),
          N2 = nrow(hogfish),
          N_psu = length(unique(hf_rvc$PSU_YEAR)),
          psu_yr = as.numeric(factor(hf_rvc$PSU_YEAR)),
          N_hab1 = length(unique(hf_rvc$HABITAT_CD)),
          hab_class1=as.numeric(factor(hf_rvc$HABITAT_CD)),
          N_strat1=length(unique(hf_rvc$STRAT)),
          stratum1=as.numeric(factor(hf_rvc$STRAT)),
          N_hab2 = length(unique(hogfish$habitat2)),
          hab_class2=as.numeric(factor(hogfish$habitat2)),
          site=as.numeric(factor(hogfish$geogr)),
          N_site=length(unique(hogfish$geogr)),
          diver=as.numeric(factor(hogfish$fish_memberid)),
          N_dv=length(unique(hogfish$fish_memberid)),
          N_mth1=length(unique(hf_rvc$MONTH)),
          N_mth2=length(unique(hogfish$month)),
          mth1=as.numeric(factor(hf_rvc$MONTH)),
          mth2=as.numeric(factor(hogfish$month)),
          dmy=as.numeric(factor(hogfish$site_dmy)),
          N_dmy=length(unique(hogfish$site_dmy)),
          my=as.numeric(factor(hogfish$m_year)),
          N_my=length(unique(hogfish$m_year)),
          K=max(hogfish$abundance2),
          X1=X1,
          Z1=ncol(X1),
          X2=X2,
          Z2=ncol(X2),
          TT=max(hf_rvc$YEAR)-min(hf_rvc$YEAR)+1,
          N_yr1=length(unique(hf_rvc$YEAR)),
          yr_index1=sort(unique(as.numeric(factor(hf_rvc$YEAR)))),
          year_id1=as.numeric(factor(hf_rvc$YEAR)),
          N_yr2=length(unique(hogfish$year)),
          yr_index2=sort(unique(as.numeric(factor(hogfish$year)))),
          year_id2=as.numeric(factor(hogfish$year)),
          R=2,
          R_year_id1=hf_rvc$REGION_YR_Id,
          R_year_id2=hogfish$R_yr_id)

file_sharedmod<- file.path(cmdstan_path(), "REEF models","shared_state_REEF_RVC_mod.stan")
mod_ss<- cmdstan_model(file_sharedmod)

fit_fk_ss<- mod_ss$sample(
  data = data_fk,
  seed = 123, 
  chains = 5, 
  parallel_chains = 5,
  iter_warmup = 200,
  iter_sampling = 400,
  refresh = 100,
  adapt_delta = 0.995,
  max_treedepth = 20 # print update every 500 iters
)


#DT model
X1<- matrix(data=c(scale(as.numeric(dt_hf_rvc$DEPTH)),scale(as.numeric(dt_hf_rvc$DEPTH^2))),ncol=2,nrow=nrow(dt_hf_rvc))
X2<- matrix(data=c(scale(as.numeric(dt_hogfish$btime)),scale(as.numeric(dt_hogfish$averagedepth)),scale(as.numeric(dt_hogfish$averagedepth)^2),scale(as.numeric(dt_hogfish$visibility)),scale(as.numeric(dt_hogfish$current)),dt_hogfish$exp_binary),ncol=6,nrow=nrow(dt_hogfish))

data_fk=list(y1 = dt_hf_rvc$NUM2,
             y2 = dt_hogfish$abundance2,
             N1 = nrow(dt_hf_rvc),
             N2 = nrow(dt_hogfish),
             N_psu = length(unique(dt_hf_rvc$PSU_YEAR)),
             psu_yr = as.numeric(factor(dt_hf_rvc$PSU_YEAR)),
             N_hab1 = length(unique(dt_hf_rvc$HABITAT_CD)),
             hab_class1=as.numeric(factor(dt_hf_rvc$HABITAT_CD)),
             N_strat1=length(unique(dt_hf_rvc$STRAT)),
             stratum1=as.numeric(factor(dt_hf_rvc$STRAT)),
             N_hab2 = length(unique(dt_hogfish$habitat2)),
             hab_class2=as.numeric(factor(dt_hogfish$habitat2)),
             site=as.numeric(factor(dt_hogfish$geogr)),
             N_site=length(unique(dt_hogfish$geogr)),
             diver=as.numeric(factor(dt_hogfish$fish_memberid)),
             N_dv=length(unique(dt_hogfish$fish_memberid)),
             N_mth1=length(unique(dt_hf_rvc$MONTH)),
             N_mth2=length(unique(dt_hogfish$month)),
             mth1=as.numeric(factor(dt_hf_rvc$MONTH)),
             mth2=as.numeric(factor(dt_hogfish$month)),
             dmy=as.numeric(factor(dt_hogfish$site_dmy)),
             N_dmy=length(unique(dt_hogfish$site_dmy)),
             my=as.numeric(factor(dt_hogfish$m_year)),
             N_my=length(unique(dt_hogfish$m_year)),
             K=max(dt_hogfish$abundance2),
             X1=X1,
             Z1=ncol(X1),
             X2=X2,
             Z2=ncol(X2),
             TT=max(dt_hf_rvc$YEAR)-min(dt_hf_rvc$YEAR)+1,
             N_yr1=length(unique(dt_hf_rvc$YEAR)),
             yr_index1=sort(unique(as.numeric(factor(dt_hf_rvc$YEAR)))),
             year_id1=as.numeric(factor(dt_hf_rvc$YEAR)),
             N_yr2=length(unique(dt_hogfish$year)),
             yr_index2=sort(unique(as.numeric(factor(dt_hogfish$year)))),
             year_id2=as.numeric(factor(dt_hogfish$year)))

file_sharedmod<- file.path(cmdstan_path(), "REEF models","shared_state_REEF_RVC_mod.stan")
mod_ss<- cmdstan_model(file_sharedmod)

fit_dt_ss<- mod_ss$sample(
  data = data_fk,
  seed = 123, 
  chains = 6, 
  parallel_chains = 6,
  iter_warmup = 200,
  iter_sampling = 500,
  refresh = 100,
  adapt_delta = 0.995,
  max_treedepth = 20 # print update every 500 iters
)


#correlated states model####
#RVC data
fk_hf_rvc=read.csv(here('data','RVC_hogfish_flkeys_99_24.csv'))
fk_hf_rvc$PSU_YEAR=paste(fk_hf_rvc$PRIMARY_SAMPLE_UNIT,fk_hf_rvc$YEAR)
sefl_hf_rvc=read.csv(here('data','RVC_hogfish_sefl_13_24.csv'))
sefl_hf_rvc$PSU_YEAR=paste(sefl_hf_rvc$PRIMARY_SAMPLE_UNIT,sefl_hf_rvc$YEAR)

sefk_hf_rvc=rbind(fk_hf_rvc,sefl_hf_rvc)

dt_hf_rvc=read.csv(here('data','RVC_hogfish_drytort_99_24.csv'))
dt_hf_rvc$PSU_YEAR=paste(dt_hf_rvc$PRIMARY_SAMPLE_UNIT,dt_hf_rvc$YEAR)

#REEF data
fk_hogfish=subset(hogfish, geogr4 %notin% c(2301,3410)) #remove sites near Tampa & Dry Tortugas
fk_hogfish=subset(fk_hogfish, year>1998) #sync with start to RVC data (1999)

dt_hogfish=subset(hogfish, geogr4==3410)
dt_hogfish=subset(dt_hogfish, year>1998) #sync with start to RVC data (1999)

fk_hogfish$month_id<- paste(fk_hogfish$year,fk_hogfish$month,sep='_')
fk_hogfish$m_year<- paste(fk_hogfish$fish_memberid,fk_hogfish$month,fk_hogfish$year,sep='_')
fk_hogfish=fk_hogfish[complete.cases(fk_hogfish$day),]

dt_hogfish$month_id<- paste(dt_hogfish$year,dt_hogfish$month,sep='_')
dt_hogfish$m_year<- paste(dt_hogfish$fish_memberid,dt_hogfish$month,dt_hogfish$year,sep='_')
dt_hogfish=dt_hogfish[complete.cases(dt_hogfish$day),]

#SE Florida & Keys

X1<- matrix(data=c(scale(as.numeric(sefk_hf_rvc$DEPTH)),scale(as.numeric(sefk_hf_rvc$DEPTH^2)),scale(as.numeric(sefk_hf_rvc$UNDERWATER_VISIBILITY)),sefk_hf_rvc$PROT),ncol=4,nrow=nrow(sefk_hf_rvc))
X2<- matrix(data=c(scale(as.numeric(fk_hogfish$btime)),scale(as.numeric(fk_hogfish$averagedepth)),scale(as.numeric(fk_hogfish$averagedepth)^2),scale(as.numeric(fk_hogfish$visibility)),scale(as.numeric(fk_hogfish$current)),fk_hogfish$exp_binary),ncol=6,nrow=nrow(fk_hogfish))

data_fk=list(y1 = sefk_hf_rvc$NUM2,
             y2 = fk_hogfish$abundance2,
             N1 = nrow(sefk_hf_rvc),
             N2 = nrow(fk_hogfish),
             N_psu = length(unique(sefk_hf_rvc$PSU_YEAR)),
             psu_yr = as.numeric(factor(sefk_hf_rvc$PSU_YEAR)),
             N_hab1 = length(unique(sefk_hf_rvc$HABITAT_CD)),
             hab_class1=as.numeric(factor(sefk_hf_rvc$HABITAT_CD)),
             N_strat1=length(unique(sefk_hf_rvc$STRAT)),
             stratum1=as.numeric(factor(sefk_hf_rvc$STRAT)),
             N_hab2 = length(unique(fk_hogfish$habitat2)),
             hab_class2=as.numeric(factor(fk_hogfish$habitat2)),
             site=as.numeric(factor(fk_hogfish$geogr)),
             N_site=length(unique(fk_hogfish$geogr)),
             diver=as.numeric(factor(fk_hogfish$fish_memberid)),
             N_dv=length(unique(fk_hogfish$fish_memberid)),
             N_mth1=length(unique(sefk_hf_rvc$MONTH)),
             N_mth2=length(unique(fk_hogfish$month)),
             mth1=as.numeric(factor(sefk_hf_rvc$MONTH)),
             mth2=as.numeric(factor(fk_hogfish$month)),
             dmy=as.numeric(factor(fk_hogfish$site_dmy)),
             N_dmy=length(unique(fk_hogfish$site_dmy)),
             my=as.numeric(factor(fk_hogfish$m_year)),
             N_my=length(unique(fk_hogfish$m_year)),
             K=max(fk_hogfish$abundance2),
             X1=X1,
             Z1=ncol(X1),
             X2=X2,
             Z2=ncol(X2),
             TT=max(sefk_hf_rvc$YEAR)-min(sefk_hf_rvc$YEAR)+1,
             year_id1=as.numeric(sefk_hf_rvc$YEAR-min(sefk_hf_rvc$YEAR)+1),
             year_id2=as.numeric(fk_hogfish$year-min(sefk_hf_rvc$YEAR)+1))

file_corrmod<- file.path(cmdstan_path(), "REEF models","correlated_states_REEF_RVC_mod.stan")
mod<- cmdstan_model(file_corrmod)

fit_fk_ss<- mod$sample(
  data = data_fk,
  seed = 123, 
  chains = 4, 
  parallel_chains = 4,
  iter_warmup = 200,
  iter_sampling = 500,
  refresh = 100,
  adapt_delta = 0.995,
  max_treedepth = 20 # print update every 500 iters
)

fit_fk_ss$save_object('./outs/rvc_reef_sefk_corrmod.RDS')
write.csv(fit_fk_ss$summary(),'./outs/corrmod_parameter_summary.csv')

fit_fk_ss=readRDS('./outs/rvc_reef_sefk_corrmod.RDS')

hf_reef_ts=ts_reef(fk_hogfish)
sefk_hf_rvc$OCC=ifelse(sefk_hf_rvc$NUM2>0,1,0)
hf_rvc_ts=ts_rvc(sefk_hf_rvc)

params=fit_fk_ss$draws(variables=c('x','a_yr1','a_yr2','cut','Cor_t'),format='draws_matrix')
  
scaled_mv_timeseries_plot(ts1=hf_rvc_ts,ts2=hf_reef_ts,params=params,path=here('outs'),TT=data_fk$TT,TT.rvc=nrow(hf_rvc_ts[complete.cases(hf_rvc_ts$n.occ),]),yr.start=1999,yr.end=2024,sp='Hogfish',reg='Southeast & Keys')

X1<- matrix(data=c(scale(as.numeric(dt_hf_rvc$DEPTH)),scale(as.numeric(dt_hf_rvc$DEPTH^2)),scale(as.numeric(dt_hf_rvc$UNDERWATER_VISIBILITY)),dt_hf_rvc$PROT),ncol=4,nrow=nrow(dt_hf_rvc))
X2<- matrix(data=c(scale(as.numeric(dt_hogfish$btime)),scale(as.numeric(dt_hogfish$averagedepth)),scale(as.numeric(dt_hogfish$averagedepth)^2),scale(as.numeric(dt_hogfish$visibility)),scale(as.numeric(dt_hogfish$current)),dt_hogfish$exp_binary),ncol=6,nrow=nrow(dt_hogfish))

data_dt=list(y1 = dt_hf_rvc$NUM2,
             y2 = dt_hogfish$abundance2,
             N1 = nrow(dt_hf_rvc),
             N2 = nrow(dt_hogfish),
             N_psu = length(unique(dt_hf_rvc$PSU_YEAR)),
             psu_yr = as.numeric(factor(dt_hf_rvc$PSU_YEAR)),
             N_hab1 = length(unique(dt_hf_rvc$HABITAT_CD)),
             hab_class1=as.numeric(factor(dt_hf_rvc$HABITAT_CD)),
             N_strat1=length(unique(dt_hf_rvc$STRAT)),
             stratum1=as.numeric(factor(dt_hf_rvc$STRAT)),
             N_hab2 = length(unique(dt_hogfish$habitat2)),
             hab_class2=as.numeric(factor(dt_hogfish$habitat2)),
             site=as.numeric(factor(dt_hogfish$geogr)),
             N_site=length(unique(dt_hogfish$geogr)),
             diver=as.numeric(factor(dt_hogfish$fish_memberid)),
             N_dv=length(unique(dt_hogfish$fish_memberid)),
             N_mth1=length(unique(dt_hf_rvc$MONTH)),
             N_mth2=length(unique(dt_hogfish$month)),
             mth1=as.numeric(factor(dt_hf_rvc$MONTH)),
             mth2=as.numeric(factor(dt_hogfish$month)),
             dmy=as.numeric(factor(dt_hogfish$site_dmy)),
             N_dmy=length(unique(dt_hogfish$site_dmy)),
             my=as.numeric(factor(dt_hogfish$m_year)),
             N_my=length(unique(dt_hogfish$m_year)),
             K=max(dt_hogfish$abundance2),
             X1=X1,
             Z1=ncol(X1),
             X2=X2,
             Z2=ncol(X2),
             TT=max(dt_hf_rvc$YEAR)-min(dt_hf_rvc$YEAR)+1,
             year_id1=as.numeric(dt_hf_rvc$YEAR-min(dt_hf_rvc$YEAR)+1),
             year_id2=as.numeric(dt_hogfish$year-min(dt_hf_rvc$YEAR)+1))

file_corrmod<- file.path(cmdstan_path(), "REEF models","correlated_states_REEF_RVC_mod.stan")
mod<- cmdstan_model(file_corrmod)

fit_dt_ss<- mod$sample(
  data = data_dt,
  seed = 123, 
  chains = 4, 
  parallel_chains = 4,
  iter_warmup = 200,
  iter_sampling = 500,
  refresh = 100,
  adapt_delta = 0.995,
  max_treedepth = 20 # print update every 500 iters
)

fit_dt_ss$save_object('./outs/rvc_reef_dt_corrmod.RDS')
write.csv(fit_dt_ss$summary(),'./outs/corrmod_dt_parameter_summary.csv')

fit_dt_ss=readRDS('./outs/rvc_reef_dt_corrmod.RDS')

hf_reef_ts=ts_reef(dt_hogfish)
dt_hf_rvc$OCC=ifelse(dt_hf_rvc$NUM2>0,1,0)
hf_rvc_ts=ts_rvc(dt_hf_rvc)

params=fit_dt_ss$draws(variables=c('x','a_yr1','a_yr2','cut','Cor_t'),format='draws_matrix')

scaled_mv_timeseries_plot(ts1=hf_rvc_ts,ts2=hf_reef_ts,params=params,path=here('outs'),TT=data_dt$TT,TT.rvc=nrow(hf_rvc_ts[complete.cases(hf_rvc_ts$n.occ),]),yr.start=1999,yr.end=2024,sp='Hogfish',reg='Dry Tortugas')

