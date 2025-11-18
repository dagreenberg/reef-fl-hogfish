library(cmdstanr);library(here)
here()
stanc_options = list("O1")
source(here('code','reef_functions.R'))
cmdstanr::set_cmdstan_path(path='C:/Users/greenbergda/Documents/.cmdstan/cmdstan-2.29.2')

#see folder 'stan' for stan model - put into cmdstan_path on local device
file_corrmod<- file.path(cmdstan_path(), "correlated_mod_states.stan")
mod<- cmdstan_model(file_corrmod)

#data####
hogfish<- read.csv(here('data','REEF_hogfish_2024.csv'))
hogfish<- subset(hogfish,geogr4 %notin% c(2301))
hogfish$habitat2<- NA
hogfish<- hogfish%>%
  mutate(
    habitat2=ifelse(habitat%in%c(0,1,11,12,7,8,9),'mixed',habitat2),
    habitat2=ifelse(habitat%in%c(4,5,6),'dropoff',habitat2),
    habitat2=ifelse(habitat%in%c(10),'artificial',habitat2),
    habitat2=ifelse(habitat%in%c(2),'highreef',habitat2),
    habitat2=ifelse(habitat%in%c(3),'lowreef',habitat2)
  )


hogfish$m_year<- paste(hogfish$month,hogfish$year,sep='_')
hogfish=hogfish[complete.cases(hogfish$day),]

#RVC estimates - mean and standard errors of index output
rvc_ests=read.csv(here('data','rvc_model_ests.csv'))

#Map of site overlap####
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

site_n = hogfish %>% group_by(geogr) %>% summarise(n=n())
site_n$geogr4=twa_sites$region.id[match(site_n$geogr,twa_sites$geogid)]
site_n$lat=twa_sites$lat_dd[match(site_n$geogr,twa_sites$geogid)]
site_n$lon=twa_sites$lon_dd[match(site_n$geogr,twa_sites$geogid)]

#RVC
fk_hf_rvc=read.csv(here('data','RVC_hogfish_flkeys_99_24.csv'))
fk_hf_rvc$PSU_YEAR=paste(fk_hf_rvc$PRIMARY_SAMPLE_UNIT,fk_hf_rvc$YEAR)
sefl_hf_rvc=read.csv(here('data','RVC_hogfish_sefl_13_24.csv'))
sefl_hf_rvc$PSU_YEAR=paste(sefl_hf_rvc$PRIMARY_SAMPLE_UNIT,sefl_hf_rvc$YEAR)
dt_hf_rvc=read.csv(here('data','RVC_hogfish_drytort_99_24.csv'))
dt_hf_rvc$PSU_YEAR=paste(dt_hf_rvc$PRIMARY_SAMPLE_UNIT,dt_hf_rvc$YEAR)
hf_rvc=rbind(fk_hf_rvc,sefl_hf_rvc,dt_hf_rvc)
hf_rvc$latlon=paste(hf_rvc$LAT_DEGREES,hf_rvc$LON_DEGREES,sep='')
sefl_hf_rvc$latlon=paste(sefl_hf_rvc$LAT_DEGREES,sefl_hf_rvc$LON_DEGREES,sep='')
fk_hf_rvc$latlon=paste(fk_hf_rvc$LAT_DEGREES,fk_hf_rvc$LON_DEGREES,sep='')
dt_hf_rvc$latlon=paste(dt_hf_rvc$LAT_DEGREES,dt_hf_rvc$LON_DEGREES,sep='')

n=hf_rvc%>%group_by(latlon)%>%summarize(n=n())

rvc_plots_sefl=distinct(sefl_hf_rvc,latlon,.keep_all=T)
rvc_plots_sefl$n=n$n[match(rvc_plots_sefl$latlon,n$latlon)]
rvc_plots_fk=distinct(fk_hf_rvc,latlon,.keep_all=T)
rvc_plots_fk$n=n$n[match(rvc_plots_fk$latlon,n$latlon)]
rvc_plots_dt=distinct(dt_hf_rvc,latlon,.keep_all=T)
rvc_plots_dt$n=n$n[match(rvc_plots_dt$latlon,n$latlon)]

main=ggplot(data = world) + 
  geom_sf(fill= 'antiquewhite') +  xlab("Longitude") + ylab("Latitude") +
  geom_point(data = rvc_plots_sefl, mapping = aes(x = LON_DEGREES, y = LAT_DEGREES), color = adjustcolor('darkcyan',alpha.f=0.1), size = log10(rvc_plots_sefl$n), alpha = 0.1)+ 
  geom_point(data = rvc_plots_fk, mapping = aes(x = LON_DEGREES, y = LAT_DEGREES), color = adjustcolor('royalblue4',alpha.f=0.1), size = log10(rvc_plots_fk$n), alpha = 0.1)+ 
  geom_point(data = rvc_plots_dt, mapping = aes(x = LON_DEGREES, y = LAT_DEGREES), color = adjustcolor('deepskyblue4',alpha.f=0.1), size = log10(rvc_plots_dt$n), alpha = 0.1)+ 
  geom_point(data = site_n, mapping = aes(x = lon, y = lat), color = adjustcolor('darkred',alpha.f=0.1), size = log10(site_n$n), alpha = 0.4)+ 
  coord_sf(xlim = c(-85, -79.5), ylim = c(24, 28), expand = T) +
  annotate(geom = 'point', x = -84.5, y = 27.6, size = log10(1)) +
  annotate(geom = 'point', x = -84.5, y = 27.3, size = log10(5)) +
  annotate(geom = 'point', x = -84.5, y = 27, size = log10(15)) +
  annotate(geom = 'point', x = -84.5, y = 26.7, size = log10(50)) +
  annotate(geom = 'point', x = -84.5, y = 26.4, size = log10(100)) + 
  annotate(geom = 'point', x = -84.5, y = 26.1, size = log10(500)) + 
  annotate(geom = 'point', x = -84.5, y = 25.8, size = log10(2e3)) +
  annotate(geom = 'text', x = -84.5,5, y = 25, label='REEF',col='darkred',hjust = 0) +
  annotate(geom = 'text', x = -84.5,5, y = 24.8, label='RVC - SE Fl.',col='darkcyan',hjust = 0) +
  annotate(geom = 'text', x = -84.5,5, y = 24.6, label='RVC - Fl. Keys',col='royalblue4',hjust = 0) +
  annotate(geom = 'text', x = -84.5,5, y = 24.4, label='RVC - Dry Tort.',col='deepskyblue4',hjust = 0) +
  annotate(geom = 'text', x = -84, y = 27.6, label='1') +
  annotate(geom = 'text', x = -84, y = 27.3, label='5') +
  annotate(geom = 'text', x = -84, y = 27, label='15') +
  annotate(geom = 'text', x = -84, y = 26.7, label='50') +
  annotate(geom = 'text', x = -84, y = 26.4, label='100') + 
  annotate(geom = 'text', x = -84, y = 26.1, label='500') + 
  annotate(geom = 'text', x = -84, y = 25.8, label='2000') +
  annotate(geom = 'text', x = -84.5, y = 27.8,label='no. surveys', size = 4)+
  theme(panel.background = element_rect(fill = 'white'),legend.title = element_blank())

pdf('./outs/site_comparison.pdf',height=6,width=8)
main
dev.off()


#Florida Keys
#Reef - subset data
fk_hogfish=subset(hogfish, geogr4 %notin% c(3410,3301,3302)) #remove sites from SE florida and Dry Tortugas
fk_hogfish$m_year<-  paste(fk_hogfish$year,fk_hogfish$month,sep='_') #create month-year index

#RVC subset data
fk_rvc_ests=subset(rvc_ests,region=='FLKEYS')
fk_rvc_ests=fk_rvc_ests[complete.cases(fk_rvc_ests$mean),]
fk_rvc_ests$year.id=match(fk_rvc_ests[,1],seq(1993,2024)) #include oldest REEF estimates - index by year

#matrix of covariates
X2<- matrix(data=c(scale(as.numeric(fk_hogfish$btime)),scale(as.numeric(fk_hogfish$averagedepth)),scale(as.numeric(fk_hogfish$averagedepth)^2),scale(as.numeric(fk_hogfish$visibility)),scale(as.numeric(fk_hogfish$current)),fk_hogfish$exp_binary),ncol=6,nrow=nrow(fk_hogfish))

#data for model
data_fk=list(y1 = fk_rvc_ests$mean,
             sd_y1 = fk_rvc_ests$sd,
             y2 = fk_hogfish$abundance2,
             N1 = nrow(fk_rvc_ests),
             N2 = nrow(fk_hogfish),
             N_hab2 = length(unique(fk_hogfish$habitat2)),
             hab_class2=as.numeric(factor(fk_hogfish$habitat2)),
             site=as.numeric(factor(fk_hogfish$geogr)),
             N_site=length(unique(fk_hogfish$geogr)),
             diver=as.numeric(factor(fk_hogfish$fish_memberid)),
             N_dv=length(unique(fk_hogfish$fish_memberid)),
             N_mth2=length(unique(fk_hogfish$month)),
             mth2=as.numeric(factor(fk_hogfish$month)),
             dmy=as.numeric(factor(fk_hogfish$site_dmy)),
             N_dmy=length(unique(fk_hogfish$site_dmy)),
             my=as.numeric(factor(fk_hogfish$m_year)),
             N_my=length(unique(fk_hogfish$m_year)),
             K=max(fk_hogfish$abundance2),
             X2=X2,
             Z2=ncol(X2),
             TT=max(fk_hogfish$year)-min(fk_hogfish$year)+1,
             year_id1=fk_rvc_ests$year.id,
             year_id2=as.numeric(fk_hogfish$year-min(fk_hogfish$year)+1))

fit_fk_ss<- mod$sample(
  data = data_fk,
  chains = 5, 
  parallel_chains = 5,
  iter_warmup = 200,
  iter_sampling = 600,
  refresh = 100,
  adapt_delta = 0.999,
  max_treedepth = 20 # print update every 500 iters
)

#save model object & parameter summary
fit_fk_ss$save_object(paste('./outs/rvc_reef_fk_corrmod_',Sys.Date(),'.RDS',sep=''))
write.csv(fit_fk_ss$summary(),paste('./outs/corrmod_fk_parameter_summary_',Sys.Date(),'.csv',sep=''))

#summary of correlation parameter
fit_fk_ss$summary(variables='Cor_t')

#draw main index parameters
params=fit_fk_ss$draws(variables=c('x','est1','a_yr2','cut','Cor_t'),format='draws_matrix')

#timeseries for reef for plot function
hf_reef_ts=ts_reef(fk_hogfish)

#plot scaled indices 
scaled_mv_timeseries_plot2(ts2=hf_reef_ts,params=params,path=here('outs'),TT=data_fk$TT,yrs.rvc=fk_rvc_ests[,1],yr.start=1993,yr.end=2024,sp='Hogfish',reg='FL Keys')

#Table of estimates - RVC estimates from original and augmented model for Florida Keys
Rvc_ests.fk=data.frame(year=seq(1993,2024),original.mean=NA,original.se=NA,state.mean=NA,state.se=NA)
Rvc_ests.fk$original.mean=fk_rvc_ests$mean[match(Rvc_ests.fk$year,fk_rvc_ests[,1])]
Rvc_ests.fk$original.se=fk_rvc_ests$sd[match(Rvc_ests.fk$year,fk_rvc_ests[,1])]
x=as.data.frame(params)[grepl('x',colnames(params))] #process states for RVC and REEF
x1=x[,1:data_fk$TT] #subset to first column of X for just RVC
Rvc_ests.fk$state.mean=apply(x1,2,mean)
Rvc_ests.fk$state.se=apply(x1,2,sd)
write.csv(Rvc_ests.fk,'./outs/florida_keys_indices_table.csv')

#SE florida
sefl_rvc_ests=subset(rvc_ests,region=='SEFL')
sefl_rvc_ests=sefl_rvc_ests[complete.cases(sefl_rvc_ests$mean),]
sefl_rvc_ests$year.id=match(sefl_rvc_ests[,1],seq(1994,2024)) #match to all years covered by REEF

sefl_hogfish=subset(hogfish, geogr4 %in% c(3301,3302)) #remove sites near Tampa & Dry Tortugas
#sefl_hogfish=subset(sefl_hogfish,year>2012)
sefl_hogfish$m_year<-  paste(sefl_hogfish$year,sefl_hogfish$month,sep='_')

X2<- matrix(data=c(scale(as.numeric(sefl_hogfish$btime)),scale(as.numeric(sefl_hogfish$averagedepth)),scale(as.numeric(sefl_hogfish$averagedepth)^2),scale(as.numeric(sefl_hogfish$visibility)),scale(as.numeric(sefl_hogfish$current)),sefl_hogfish$exp_binary),ncol=6,nrow=nrow(sefl_hogfish))

data_sefl=list(y1 = sefl_rvc_ests$mean,
             sd_y1 = sefl_rvc_ests$sd,
             y2 = sefl_hogfish$abundance2,
             N1 = nrow(sefl_rvc_ests),
             N2 = nrow(sefl_hogfish),
             N_hab2 = length(unique(sefl_hogfish$habitat2)),
             hab_class2=as.numeric(factor(sefl_hogfish$habitat2)),
             site=as.numeric(factor(sefl_hogfish$geogr)),
             N_site=length(unique(sefl_hogfish$geogr)),
             diver=as.numeric(factor(sefl_hogfish$fish_memberid)),
             N_dv=length(unique(sefl_hogfish$fish_memberid)),
             N_mth2=length(unique(sefl_hogfish$month)),
             mth2=as.numeric(factor(sefl_hogfish$month)),
             dmy=as.numeric(factor(sefl_hogfish$site_dmy)),
             N_dmy=length(unique(sefl_hogfish$site_dmy)),
             my=as.numeric(factor(sefl_hogfish$m_year)),
             N_my=length(unique(sefl_hogfish$m_year)),
             K=max(sefl_hogfish$abundance2),
             X2=X2,
             Z2=ncol(X2),
             TT=max(sefl_hogfish$year)-min(sefl_hogfish$year)+1,
             year_id1=sefl_rvc_ests$year.id,
             year_id2=as.numeric(sefl_hogfish$year-min(sefl_hogfish$year)+1))

fit_sefl_ss<- mod$sample(
  data = data_sefl,
  seed = 123, 
  chains = 5, 
  parallel_chains = 5,
  iter_warmup = 200,
  iter_sampling = 500,
  refresh = 100,
  adapt_delta = 0.999,
  max_treedepth = 20 # print update every 500 iters
)

#save model fit & parameter sumamry
fit_sefl_ss$save_object(paste('./outs/rvc_reef_sefl_corrmod_',Sys.Date(),'.RDS',sep=''))
write.csv(fit_sefl_ss$summary(),paste('./outs/corrmod_sefl_parameter_summary_',Sys.Date(),'.csv',sep=''))

#summary of correlation parameter
fit_sefl_ss$summary(variables='Cor_t')

#extract key index parameters
params=fit_sefl_ss$draws(variables=c('x','est1','a_yr2','cut','Cor_t'),format='draws_matrix')

#reef timeseries
hf_reef_ts=ts_reef(sefl_hogfish)

#scaled indices plot
scaled_mv_timeseries_plot2(ts2=hf_reef_ts,params=params,path=here('outs'),TT=data_sefl$TT,yrs.rvc=sefl_rvc_ests[,1],yr.start=1994,yr.end=2024,sp='Hogfish',reg='SE Florida')

#Table of estimates - RVC estimates from original and augmented model for Florida Keys
Rvc_ests.sefl=data.frame(year=seq(1994,2024),original.mean=NA,original.se=NA,state.mean=NA,state.se=NA)
Rvc_ests.sefl$original.mean=sefl_rvc_ests$mean[match(Rvc_ests.sefl$year,sefl_rvc_ests[,1])]
Rvc_ests.sefl$original.se=sefl_rvc_ests$sd[match(Rvc_ests.sefl$year,sefl_rvc_ests[,1])]
x=as.data.frame(params)[grepl('x',colnames(params))] #process states for RVC and REEF
x1=x[,1:data_sefl$TT] #subset to first column of X for just RVC
Rvc_ests.sefl$state.mean=apply(x1,2,mean)
Rvc_ests.sefl$state.se=apply(x1,2,sd)
write.csv(Rvc_ests.sefl,'./outs/se_florida_indices_table.csv')
