library(cmdstanr);library(loo);library(bayesplot);library(dplyr);library(tidyverse);library(stringr);library(lubridate);library(rlist)
library(here)
here()
stanc_options = list("O1")
source(here('code','reef_functions.R'))
cmdstanr::set_cmdstan_path(path='C:/Users/greenbergda/Documents/.cmdstan/cmdstan-2.29.2')

#load data####
hogfish<- read.csv(here('Hogfish_2024.csv'))
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
