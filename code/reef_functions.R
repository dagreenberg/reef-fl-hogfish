'%notin%' <- Negate(`%in%`)

reef_filter_sp = function(R,GZ,sp,invert){ #function to trim dataframe and add in occurrence data by species
  occ_list<- list()
  TempDat<-R %>% subset(geogr4 %in% GZ) %>% select('formid','speciesid','abundance',everything())
  if(invert==0){
    TempDat<- subset(TempDat,type!=2)
  }
  if(invert==1){
    TempDat<- subset(TempDat,type!=1)
  }
  if(invert==2){
    TempDat<- subset(TempDat,type!=1)
  }
  TempDat_sp<- TempDat %>% subset(speciesid==sp)
  TempDat_site<- TempDat_sp %>% group_by(geogr) %>% summarize(n.year=length(unique(year)),n=n()) %>% subset(n.year>4&n>=20)
  
  TempDat2<- subset(TempDat, geogr %in% TempDat_site$geogr)
  
  Zeros<- tidyr::complete(TempDat2,formid,tidyr::nesting(speciesid),
                   fill=list(abundance=0)) %>%
    anti_join(TempDat2) #Creates a dataset where each species is represented in every survey (adding in zero observations)
  
  m<- match(Zeros$formid,TempDat2$formid) #Match the zero data frame to the rest of the data
  Zeros[,4:ncol(Zeros)]<- TempDat2[m,4:ncol(TempDat2)] #Replicate the survey-level data (except abundance)
  Zeros<- subset(Zeros,speciesid==sp)
  TempDat2<- subset(TempDat2,speciesid==sp)
  TempDat3<- rbind(TempDat2,Zeros) %>% arrange(formid,speciesid) #Combine presence and absence data for all species
  TempDat3$occ=ifelse(TempDat3$abundance>0,1,0) #Code presence/absence based on abundance
  TempDat3<- TempDat3 %>% mutate(exp_binary=ifelse(exp=='E',1,0))
  TempDat3$abundance2<- TempDat3$abundance+1 #For modelling abundance categories, has to start at 1
  TempDat3$abund_trans <-NA #Transform abundance categories into the minimum counts
  TempDat3<- TempDat3%>%
    mutate(
      abund_trans=ifelse(abundance==0,0,abund_trans),
      abund_trans=ifelse(abundance==1,1,abund_trans),
      abund_trans=ifelse(abundance==2,2,abund_trans),
      abund_trans=ifelse(abundance==3,11,abund_trans),
      abund_trans=ifelse(abundance==4,101,abund_trans)
    )
 
   
  surveyors<- TempDat3 %>% group_by(fish_memberid) %>% summarize(n=length(unique(formid))) #Determine survey effort of each diver in the region
  surveyors_trim<- subset(surveyors,n>=9) #only keep surveys by members with 10 or more dives
  
  TempDat4<- subset(TempDat3, fish_memberid %in% surveyors_trim$fish_memberid) #Trim out surveys from surveys with less than 5
  
  TempDat_site2<- TempDat4 %>% group_by(geogr) %>% summarize(n=n(),n.yr=n_distinct(year)) %>% subset(n>49&n.yr>3)
  TempDat5<- subset(TempDat4, geogr %in% TempDat_site2$geogr)

  survs<- TempDat5 %>% distinct(formid,.keep_all = T) %>% group_by(geogr,day,month,year) %>% summarize(n=n()) %>% mutate(site_dmy=paste(geogr,day,month,year,sep='')) %>% subset(n>1)
  TempDat6<- TempDat5 %>% mutate(site_dmy=ifelse(paste(geogr,day,month,year,sep='') %in% survs$site_dmy,paste(geogr,day,month,year,sep=''),'baseline'))
  
  return(TempDat6)
}

reef_filter_geog = function(R,geog,sp,invert){ #function to trim dataframe and add in occurrence data by species
  occ_list<- list()
  TempDat<-R %>% subset(geogr %in% geog) %>% select('formid','speciesid','abundance',everything())
  if(invert==0){
    TempDat<- subset(TempDat,type!=2)
  }
  if(invert==1){
    TempDat<- subset(TempDat,type!=1)
  }
  if(invert==2){
    TempDat<- subset(TempDat,type!=1)
  }
  TempDat_sp<- TempDat %>% subset(speciesid==sp)
  TempDat_site<- TempDat_sp %>% group_by(geogr) %>% summarize(n=n(),n.year=length(unique(year))) %>% subset(n.year>2)
  
  TempDat2<- subset(TempDat, geogr %in% TempDat_site$geogr)
  Zeros<- tidyr::complete(TempDat2,formid,nesting(speciesid),
                   fill=list(abundance=0)) %>%
    anti_join(TempDat) #Creates a dataset where each species is represented in every survey (adding in zero observations)
  m<- match(Zeros$formid,TempDat2$formid) #Match the zero data frame to the rest of the data
  Zeros[,4:ncol(Zeros)]<- TempDat2[m,4:ncol(TempDat2)] #Replicate the survey-level data (except abundance)
  TempDat3<- rbind(TempDat2,Zeros) %>% arrange(formid,speciesid) #Combine presence and absence data for all species
  TempDat3$occ=ifelse(TempDat3$abundance>0,1,0) #Code presence/absence based on abundance
  TempDat3<- TempDat3 %>% mutate(exp_binary=ifelse(exp=='E',1,0))
  TempDat3$abundance2<- TempDat3$abundance+1 #For modelling abundance categories, has to start at 1
  TempDat3$abund_trans <-NA #Transform abundance categories into the minimum counts
  TempDat3<- TempDat3%>%
    mutate(
      abund_trans=ifelse(abundance==0,0,abund_trans),
      abund_trans=ifelse(abundance==1,1,abund_trans),
      abund_trans=ifelse(abundance==2,2,abund_trans),
      abund_trans=ifelse(abundance==3,11,abund_trans),
      abund_trans=ifelse(abundance==4,101,abund_trans)
    )
  surveyors<- TempDat3 %>% group_by(fish_memberid) %>% summarize(n=length(unique(formid))) #Determine survey effort of each diver in the region
 # surveyors_trim<- subset(surveyors,n>=3) #only keep surveys by members with 5 or more dives
  
#  TempDat4<- subset(TempDat3, fish_memberid %in% surveyors_trim$fish_memberid) #Trim out surveys from surveys with less than 5
  site_subsets<-TempDat3 %>% group_by(geogr) %>% summarize(n=n_distinct(formid)) %>% subset(n>=3) #Calculate surveys per site
  
  TempDat5<- TempDat3 %>% subset(geogr %in% site_subsets$geogr) #Only keep well surveyed sites (n>=5)
  survs<- TempDat5 %>% distinct(formid,.keep_all = T) %>% group_by(geogr,day,month,year) %>% summarize(n=n()) %>% mutate(site_dmy=paste(geogr,day,month,year,sep='')) %>% subset(n>1)
  TempDat5<- TempDat5 %>% mutate(site_dmy=ifelse(paste(geogr,day,month,year,sep='') %in% survs$site_dmy,paste(geogr,day,month,year,sep=''),'baseline'))
  
  
  occ_dat<- subset(TempDat5,speciesid==sp) #Subset out each species in the provided dataframe
  return(occ_dat)
}


State_space_timeseries_plot<- function(sp,GZ,params,TT,ts,cols,n.iter){
  #  pdf(paste('./figures/',paste(sp,GZ,'state_space',sep='_'),'.pdf',sep=''),width=8,height=6)
  
  cuts=as.data.frame(params)[grepl('cut',colnames(params))];
  x=as.data.frame(params)[grepl('x',colnames(params))];#remove x0
  y=as.data.frame(params)[grepl('a_yr',colnames(params))]; 
  
  lambda_mat<- list()  
  for(i in 1:TT){
    reef_coef<- data.frame(p_0.y=NA,p_1.y=NA,p_2.y=NA,p_11.y=NA,p_101.y=NA,p_0.x=NA,p_1.x=NA,p_2.x=NA,p_11.x=NA,p_101.x=NA,lambda.y=NA,lambda.x=NA,n.iter=seq(1,n.iter))
    
    if(ncol(cuts)==2){
        reef_coef[,1]<- plogis(cuts[,1]-y[,i])
        reef_coef[,2]<-plogis(cuts[,2]-y[,i])-plogis(cuts[,1]-y[,i])
        reef_coef[,3]<-1-plogis(cuts[,2]-y[,i])
        reef_coef[,4]<- 0
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(cuts[,1]-(x[,i]))
        reef_coef[,7]<-plogis(cuts[,2]-(x[,i]))-plogis(cuts[,1]-(x[,i]))
        reef_coef[,8]<-1-plogis(cuts[,2]-(x[,i]))
        reef_coef[,9]<- 0
        reef_coef[,10]<- 0
        reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    if(ncol(cuts)==3){
        reef_coef[,1]<- plogis(cuts[,1]-y[,i])
        reef_coef[,2]<-plogis(cuts[,2]-y[,i])-plogis(cuts[,1]-y[,i])
        reef_coef[,3]<-plogis(cuts[,3]-y[,i])-plogis(cuts[,2]-y[,i])
        reef_coef[,4]<- 1-plogis(cuts[,3]-y[,i])
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(cuts[,1]-(x[,i]))
        reef_coef[,7]<-plogis(cuts[,2]-(x[,i]))-plogis(cuts[,1]-(x[,i]))
        reef_coef[,8]<-plogis(cuts[,3]-(x[,i]))-plogis(cuts[,2]-(x[,i]))
        reef_coef[,9]<- 1-plogis(cuts[,3]-(x[,i]))
        reef_coef[,10]<- 0
        reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    
    if(ncol(cuts)==4){
        reef_coef[,1]<- plogis(cuts[,1]-y[,i])
        reef_coef[,2]<-plogis(cuts[,2]-y[,i])-plogis(cuts[,1]-y[,i])
        reef_coef[,3]<-plogis(cuts[,3]-y[,i])-plogis(cuts[,2]-y[,i])
        reef_coef[,4]<- plogis(cuts[,4]-y[,i])-plogis(cuts[,3]-y[,i])
        reef_coef[,5]<- 1- plogis(cuts[,4]-y[,i])
        reef_coef[,6]=plogis(cuts[,1]-(x[,i]))
        reef_coef[,7]<-plogis(cuts[,2]-(x[,i]))-plogis(cuts[,1]-(x[,i]))
        reef_coef[,8]<-plogis(cuts[,3]-(x[,i]))-plogis(cuts[,2]-(x[,i]))
        reef_coef[,9]<- plogis(cuts[,4]-(x[,i]))-plogis(cuts[,3]-(x[,i]))
        reef_coef[,10]<- 1-plogis(cuts[,4]-(x[,i]))
       reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    lambda_mat[[i]]=reef_coef
    
  }  
#  median.x=median(do.call(rbind, lapply(lambda_mat, data.frame, stringsAsFactors=FALSE))$lambda.x)
 # median.y=median(do.call(rbind, lapply(lambda_mat, data.frame, stringsAsFactors=FALSE))$lambda.y)
  
  x_mat<- data.frame(median.reef=NA,l.95.rf=NA,u.95.rf=NA)
  for(i in 1:TT){
    x_mat[i,1]=median(lambda_mat[[i]]$lambda.x)
    x_mat[i,2]=quantile(lambda_mat[[i]]$lambda.x,0.1)
    x_mat[i,3]=quantile(lambda_mat[[i]]$lambda.x,0.9)
  }
  
  y_mat<- data.frame(year=seq(min(ts$year),max(ts$year)),median.reef=NA,l.95.rf=NA,u.95.rf=NA)
  
  for(i in 1:TT){
    y_mat[i,2]=median(lambda_mat[[i]]$lambda.y)
    y_mat[i,3]=quantile(lambda_mat[[i]]$lambda.y,0.1)
    y_mat[i,4]=quantile(lambda_mat[[i]]$lambda.y,0.9)
  }
  
  par(xpd=T)
  plot(y_mat$median.reef~y_mat$year,type='n',ylim=c(min(x_mat),max(c(max(na.omit(y_mat[,2])),max(x_mat)))),col='darkblue',bty='l',ylab=expression('Mean count per survey'),xlab='Year',main=paste(sp,GZ,sep=', '))
  lines(x_mat[,1]~y_mat$year,lty=5,lwd=2,col=cols[1])
  x<- c(y_mat$year, rev(y_mat$year))
  y1<- c(x_mat[,2], rev(x_mat[,3]))
  polygon(x, y1, col = adjustcolor(cols[1], alpha = 0.2), border=NA) # Add uncertainty polygon
  
  lines(y_mat$median.reef~y_mat$year,col=cols[2],lwd=2)
  points(y_mat$median.reef~y_mat$year,col='white',pch=21,bg=cols[2],cex=1.5)
  text(y=rep(min(x_mat),nrow(ts)),x=ts$year,ts$n.survs,cex=0.5)
  #  dev.off(paste(paste(sp,GZ,sep='_'),'.pdf',sep=''))
  #  dev.off()
}

State_space_timeseries_plot_rev<- function(sp,GZ,params,TT,ts,cols,n.iter){
  #  pdf(paste('./figures/',paste(sp,GZ,'state_space',sep='_'),'.pdf',sep=''),width=8,height=6)
  
  cuts=as.data.frame(params)[grepl('cut',colnames(params))];
  x=as.data.frame(params)[grepl('x',colnames(params))];#remove x0
  y=as.data.frame(params)[grepl('a_yr',colnames(params))]; 
  
  lambda_mat<- list()  
  for(i in 1:TT){
    reef_coef<- data.frame(p_0.y=NA,p_1.y=NA,p_2.y=NA,p_11.y=NA,p_101.y=NA,p_0.x=NA,p_1.x=NA,p_2.x=NA,p_11.x=NA,p_101.x=NA,lambda.y=NA,lambda.x=NA,n.iter=seq(1,n.iter))
    
    if(ncol(cuts)==2){
      reef_coef[,1]<- plogis(cuts[,1]-y[,i])
      reef_coef[,2]<-plogis(cuts[,2]-y[,i])-plogis(cuts[,1]-y[,i])
      reef_coef[,3]<-1-plogis(cuts[,2]-y[,i])
      reef_coef[,4]<- 0
      reef_coef[,5]<- 0
      reef_coef[,6]=plogis(cuts[,1]-(x[,i]))
      reef_coef[,7]<-plogis(cuts[,2]-(x[,i]))-plogis(cuts[,1]-(x[,i]))
      reef_coef[,8]<-1-plogis(cuts[,2]-(x[,i]))
      reef_coef[,9]<- 0
      reef_coef[,10]<- 0
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    if(ncol(cuts)==3){
      reef_coef[,1]<- plogis(cuts[,1]-y[,i])
      reef_coef[,2]<-plogis(cuts[,2]-y[,i])-plogis(cuts[,1]-y[,i])
      reef_coef[,3]<-plogis(cuts[,3]-y[,i])-plogis(cuts[,2]-y[,i])
      reef_coef[,4]<- 1-plogis(cuts[,3]-y[,i])
      reef_coef[,5]<- 0
      reef_coef[,6]=plogis(cuts[,1]-(x[,i]))
      reef_coef[,7]<-plogis(cuts[,2]-(x[,i]))-plogis(cuts[,1]-(x[,i]))
      reef_coef[,8]<-plogis(cuts[,3]-(x[,i]))-plogis(cuts[,2]-(x[,i]))
      reef_coef[,9]<- 1-plogis(cuts[,3]-(x[,i]))
      reef_coef[,10]<- 0
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    
    if(ncol(cuts)==4){
      reef_coef[,1]<- plogis(cuts[,1]-y[,i])
      reef_coef[,2]<-plogis(cuts[,2]-y[,i])-plogis(cuts[,1]-y[,i])
      reef_coef[,3]<-plogis(cuts[,3]-y[,i])-plogis(cuts[,2]-y[,i])
      reef_coef[,4]<- plogis(cuts[,4]-y[,i])-plogis(cuts[,3]-y[,i])
      reef_coef[,5]<- 1- plogis(cuts[,4]-y[,i])
      reef_coef[,6]=plogis(cuts[,1]-(x[,i]))
      reef_coef[,7]<-plogis(cuts[,2]-(x[,i]))-plogis(cuts[,1]-(x[,i]))
      reef_coef[,8]<-plogis(cuts[,3]-(x[,i]))-plogis(cuts[,2]-(x[,i]))
      reef_coef[,9]<- plogis(cuts[,4]-(x[,i]))-plogis(cuts[,3]-(x[,i]))
      reef_coef[,10]<- 1-plogis(cuts[,4]-(x[,i]))
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    lambda_mat[[i]]=reef_coef
    
  }  
  #  median.x=median(do.call(rbind, lapply(lambda_mat, data.frame, stringsAsFactors=FALSE))$lambda.x)
  # median.y=median(do.call(rbind, lapply(lambda_mat, data.frame, stringsAsFactors=FALSE))$lambda.y)
  
  x_mat<- data.frame(median.reef=NA,l.95.rf=NA,u.95.rf=NA)
  for(i in 1:TT){
    x_mat[i,1]=median(lambda_mat[[i]]$lambda.x)
    x_mat[i,2]=quantile(lambda_mat[[i]]$lambda.x,0.2)
    x_mat[i,3]=quantile(lambda_mat[[i]]$lambda.x,0.8)
  }
  
  y_mat<- data.frame(year=seq(min(ts$year),max(ts$year)),median.reef=NA,l.95.rf=NA,u.95.rf=NA)
  
  for(i in 1:TT){
    y_mat[i,2]=median(lambda_mat[[i]]$lambda.y)
    y_mat[i,3]=quantile(lambda_mat[[i]]$lambda.y,0.2)
    y_mat[i,4]=quantile(lambda_mat[[i]]$lambda.y,0.8)
  }
  
  par(xpd=T)
  plot(x_mat$median.reef~y_mat$year,type='n',ylim=c(min(x_mat),max(c(max(na.omit(x_mat[,2])),max(x_mat)))),xlim=c(min(y_mat$year),max(y_mat$year)),col='darkblue',bty='l',ylab=expression('Mean count per survey'),xlab='Year')
  lines(x_mat[,1]~y_mat$year,lwd=2,col=cols[1])
  x<- c(y_mat$year, rev(y_mat$year))
  y1<- c(x_mat[,2], rev(x_mat[,3]))
  polygon(x, y1, col = adjustcolor(cols[1], alpha = 0.2), border=NA) # Add uncertainty polygon
  
  points(x_mat$median.reef~y_mat$year,col='white',pch=21,bg=cols[2],cex=1.5)
#  text(y=rep(min(x_mat),nrow(ts)),x=ts$year,ts$n.survs,cex=0.5)
  #  dev.off(paste(paste(sp,GZ,sep='_'),'.pdf',sep=''))
  #  dev.off()
}

State_space_timeseries_plot_alt2<- function(sp,GZ,params,TT,ts,cols,n.iter){
  #  pdf(paste('./figures/',paste(sp,GZ,'state_space',sep='_'),'.pdf',sep=''),width=8,height=6)
  
  cuts=as.data.frame(params)[grepl('cut',colnames(params))];
  x=as.data.frame(params)[grepl('x_mu',colnames(params))];#mean among all sites
  xs=as.data.frame(params)[grepl('x_site_t',colnames(params))];#site specific trajectories
  
  lambda_mat<- list()  
  for(i in 1:TT){
    reef_coef<- data.frame(p_0.x=NA,p_1.x=NA,p_2.x=NA,p_11.x=NA,p_101.x=NA,lambda.x=NA,n.iter=seq(1,n.iter))
    
    if(ncol(cuts)==2){
      reef_coef[,1]=plogis(cuts[,1]-(x[,i]))
      reef_coef[,2]<-plogis(cuts[,2]-(x[,i]))-plogis(cuts[,1]-(x[,i]))
      reef_coef[,3]<-1-plogis(cuts[,2]-(x[,i]))
      reef_coef[,4]<- 0
      reef_coef[,5]<- 0
      reef_coef[,6]<- apply(reef_coef[,1:5],1,abund_tranfs)
    }
    if(ncol(cuts)==3){
      reef_coef[,1]=plogis(cuts[,1]-(x[,i]))
      reef_coef[,2]<-plogis(cuts[,2]-(x[,i]))-plogis(cuts[,1]-(x[,i]))
      reef_coef[,3]<-plogis(cuts[,3]-(x[,i]))-plogis(cuts[,2]-(x[,i]))
      reef_coef[,4]<- 1-plogis(cuts[,3]-(x[,i]))
      reef_coef[,5]<- 0
      reef_coef[,6]<- apply(reef_coef[,1:5],1,abund_tranfs)
    }
    
    if(ncol(cuts)==4){
      reef_coef[,1]=plogis(cuts[,1]-(x[,i]))
      reef_coef[,2]<-plogis(cuts[,2]-(x[,i]))-plogis(cuts[,1]-(x[,i]))
      reef_coef[,3]<-plogis(cuts[,3]-(x[,i]))-plogis(cuts[,2]-(x[,i]))
      reef_coef[,4]<- plogis(cuts[,4]-(x[,i]))-plogis(cuts[,3]-(x[,i]))
      reef_coef[,5]<- 1-plogis(cuts[,4]-(x[,i]))
      reef_coef[,6]<- apply(reef_coef[,1:5],1,abund_tranfs)
    }
    lambda_mat[[i]]=reef_coef
    
  }  
  #  median.x=median(do.call(rbind, lapply(lambda_mat, data.frame, stringsAsFactors=FALSE))$lambda.x)
  # median.y=median(do.call(rbind, lapply(lambda_mat, data.frame, stringsAsFactors=FALSE))$lambda.y)
  
  x_mat<- data.frame(year=seq(min(ts$year),max(ts$year)),median.reef=NA,l.95.rf=NA,u.95.rf=NA)
  for(i in 1:TT){
    x_mat[i,2]=median(lambda_mat[[i]]$lambda.x)
    x_mat[i,3]=quantile(lambda_mat[[i]]$lambda.x,0.1)
    x_mat[i,4]=quantile(lambda_mat[[i]]$lambda.x,0.9)
  }
  
  
  par(xpd=T)
  plot(x_mat$median.reef~x_mat$year,type='n',ylim=c(min(x_mat[,2:4]),max((x_mat[,2:4]))),col='darkblue',bty='l',ylab=expression('Mean count per survey'),xlab='Year')
  lines(x_mat[,2]~ts$year,lty=5,lwd=2,col=cols[1])
  x<- c(x_mat$year, rev(x_mat$year))
  y1<- c(x_mat[,3], rev(x_mat[,4]))
  polygon(x, y1, col = adjustcolor(cols[1], alpha = 0.2), border=NA) # Add uncertainty polygon
  
  lines(y_mat$median.reef~y_mat$year,col=cols[2],lwd=2)
  points(y_mat$median.reef~y_mat$year,col='white',pch=21,bg=cols[2],cex=1.5)
  text(y=rep(min(x_mat),nrow(ts)),x=ts$year,ts$n.survs,cex=0.5)
  #  dev.off(paste(paste(sp,GZ,sep='_'),'.pdf',sep=''))
  #  dev.off()
}

State_space_timeseries_plot_old<- function(sp,GZ,params,TT,ts,cols,n.iter){
  #  pdf(paste('./figures/',paste(sp,GZ,'state_space',sep='_'),'.pdf',sep=''),width=8,height=6)
  
  cuts=params$c
  x=params$x #remove x0
  y=params$a_yr
  
  lambda_mat<- list()  
  for(i in 1:TT){
    reef_coef<- data.frame(p_0.y=NA,p_1.y=NA,p_2.y=NA,p_11.y=NA,p_101.y=NA,p_0.x=NA,p_1.x=NA,p_2.x=NA,p_11.x=NA,p_101.x=NA,lambda.y=NA,lambda.x=NA,n.iter=seq(1,n.iter))
    
    if(ncol(cuts)==2){
      reef_coef[,1]<- plogis(cuts[,1]-y[,i])
      reef_coef[,2]<-plogis(cuts[,2]-y[,i])-plogis(cuts[,1]-y[,i])
      reef_coef[,3]<-1-plogis(cuts[,2]-y[,i])
      reef_coef[,4]<- 0
      reef_coef[,5]<- 0
      reef_coef[,6]=plogis(cuts[,1]-(x[,i]))
      reef_coef[,7]<-plogis(cuts[,2]-(x[,i]))-plogis(cuts[,1]-(x[,i]))
      reef_coef[,8]<-1-plogis(cuts[,2]-(x[,i]))
      reef_coef[,9]<- 0
      reef_coef[,10]<- 0
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    if(ncol(cuts)==3){
      reef_coef[,1]<- plogis(cuts[,1]-y[,i])
      reef_coef[,2]<-plogis(cuts[,2]-y[,i])-plogis(cuts[,1]-y[,i])
      reef_coef[,3]<-plogis(cuts[,3]-y[,i])-plogis(cuts[,2]-y[,i])
      reef_coef[,4]<- 1-plogis(cuts[,3]-y[,i])
      reef_coef[,5]<- 0
      reef_coef[,6]=plogis(cuts[,1]-(x[,i]))
      reef_coef[,7]<-plogis(cuts[,2]-(x[,i]))-plogis(cuts[,1]-(x[,i]))
      reef_coef[,8]<-plogis(cuts[,3]-(x[,i]))-plogis(cuts[,2]-(x[,i]))
      reef_coef[,9]<- 1-plogis(cuts[,3]-(x[,i]))
      reef_coef[,10]<- 0
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    
    if(ncol(cuts)==4){
      reef_coef[,1]<- plogis(cuts[,1]-y[,i])
      reef_coef[,2]<-plogis(cuts[,2]-y[,i])-plogis(cuts[,1]-y[,i])
      reef_coef[,3]<-plogis(cuts[,3]-y[,i])-plogis(cuts[,2]-y[,i])
      reef_coef[,4]<- plogis(cuts[,4]-y[,i])-plogis(cuts[,3]-y[,i])
      reef_coef[,5]<- 1- plogis(cuts[,4]-y[,i])
      reef_coef[,6]=plogis(cuts[,1]-(x[,i]))
      reef_coef[,7]<-plogis(cuts[,2]-(x[,i]))-plogis(cuts[,1]-(x[,i]))
      reef_coef[,8]<-plogis(cuts[,3]-(x[,i]))-plogis(cuts[,2]-(x[,i]))
      reef_coef[,9]<- plogis(cuts[,4]-(x[,i]))-plogis(cuts[,3]-(x[,i]))
      reef_coef[,10]<- 1-plogis(cuts[,4]-(x[,i]))
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    lambda_mat[[i]]=reef_coef
    
  }  
  #  median.x=median(do.call(rbind, lapply(lambda_mat, data.frame, stringsAsFactors=FALSE))$lambda.x)
  # median.y=median(do.call(rbind, lapply(lambda_mat, data.frame, stringsAsFactors=FALSE))$lambda.y)
  
  x_mat<- data.frame(median.reef=NA,l.95.rf=NA,u.95.rf=NA)
  for(i in 1:TT){
    x_mat[i,1]=median(lambda_mat[[i]]$lambda.x)
    x_mat[i,2]=quantile(lambda_mat[[i]]$lambda.x,0.1)
    x_mat[i,3]=quantile(lambda_mat[[i]]$lambda.x,0.9)
  }
  
  y_mat<- data.frame(year=seq(min(ts$year),max(ts$year)),median.reef=NA,l.95.rf=NA,u.95.rf=NA)
  
  for(i in 1:TT){
    y_mat[i,2]=median(lambda_mat[[i]]$lambda.y)
    y_mat[i,3]=quantile(lambda_mat[[i]]$lambda.y,0.1)
    y_mat[i,4]=quantile(lambda_mat[[i]]$lambda.y,0.9)
  }
  
  par(xpd=T)
  plot(y_mat$median.reef~y_mat$year,type='n',ylim=c(min(x_mat),max(c(max(na.omit(y_mat[,2])),max(x_mat)))),col='darkblue',bty='l',ylab=expression('Mean count per survey'),xlab='Year',main=paste(sp,GZ,sep=', '))
  lines(x_mat[,1]~y_mat$year,lty=5,lwd=2,col=cols[1])
  x<- c(y_mat$year, rev(y_mat$year))
  y1<- c(x_mat[,2], rev(x_mat[,3]))
  polygon(x, y1, col = adjustcolor(cols[1], alpha = 0.2), border=NA) # Add uncertainty polygon
  
  lines(y_mat$median.reef~y_mat$year,col=cols[2],lwd=2)
  points(y_mat$median.reef~y_mat$year,col='white',pch=21,bg=cols[2],cex=1.5)
  text(y=rep(min(x_mat),nrow(ts)),x=ts$year,ts$n.survs,cex=0.5)
  #  dev.off(paste(paste(sp,GZ,sep='_'),'.pdf',sep=''))
  #  dev.off()
}

est_ts_plot<- function(sp,grp,params,TT,ts,n.groups,col.palette1,col.palette2,n.iter,log=0){
  lambda_mat<- list()
  x_mat<- list()
  y_mat<- list()
  
  cuts=as.data.frame(params)[grepl('cut',colnames(params))]
  x=as.data.frame(params)[grepl('x',colnames(params))]; x=x[,-1] #remove x0
  y=as.data.frame(params)[grepl('a_yr',colnames(params))]; 
  yl=list();xl=list()
  for(r in 1:n.groups){
    yl[[r]]=y[grepl(paste(',',r,sep=''),colnames(y))]
    xl[[r]]=x[grepl(paste(',',r,sep=''),colnames(x))]
  }
 
  for(r in 1:n.groups){
    ts_estimates<- list()
    for(i in 1:TT){
      reef_coef<- data.frame(p_0.y=NA,p_1.y=NA,p_2.y=NA,p_11.y=NA,p_101.y=NA,p_0.x=NA,p_1.x=NA,p_2.x=NA,p_11.x=NA,p_101.x=NA,lambda.y=NA,lambda.x=NA,n.iter=seq(1,n.iter))
      
      if(ncol(cuts)==2){
        reef_coef[,1]<- plogis(cuts[,1]-yl[[r]][,i])
        reef_coef[,2]<-plogis(cuts[,2]-yl[[r]][,i])-plogis(cuts[,1]-yl[[r]][,i])
        reef_coef[,3]<-1-plogis(cuts[,2]-yl[[r]][,i])
        reef_coef[,4]<- 0
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(cuts[,1]-(xl[[r]][,i]))
        reef_coef[,7]<-plogis(cuts[,2]-(xl[[r]][,i]))-plogis(cuts[,1]-(xl[[r]][,i]))
        reef_coef[,8]<-1-plogis(cuts[,2]-(xl[[r]][,i]))
        reef_coef[,9]<- 0
        reef_coef[,10]<- 0
        reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
      }
      if(ncol(cuts)==3){
        reef_coef[,1]<- plogis(cuts[,1]-yl[[r]][,i])
        reef_coef[,2]<-plogis(cuts[,2]-yl[[r]][,i])-plogis(cuts[,1]-yl[[r]][,i])
        reef_coef[,3]<-plogis(cuts[,3]-yl[[r]][,i])-plogis(cuts[,2]-yl[[r]][,i])
        reef_coef[,4]<- 1-plogis(cuts[,3]-yl[[r]][,i])
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(cuts[,1]-(xl[[r]][,i]))
        reef_coef[,7]<-plogis(cuts[,2]-(xl[[r]][,i]))-plogis(cuts[,1]-(xl[[r]][,i]))
        reef_coef[,8]<-plogis(cuts[,3]-(xl[[r]][,i]))-plogis(cuts[,2]-(xl[[r]][,i]))
        reef_coef[,9]<- 1-plogis(cuts[,3]-(xl[[r]][,i]))
        reef_coef[,10]<- 0
        reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
      }
      
      if(ncol(cuts)==4){
        reef_coef[,1]<- plogis(cuts[,1]-yl[[r]][,i])
        reef_coef[,2]<-plogis(cuts[,2]-yl[[r]][,i])-plogis(cuts[,1]-yl[[r]][,i])
        reef_coef[,3]<-plogis(cuts[,3]-yl[[r]][,i])-plogis(cuts[,2]-yl[[r]][,i])
        reef_coef[,4]<- plogis(cuts[,4]-yl[[r]][,i])-plogis(cuts[,3]-yl[[r]][,i])
        reef_coef[,5]<- 1- plogis(cuts[,4]-yl[[r]][,i])
        reef_coef[,6]=plogis(cuts[,1]-(xl[[r]][,i]))
        reef_coef[,7]<-plogis(cuts[,2]-(xl[[r]][,i]))-plogis(cuts[,1]-(xl[[r]][,i]))
        reef_coef[,8]<-plogis(cuts[,3]-(xl[[r]][,i]))-plogis(cuts[,2]-(xl[[r]][,i]))
        reef_coef[,9]<- plogis(cuts[,4]-(xl[[r]][,i]))-plogis(cuts[,3]-(xl[[r]][,i]))
        reef_coef[,10]<- 1-plogis(cuts[,4]-(xl[[r]][,i]))
        reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
      }
      ts_estimates[[i]]=reef_coef
      }
    lambda_mat[[r]]=ts_estimates
  } 
  
  for(r in 1:n.groups){  
    x_mat[[r]]<- data.frame(median.reef=NA,l.80.rf=NA,u.80.rf=NA,l.90.rf=NA,u.90.rf=NA)
    for(i in 1:TT){
      x_mat[[r]][i,1]=median(lambda_mat[[r]][[i]]$lambda.x)
      x_mat[[r]][i,2]=quantile(lambda_mat[[r]][[i]]$lambda.x,0.1)
      x_mat[[r]][i,3]=quantile(lambda_mat[[r]][[i]]$lambda.x,0.9)
      x_mat[[r]][i,4]=quantile(lambda_mat[[r]][[i]]$lambda.x,0.05)
      x_mat[[r]][i,5]=quantile(lambda_mat[[r]][[i]]$lambda.x,0.95)
    }
    
    y_mat[[r]]<- data.frame(year=seq(min(ts$year),max(ts$year)),median.reef=NA,l.80.rf=NA,u.80.rf=NA,l.95.rf=NA,u.95.rf=NA)
    for(i in 1:TT){
      y_mat[[r]][i,2]=median(lambda_mat[[r]][[i]]$lambda.y)
      y_mat[[r]][i,3]=quantile(lambda_mat[[r]][[i]]$lambda.y,0.1)
      y_mat[[r]][i,4]=quantile(lambda_mat[[r]][[i]]$lambda.y,0.9)
      y_mat[[r]][i,5]=quantile(lambda_mat[[r]][[i]]$lambda.y,0.05)
      y_mat[[r]][i,6]=quantile(lambda_mat[[r]][[i]]$lambda.y,0.95)
    }
  }
  x_mat_full<- do.call(rbind, lapply(x_mat, data.frame, stringsAsFactors=FALSE))
  y_mat_full<- do.call(rbind, lapply(y_mat, data.frame, stringsAsFactors=FALSE))
  par(xpd=T,mar=c(4.1,4.1,3.1,3.1))
  if(log==1){
    plot(y_mat[[1]]$median.reef~y_mat[[1]]$year,type='n',ylim=c(min(x_mat_full),max(c(max(y_mat_full[,2]),max(x_mat_full)))),col='darkblue',bty='l',ylab=expression('Mean count per survey'),xlab='Year',log='y')
  }
  if(log==0){
    plot(y_mat[[1]]$median.reef~y_mat[[1]]$year,type='n',ylim=c(min(x_mat_full),max(c(max(y_mat_full[,2]),max(x_mat_full)))),col='darkblue',bty='l',ylab=expression('Mean count per survey'),xlab='Year')
  }
  for(r in 1:n.groups){
    lines(x_mat[[r]][,1]~y_mat[[r]]$year,lty=5,lwd=2,col=col.palette1[r])
    x<- c(y_mat[[r]]$year, rev(y_mat[[r]]$year))
    y.90<- c(x_mat[[r]][,4], rev(x_mat[[r]][,5]))
#    y.80<- c(x_mat[[r]][,2], rev(x_mat[[r]][,3]))
    polygon(x, y.90, col = adjustcolor(col.palette1[r], alpha = 0.2), border=NA) # Add uncertainty polygon
#    polygon(x, y.80, col = adjustcolor(col.palette1[r], alpha = 0.2), border=NA) # Add uncertainty polygon
    
    lines(y_mat[[r]]$median.reef~y_mat[[r]]$year,col=col.palette2[r],lwd=2)
    points(y_mat[[r]]$median.reef~y_mat[[r]]$year,col='white',pch=21,bg=col.palette2[r],cex=1.5)
  }
  text(y=rep(par("usr")[3]+(max(x_mat_full)-min(x_mat_full))*0.015,nrow(ts)),x=ts$year,ts$n.survs,cex=0.5)
  text(y=par("usr")[3]+(max(x_mat_full)-min(x_mat_full))*0.045,x=par("usr")[2]+0.2,'Surveys',cex=0.6,pos=1,xpd=T,font=1)
  for(r in 1:n.groups){
    lines(y=rep(par("usr")[4]-(max(x_mat_full)-min(x_mat_full))*0.03*r,2),x=c(par("usr")[2]-0.5,par("usr")[2]+0.5),col=col.palette2[r],xpd=T,lwd=2)
  }
  for(r in 1:n.groups){
    text(y=par("usr")[4]-(max(x_mat_full)-min(x_mat_full))*0.03*r,x=par("usr")[2]-2.5,grp[r],col=col.palette2[r],xpd=T,lwd=2,cex=0.6)
  }
}

est_ts_plot_scaled<- function(sp,grp,params,TT,ts,n.groups,col.palette1,col.palette2,log=0,n.iter){
  lambda_mat<- list()
  x_mat<- list()
  y_mat<- list()
  
  cuts=as.data.frame(params)[grepl('cut',colnames(params))]
  x=as.data.frame(params)[grepl('x',colnames(params))]; x=x[,-1] #remove x0
  y=as.data.frame(params)[grepl('a_yr',colnames(params))]; 
  yl=list();xl=list()
  
  median.x= numeric(n.groups)
  median.y= numeric(n.groups)
  for(r in 1:n.groups){
    yl[[r]]=y[grepl(paste(',',r,sep=''),colnames(y))]
    xl[[r]]=x[grepl(paste(',',r,sep=''),colnames(x))]
  }
  
  for(r in 1:n.groups){
    ts_estimates<- list()
    for(i in 1:TT){
      reef_coef<- data.frame(p_0.y=NA,p_1.y=NA,p_2.y=NA,p_11.y=NA,p_101.y=NA,p_0.x=NA,p_1.x=NA,p_2.x=NA,p_11.x=NA,p_101.x=NA,lambda.y=NA,lambda.x=NA,n.iter=seq(1,n.iter))
      
      if(ncol(cuts)==2){
        reef_coef[,1]<- plogis(cuts[,1]-yl[[r]][,i])
        reef_coef[,2]<-plogis(cuts[,2]-yl[[r]][,i])-plogis(cuts[,1]-yl[[r]][,i])
        reef_coef[,3]<-1-plogis(cuts[,2]-yl[[r]][,i])
        reef_coef[,4]<- 0
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(cuts[,1]-(xl[[r]][,i]))
        reef_coef[,7]<-plogis(cuts[,2]-(xl[[r]][,i]))-plogis(cuts[,1]-(xl[[r]][,i]))
        reef_coef[,8]<-1-plogis(cuts[,2]-(xl[[r]][,i]))
        reef_coef[,9]<- 0
        reef_coef[,10]<- 0
        reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
      }
      if(ncol(cuts)==3){
        reef_coef[,1]<- plogis(cuts[,1]-yl[[r]][,i])
        reef_coef[,2]<-plogis(cuts[,2]-yl[[r]][,i])-plogis(cuts[,1]-yl[[r]][,i])
        reef_coef[,3]<-plogis(cuts[,3]-yl[[r]][,i])-plogis(cuts[,2]-yl[[r]][,i])
        reef_coef[,4]<- 1-plogis(cuts[,3]-yl[[r]][,i])
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(cuts[,1]-(xl[[r]][,i]))
        reef_coef[,7]<-plogis(cuts[,2]-(xl[[r]][,i]))-plogis(cuts[,1]-(xl[[r]][,i]))
        reef_coef[,8]<-plogis(cuts[,3]-(xl[[r]][,i]))-plogis(cuts[,2]-(xl[[r]][,i]))
        reef_coef[,9]<- 1-plogis(cuts[,3]-(xl[[r]][,i]))
        reef_coef[,10]<- 0
        reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
      }
      
      if(ncol(cuts)==4){
        reef_coef[,1]<- plogis(cuts[,1]-yl[[r]][,i])
        reef_coef[,2]<-plogis(cuts[,2]-yl[[r]][,i])-plogis(cuts[,1]-yl[[r]][,i])
        reef_coef[,3]<-plogis(cuts[,3]-yl[[r]][,i])-plogis(cuts[,2]-yl[[r]][,i])
        reef_coef[,4]<- plogis(cuts[,4]-yl[[r]][,i])-plogis(cuts[,3]-yl[[r]][,i])
        reef_coef[,5]<- 1- plogis(cuts[,4]-yl[[r]][,i])
        reef_coef[,6]=plogis(cuts[,1]-(xl[[r]][,i]))
        reef_coef[,7]<-plogis(cuts[,2]-(xl[[r]][,i]))-plogis(cuts[,1]-(xl[[r]][,i]))
        reef_coef[,8]<-plogis(cuts[,3]-(xl[[r]][,i]))-plogis(cuts[,2]-(xl[[r]][,i]))
        reef_coef[,9]<- plogis(cuts[,4]-(xl[[r]][,i]))-plogis(cuts[,3]-(xl[[r]][,i]))
        reef_coef[,10]<- 1-plogis(cuts[,4]-(xl[[r]][,i]))
        reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
      }
      ts_estimates[[i]]=reef_coef
    }
    lambda_mat[[r]]=ts_estimates
    median.x[r]=median(do.call(rbind, lapply(lambda_mat[[r]], data.frame, stringsAsFactors=FALSE))$lambda.x)
    median.y[r]=median(do.call(rbind, lapply(lambda_mat[[r]], data.frame, stringsAsFactors=FALSE))$lambda.y)
  } 
  
  for(r in 1:n.groups){  
    x_mat[[r]]<- data.frame(median.reef=NA,l.80.rf=NA,u.80.rf=NA,l.90.rf=NA,u.90.rf=NA)
    for(i in 1:TT){
      x_mat[[r]][i,1]=median(lambda_mat[[r]][[i]]$lambda.x/median.x[r])
      x_mat[[r]][i,2]=quantile(lambda_mat[[r]][[i]]$lambda.x/median.x[r],0.1)
      x_mat[[r]][i,3]=quantile(lambda_mat[[r]][[i]]$lambda.x/median.x[r],0.9)
      x_mat[[r]][i,4]=quantile(lambda_mat[[r]][[i]]$lambda.x/median.x[r],0.05)
      x_mat[[r]][i,5]=quantile(lambda_mat[[r]][[i]]$lambda.x/median.x[r],0.95)
    }
    
    y_mat[[r]]<- data.frame(year=seq(min(ts$year),max(ts$year)),median.reef=NA,l.80.rf=NA,u.80.rf=NA,l.95.rf=NA,u.95.rf=NA)
    for(i in 1:TT){
      y_mat[[r]][i,2]=median(lambda_mat[[r]][[i]]$lambda.y/median.y[r])
      y_mat[[r]][i,3]=quantile(lambda_mat[[r]][[i]]$lambda.y/median.y[r],0.1)
      y_mat[[r]][i,4]=quantile(lambda_mat[[r]][[i]]$lambda.y/median.y[r],0.9)
      y_mat[[r]][i,5]=quantile(lambda_mat[[r]][[i]]$lambda.y/median.y[r],0.05)
      y_mat[[r]][i,6]=quantile(lambda_mat[[r]][[i]]$lambda.y/median.y[r],0.95)
    }
  }
  x_mat_full<- do.call(rbind, lapply(x_mat, data.frame, stringsAsFactors=FALSE))
  y_mat_full<- do.call(rbind, lapply(y_mat, data.frame, stringsAsFactors=FALSE))
  par(xpd=T,mar=c(4.1,4.1,3.1,3.1))
  if(log==1){
    plot(y_mat[[1]]$median.reef~y_mat[[1]]$year,type='n',ylim=c(min(x_mat_full),max(c(max(y_mat_full[,2]),max(x_mat_full)))),col='darkblue',bty='l',ylab=expression('Scaled relative abundance'),xlab='Year',main=sp,log='y')
  }
  if(log==0){
    plot(y_mat[[1]]$median.reef~y_mat[[1]]$year,type='n',ylim=c(min(x_mat_full),max(c(max(y_mat_full[,2]),max(x_mat_full)))),col='darkblue',bty='l',ylab=expression('Scaled relative abundance'),xlab='Year',main=sp)
  }
  for(r in 1:n.groups){
    lines(x_mat[[r]][,1]~y_mat[[r]]$year,lty=5,lwd=2,col=col.palette1[r])
    x<- c(y_mat[[r]]$year, rev(y_mat[[r]]$year))
    y.90<- c(x_mat[[r]][,4], rev(x_mat[[r]][,5]))
    #    y.80<- c(x_mat[[r]][,2], rev(x_mat[[r]][,3]))
    polygon(x, y.90, col = adjustcolor(col.palette1[r], alpha = 0.2), border=NA) # Add uncertainty polygon
    #    polygon(x, y.80, col = adjustcolor(col.palette1[r], alpha = 0.2), border=NA) # Add uncertainty polygon
    
    lines(y_mat[[r]]$median.reef~y_mat[[r]]$year,col=col.palette2[r],lwd=2)
    points(y_mat[[r]]$median.reef~y_mat[[r]]$year,col='white',pch=21,bg=col.palette2[r],cex=1.5)
  }
  text(y=rep(par("usr")[3]+(max(x_mat_full)-min(x_mat_full))*0.015,nrow(ts)),x=ts$year,ts$n.survs,cex=0.5)
  text(y=par("usr")[3]+(max(x_mat_full)-min(x_mat_full))*0.045,x=par("usr")[2]+0.2,'Surveys',cex=0.6,pos=1,xpd=T,font=1)
  for(r in 1:n.groups){
    lines(y=rep(par("usr")[4]-(max(x_mat_full)-min(x_mat_full))*0.03*r,2),x=c(par("usr")[2]-0.5,par("usr")[2]+0.5),col=col.palette2[r],xpd=T,lwd=2)
  }
  for(r in 1:n.groups){
    text(y=par("usr")[4]-(max(x_mat_full)-min(x_mat_full))*0.03*r,x=par("usr")[2]-2.5,grp[r],col=col.palette2[r],xpd=T,lwd=2,cex=0.6)
  }
  
}

Multi_state_space_timeseries_plot<- function(sp,grp,params,TT,ts,n.groups,col.palette1,col.palette2,log=0,n.iter){
  lambda_mat<- list()
  x_mat<- list()
  y_mat<- list()

  cuts=as.data.frame(params)[grepl('cut',colnames(params))];
  x=as.data.frame(params)[grepl('x',colnames(params))]
  x=x[,-c(1:n.groups)] #remove x0s
  y=as.data.frame(params)[grepl('a_yr',colnames(params))]; 
  x_r=list();y_r=list()
  for(r in 1:n.groups){
    x_r[[r]]=x[grepl(paste(',',r,']',sep=''),colnames(x))]
    y_r[[r]]=y[grepl(paste(',',r,']',sep=''),colnames(y))]
  }
  
  
  for(r in 1:n.groups){
    ts_estimates<- list()
    for(i in 1:TT){
      reef_coef<- data.frame(p_0.y=NA,p_1.y=NA,p_2.y=NA,p_11.y=NA,p_101.y=NA,p_0.x=NA,p_1.x=NA,p_2.x=NA,p_11.x=NA,p_101.x=NA,lambda.y=NA,lambda.x=NA,n.iter=seq(1,n.iter))
      
      if(ncol(cuts)==2){
        reef_coef[,1]<- plogis(cuts[,1]-y_r[[r]][,i])
        reef_coef[,2]<-plogis(cuts[,2]-y_r[[r]][,i])-plogis(cuts[,1]-y_r[[r]][,i])
        reef_coef[,3]<-1-plogis(cuts[,2]-y_r[[r]][,i])
        reef_coef[,4]<- 0
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(cuts[,1]-(x_r[[r]][,i]))
        reef_coef[,7]<-plogis(cuts[,2]-(x_r[[r]][,i]))-plogis(cuts[,1]-(x_r[[r]][,i]))
        reef_coef[,8]<-1-plogis(cuts[,2]-(x_r[[r]][,i]))
        reef_coef[,9]<- 0
        reef_coef[,10]<- 0
        reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
      }
      if(ncol(cuts)==3){
        reef_coef[,1]<- plogis(cuts[,1]-y_r[[r]][,i])
        reef_coef[,2]<-plogis(cuts[,2]-y_r[[r]][,i])-plogis(cuts[,1]-y_r[[r]][,i])
        reef_coef[,3]<-plogis(cuts[,3]-y_r[[r]][,i])-plogis(cuts[,2]-y_r[[r]][,i])
        reef_coef[,4]<- 1-plogis(cuts[,3]-y_r[[r]][,i])
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(cuts[,1]-(x_r[[r]][,i]))
        reef_coef[,7]<-plogis(cuts[,2]-(x_r[[r]][,i]))-plogis(cuts[,1]-(x_r[[r]][,i]))
        reef_coef[,8]<-plogis(cuts[,3]-(x_r[[r]][,i]))-plogis(cuts[,2]-(x_r[[r]][,i]))
        reef_coef[,9]<- 1-plogis(cuts[,3]-(x_r[[r]][,i]))
        reef_coef[,10]<- 0
        reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
      }
      
      if(ncol(cuts)==4){
        reef_coef[,1]<- plogis(cuts[,1]-y_r[[r]][,i])
        reef_coef[,2]<-plogis(cuts[,2]-y_r[[r]][,i])-plogis(cuts[,1]-y_r[[r]][,i])
        reef_coef[,3]<-plogis(cuts[,3]-y_r[[r]][,i])-plogis(cuts[,2]-y_r[[r]][,i])
        reef_coef[,4]<- plogis(cuts[,4]-y_r[[r]][,i])-plogis(cuts[,3]-y_r[[r]][,i])
        reef_coef[,5]<- 1- plogis(cuts[,4]-y_r[[r]][,i])
        reef_coef[,6]=plogis(cuts[,1]-(x_r[[r]][,i]))
        reef_coef[,7]<-plogis(cuts[,2]-(x_r[[r]][,i]))-plogis(cuts[,1]-(x_r[[r]][,i]))
        reef_coef[,8]<-plogis(cuts[,3]-(x_r[[r]][,i]))-plogis(cuts[,2]-(x_r[[r]][,i]))
        reef_coef[,9]<- plogis(cuts[,4]-(x_r[[r]][,i]))-plogis(cuts[,3]-(x_r[[r]][,i]))
        reef_coef[,10]<- 1-plogis(cuts[,4]-(x_r[[r]][,i]))
        reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
      }
      ts_estimates[[i]]=reef_coef
    }
    lambda_mat[[r]]=ts_estimates
  } 
  
  
  for(r in 1:n.groups){  
    x_mat[[r]]<- data.frame(median.reef=NA,l.80.rf=NA,u.80.rf=NA,l.90.rf=NA,u.90.rf=NA)
    for(i in 1:TT){
      x_mat[[r]][i,1]=median(lambda_mat[[r]][[i]]$lambda.x)
      x_mat[[r]][i,2]=quantile(lambda_mat[[r]][[i]]$lambda.x,0.1)
      x_mat[[r]][i,3]=quantile(lambda_mat[[r]][[i]]$lambda.x,0.9)
      x_mat[[r]][i,4]=quantile(lambda_mat[[r]][[i]]$lambda.x,0.05)
      x_mat[[r]][i,5]=quantile(lambda_mat[[r]][[i]]$lambda.x,0.95)
    }
    
    y_mat[[r]]<- data.frame(year=seq(min(ts$year),max(ts$year)),median.reef=NA,l.80.rf=NA,u.80.rf=NA,l.95.rf=NA,u.95.rf=NA)
    for(i in 1:TT){
      y_mat[[r]][i,2]=median(lambda_mat[[r]][[i]]$lambda.y)
      y_mat[[r]][i,3]=quantile(lambda_mat[[r]][[i]]$lambda.y,0.1)
      y_mat[[r]][i,4]=quantile(lambda_mat[[r]][[i]]$lambda.y,0.9)
      y_mat[[r]][i,5]=quantile(lambda_mat[[r]][[i]]$lambda.y,0.05)
      y_mat[[r]][i,6]=quantile(lambda_mat[[r]][[i]]$lambda.y,0.95)
    }
  }
  x_mat_full<- do.call(rbind, lapply(x_mat, data.frame, stringsAsFactors=FALSE))
  y_mat_full<- do.call(rbind, lapply(y_mat, data.frame, stringsAsFactors=FALSE))
  par(xpd=T,mar=c(4.1,4.1,3.1,3.1))
  if(log==1){
    plot(y_mat[[1]]$median.reef~y_mat[[1]]$year,type='n',ylim=c(min(x_mat_full),max(c(max(y_mat_full[,2]),max(x_mat_full)))),col='darkblue',bty='l',ylab=expression('Mean count per survey'),xlab='Year',main=sp,log='y')
  }
  if(log==0){
    plot(y_mat[[1]]$median.reef~y_mat[[1]]$year,type='n',ylim=c(min(x_mat_full),max(c(max(y_mat_full[,2]),max(x_mat_full)))),col='darkblue',bty='l',ylab=expression('Mean count per survey'),xlab='Year',main=sp)
  }
  for(r in 1:n.groups){
    lines(x_mat[[r]][,1]~y_mat[[r]]$year,lty=5,lwd=2,col=col.palette1[r])
    x<- c(y_mat[[r]]$year, rev(y_mat[[r]]$year))
    y.90<- c(x_mat[[r]][,4], rev(x_mat[[r]][,5]))
    #    y.80<- c(x_mat[[r]][,2], rev(x_mat[[r]][,3]))
    polygon(x, y.90, col = adjustcolor(col.palette1[r], alpha = 0.2), border=NA) # Add uncertainty polygon
    #    polygon(x, y.80, col = adjustcolor(col.palette1[r], alpha = 0.2), border=NA) # Add uncertainty polygon
    
    lines(y_mat[[r]]$median.reef~y_mat[[r]]$year,col=col.palette2[r],lwd=2)
    points(y_mat[[r]]$median.reef~y_mat[[r]]$year,col='white',pch=21,bg=col.palette2[r],cex=1.5)
  }
  text(y=rep(par("usr")[3]+(max(x_mat_full)-min(x_mat_full))*0.015,nrow(ts)),x=ts$year,ts$n.survs,cex=0.5)
  text(y=par("usr")[3]+(max(x_mat_full)-min(x_mat_full))*0.045,x=par("usr")[2]+0.2,'Surveys',cex=0.6,pos=1,xpd=T,font=1)
  for(r in 1:n.groups){
    lines(y=rep(par("usr")[4]-(max(x_mat_full)-min(x_mat_full))*0.03*r,2),x=c(par("usr")[2]-0.5,par("usr")[2]+0.5),col=col.palette2[r],xpd=T,lwd=2)
  }
  for(r in 1:n.groups){
    text(y=par("usr")[4]-(max(x_mat_full)-min(x_mat_full))*0.03*r,x=par("usr")[2]-2.5,grp[r],col=col.palette2[r],xpd=T,lwd=2,cex=0.6)
  }
}

Multi_state_space_timeseries_plot_scaled<- function(sp,grp,params1,TT,ts,n.groups,col.palette1,col.palette2,log){
  lambda_mat<- list()
  x_mat<- list()
  y_mat<- list()
  median.x= numeric(n.groups)
  median.y= numeric(n.groups)
  for(r in 1:n.groups){
    ts_estimates<- list()
    for(i in 1:TT){
      reef_coef<- data.frame(p_0.y=NA,p_1.y=NA,p_2.y=NA,p_11.y=NA,p_101.y=NA,p_0.x=NA,p_1.x=NA,p_2.x=NA,p_11.x=NA,p_101.x=NA,lambda.y=NA,lambda.x=NA,iter=seq(1,nrow(params1$c)))
      
      if(ncol(params1$c)==2){
        reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,,r][,i])-plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,3]<-1-plogis(params1$c[,2]-params1$a_yr[,,r][,i])
        reef_coef[,4]<- 0
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,,r][,i])-plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,8]<-1-plogis(params1$c[,2]-params1$x[,,r][,i])
        reef_coef[,9]<- 0
        reef_coef[,10]<- 0
        
        reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
      }
      
      if(ncol(params1$c)==3){
        reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,,r][,i])-plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,3]<-plogis(params1$c[,3]-params1$a_yr[,,r][,i])-plogis(params1$c[,2]-params1$a_yr[,,r][,i])
        reef_coef[,4]<- 1-plogis(params1$c[,3]-params1$a_yr[,,r][,i])
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,,r][,i])-plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,,r][,i])-plogis(params1$c[,2]-params1$x[,,r][,i])
        reef_coef[,9]<- 1-plogis(params1$c[,3]-params1$x[,,r][,i])
        reef_coef[,10]<- 0
        
        reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
      }
      
      if(ncol(params1$c)==4){
        reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,,r][,i])-plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,3]<-plogis(params1$c[,3]-params1$a_yr[,,r][,i])-plogis(params1$c[,2]-params1$a_yr[,,r][,i])
        reef_coef[,4]<- plogis(params1$c[,4]-params1$a_yr[,,r][,i])-plogis(params1$c[,3]-params1$a_yr[,,r][,i])
        reef_coef[,5]<- 1- plogis(params1$c[,4]-params1$a_yr[,,r][,i])
        reef_coef[,6]=plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,,r][,i])-plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,,r][,i])-plogis(params1$c[,2]-params1$x[,,r][,i])
        reef_coef[,9]<- plogis(params1$c[,4]-params1$x[,,r][,i])-plogis(params1$c[,3]-params1$x[,,r][,i])
        reef_coef[,10]<- 1-plogis(params1$c[,4]-params1$x[,,r][,i])
        
        reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
      }
      ts_estimates[[i]]=reef_coef
    }
    lambda_mat[[r]]=ts_estimates
    median.x[r]=median(do.call(rbind, lapply(lambda_mat[[r]], data.frame, stringsAsFactors=FALSE))$lambda.x)
    median.y[r]=median(do.call(rbind, lapply(lambda_mat[[r]], data.frame, stringsAsFactors=FALSE))$lambda.y)
  } 
  
  
  for(r in 1:n.groups){  
    x_mat[[r]]<- data.frame(median.reef=NA,l.80.rf=NA,u.80.rf=NA,l.90.rf=NA,u.90.rf=NA)
    for(i in 1:TT){
      x_mat[[r]][i,1]=median(lambda_mat[[r]][[i]]$lambda.x/median.x[r])
      x_mat[[r]][i,2]=quantile(lambda_mat[[r]][[i]]$lambda.x/median.x[r],0.1)
      x_mat[[r]][i,3]=quantile(lambda_mat[[r]][[i]]$lambda.x/median.x[r],0.9)
      x_mat[[r]][i,4]=quantile(lambda_mat[[r]][[i]]$lambda.x/median.x[r],0.05)
      x_mat[[r]][i,5]=quantile(lambda_mat[[r]][[i]]$lambda.x/median.x[r],0.95)
    }
    
    y_mat[[r]]<- data.frame(year=seq(min(ts$year),max(ts$year)),median.reef=NA,l.80.rf=NA,u.80.rf=NA,l.95.rf=NA,u.95.rf=NA)
    for(i in 1:TT){
      y_mat[[r]][i,2]=median(lambda_mat[[r]][[i]]$lambda.y/median.y[r])
      y_mat[[r]][i,3]=quantile(lambda_mat[[r]][[i]]$lambda.y/median.y[r],0.1)
      y_mat[[r]][i,4]=quantile(lambda_mat[[r]][[i]]$lambda.y/median.y[r],0.9)
      y_mat[[r]][i,5]=quantile(lambda_mat[[r]][[i]]$lambda.y/median.y[r],0.05)
      y_mat[[r]][i,6]=quantile(lambda_mat[[r]][[i]]$lambda.y/median.y[r],0.95)
    }
  }
  x_mat_full<- do.call(rbind, lapply(x_mat, data.frame, stringsAsFactors=FALSE))
  y_mat_full<- do.call(rbind, lapply(y_mat, data.frame, stringsAsFactors=FALSE))
  par(xpd=T,mar=c(4.1,4.1,3.1,3.1))
  if(log==1){
    plot(y_mat[[1]]$median.reef~y_mat[[1]]$year,type='n',ylim=c(min(x_mat_full),max(c(max(y_mat_full[,2]),max(x_mat_full)))),col='darkblue',bty='l',ylab=expression('Scaled relative abundance'),xlab='Year',main=sp,log='y')
  }
  if(log==0){
    plot(y_mat[[1]]$median.reef~y_mat[[1]]$year,type='n',ylim=c(min(x_mat_full),max(c(max(y_mat_full[,2]),max(x_mat_full)))),col='darkblue',bty='l',ylab=expression('Scaled relative abundance'),xlab='Year',main=sp)
  }
  for(r in 1:n.groups){
    lines(x_mat[[r]][,1]~y_mat[[r]]$year,lty=5,lwd=2,col=col.palette1[r])
    x<- c(y_mat[[r]]$year, rev(y_mat[[r]]$year))
    y.90<- c(x_mat[[r]][,4], rev(x_mat[[r]][,5]))
    #    y.80<- c(x_mat[[r]][,2], rev(x_mat[[r]][,3]))
    polygon(x, y.90, col = adjustcolor(col.palette1[r], alpha = 0.2), border=NA) # Add uncertainty polygon
    #    polygon(x, y.80, col = adjustcolor(col.palette1[r], alpha = 0.2), border=NA) # Add uncertainty polygon
    
    lines(y_mat[[r]]$median.reef~y_mat[[r]]$year,col=col.palette2[r],lwd=2)
    points(y_mat[[r]]$median.reef~y_mat[[r]]$year,col='white',pch=21,bg=col.palette2[r],cex=1.5)
  }
  text(y=rep(par("usr")[3]+(max(x_mat_full)-min(x_mat_full))*0.015,nrow(ts)),x=ts$year,ts$n.survs,cex=0.5)
  text(y=par("usr")[3]+(max(x_mat_full)-min(x_mat_full))*0.045,x=par("usr")[2]+0.2,'Surveys',cex=0.6,pos=1,xpd=T,font=1)
  for(r in 1:n.groups){
    lines(y=rep(par("usr")[4]-(max(x_mat_full)-min(x_mat_full))*0.03*r,2),x=c(par("usr")[2]-0.5,par("usr")[2]+0.5),col=col.palette2[r],xpd=T,lwd=2)
  }
  for(r in 1:n.groups){
    text(y=par("usr")[4]-(max(x_mat_full)-min(x_mat_full))*0.03*r,x=par("usr")[2]-2.5,grp[r],col=col.palette2[r],xpd=T,lwd=2,cex=0.6)
  }

}

Multi_state_space_timeseries_plot_pdf<- function(sp,grp,params1,TT,ts,n.groups,col.palette1,col.palette2){
  pdf(file.path(path,paste(sprintf("%02d",i),'_',sp,'_',GZ,'.pdf',sep='')),width=8,height=6)
  lambda_mat<- list()
  x_mat<- list()
  y_mat<- list()
  for(r in 1:n.groups){
    ts_estimates<- list()
    for(i in 1:TT){
      reef_coef<- data.frame(p_0.y=NA,p_1.y=NA,p_2.y=NA,p_11.y=NA,p_101.y=NA,p_0.x=NA,p_1.x=NA,p_2.x=NA,p_11.x=NA,p_101.x=NA,lambda.y=NA,lambda.x=NA,iter=seq(1,nrow(params$c)))
      
      if(ncol(params1$c)==2){
        reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,,r][,i])-plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,3]<-1-plogis(params1$c[,2]-params1$a_yr[,,r][,i])
        reef_coef[,4]<- 0
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,,r][,i])-plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,8]<-1-plogis(params1$c[,2]-params1$x[,,r][,i])
        reef_coef[,9]<- 0
        reef_coef[,10]<- 0
        
        reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
      }
      
      if(ncol(params1$c)==3){
        reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,,r][,i])-plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,3]<-plogis(params1$c[,3]-params1$a_yr[,,r][,i])-plogis(params1$c[,2]-params1$a_yr[,,r][,i])
        reef_coef[,4]<- 1-plogis(params1$c[,3]-params1$a_yr[,,r][,i])
        reef_coef[,5]<- 0
        reef_coef[,6]=plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,,r][,i])-plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,,r][,i])-plogis(params1$c[,2]-params1$x[,,r][,i])
        reef_coef[,9]<- 1-plogis(params1$c[,3]-params1$x[,,r][,i])
        reef_coef[,10]<- 0
        
        reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
      }
      
      if(ncol(params1$c)==4){
        reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,,r][,i])-plogis(params1$c[,1]-params1$a_yr[,,r][,i])
        reef_coef[,3]<-plogis(params1$c[,3]-params1$a_yr[,,r][,i])-plogis(params1$c[,2]-params1$a_yr[,,r][,i])
        reef_coef[,4]<- plogis(params1$c[,4]-params1$a_yr[,,r][,i])-plogis(params1$c[,3]-params1$a_yr[,,r][,i])
        reef_coef[,5]<- 1- plogis(params1$c[,4]-params1$a_yr[,,r][,i])
        reef_coef[,6]=plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,,r][,i])-plogis(params1$c[,1]-params1$x[,,r][,i])
        reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,,r][,i])-plogis(params1$c[,2]-params1$x[,,r][,i])
        reef_coef[,9]<- plogis(params1$c[,4]-params1$x[,,r][,i])-plogis(params1$c[,3]-params1$x[,,r][,i])
        reef_coef[,10]<- 1-plogis(params1$c[,4]-params1$x[,,r][,i])
        
        reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
        reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
      }
      ts_estimates[[i]]=reef_coef
    }
    lambda_mat[[r]]=ts_estimates
  } 
  
  
  for(r in 1:n.groups){  
    x_mat[[r]]<- data.frame(median.reef=NA,l.80.rf=NA,u.80.rf=NA,l.90.rf=NA,u.90.rf=NA)
    for(i in 1:TT){
      x_mat[[r]][i,1]=median(lambda_mat[[r]][[i]]$lambda.x)
      x_mat[[r]][i,2]=quantile(lambda_mat[[r]][[i]]$lambda.x,0.1)
      x_mat[[r]][i,3]=quantile(lambda_mat[[r]][[i]]$lambda.x,0.9)
      x_mat[[r]][i,4]=quantile(lambda_mat[[r]][[i]]$lambda.x,0.05)
      x_mat[[r]][i,5]=quantile(lambda_mat[[r]][[i]]$lambda.x,0.95)
    }
    
    y_mat[[r]]<- data.frame(year=seq(min(ts$year),max(ts$year)),median.reef=NA,l.80.rf=NA,u.80.rf=NA,l.95.rf=NA,u.95.rf=NA)
    for(i in 1:TT){
      y_mat[[r]][i,2]=median(lambda_mat[[r]][[i]]$lambda.y)
      y_mat[[r]][i,3]=quantile(lambda_mat[[r]][[i]]$lambda.y,0.1)
      y_mat[[r]][i,4]=quantile(lambda_mat[[r]][[i]]$lambda.y,0.9)
      y_mat[[r]][i,5]=quantile(lambda_mat[[r]][[i]]$lambda.y,0.05)
      y_mat[[r]][i,6]=quantile(lambda_mat[[r]][[i]]$lambda.y,0.95)
    }
  }
  x_mat_full<- do.call(rbind, lapply(x_mat, data.frame, stringsAsFactors=FALSE))
  y_mat_full<- do.call(rbind, lapply(y_mat, data.frame, stringsAsFactors=FALSE))
  par(xpd=T,mar=c(4.1,4.1,3.1,3.1))
  plot(y_mat[[1]]$median.reef~y_mat[[1]]$year,type='n',ylim=c(min(x_mat_full),max(c(max(y_mat_full[,2]),max(x_mat_full)))),col='darkblue',bty='l',ylab=expression('Mean count per survey'),xlab='Year',main=sp)
  for(r in 1:n.groups){
    lines(x_mat[[r]][,1]~y_mat[[r]]$year,lty=5,lwd=2,col=col.palette1[r])
    x<- c(y_mat[[r]]$year, rev(y_mat[[r]]$year))
    y.90<- c(x_mat[[r]][,4], rev(x_mat[[r]][,5]))
    #    y.80<- c(x_mat[[r]][,2], rev(x_mat[[r]][,3]))
    polygon(x, y.90, col = adjustcolor(col.palette1[r], alpha = 0.1), border=NA) # Add uncertainty polygon
    #    polygon(x, y.80, col = adjustcolor(col.palette1[r], alpha = 0.2), border=NA) # Add uncertainty polygon
    
    lines(y_mat[[r]]$median.reef~y_mat[[r]]$year,col=col.palette2[r],lwd=2)
    points(y_mat[[r]]$median.reef~y_mat[[r]]$year,col='white',pch=21,bg=col.palette2[r],cex=1.5)
  }
  text(y=rep(par("usr")[3]+(max(x_mat_full)-min(x_mat_full))*0.015,nrow(ts)),x=ts$year,ts$n.survs,cex=0.5)
  text(y=par("usr")[3]+(max(x_mat_full)-min(x_mat_full))*0.045,x=par("usr")[2]+0.2,'Surveys',cex=0.6,pos=1,xpd=T,font=1)
  for(r in 1:n.groups){
    lines(y=rep(par("usr")[4]-(max(x_mat_full)-min(x_mat_full))*0.03*r,2),x=c(par("usr")[2]-0.5,par("usr")[2]+0.5),col=col.palette2[r],xpd=T,lwd=2)
  }
  for(r in 1:n.groups){
    text(y=par("usr")[4]-(max(x_mat_full)-min(x_mat_full))*0.03*r,x=par("usr")[2]-2.5,grp[r],col=col.palette2[r],xpd=T,lwd=2,cex=0.6)
  }
  
  dev.off(file.path(path,paste(sprintf("%02d",i),'_',sp,'_',GZ,'.pdf',sep='')))
  dev.off()
}

State_space_timeseries_plot_pdf<- function(sp,GZ,params1,TT,ts,path,i){
  pdf(file.path(path,paste(sprintf("%02d",i),'_',sp,'_',GZ,'.pdf',sep='')),width=8,height=6)
  lambda_mat<- list()  
  for(i in 1:TT){
    reef_coef<- data.frame(p_0.y=NA,p_1.y=NA,p_2.y=NA,p_11.y=NA,p_101.y=NA,p_0.x=NA,p_1.x=NA,p_2.x=NA,p_11.x=NA,p_101.x=NA,lambda.y=NA,lambda.x=NA,iter=seq(1,nrow(params1$c)))
    
    if(ncol(params1$c)==2){
      reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,i])
      reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,i])-plogis(params1$c[,1]-params1$a_yr[,i])
      reef_coef[,3]<-1-plogis(params1$c[,2]-params1$a_yr[,i])
      reef_coef[,4]<- 0
      reef_coef[,5]<- 0
      reef_coef[,6]=plogis(params1$c[,1]-params1$x[,i])
      reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,i])-plogis(params1$c[,1]-params1$x[,i])
      reef_coef[,8]<-1-plogis(params1$c[,2]-params1$x[,i])
      reef_coef[,9]<- 0
      reef_coef[,10]<- 0
      
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    
    if(ncol(params1$c)==3){
      reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,i])
      reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,i])-plogis(params1$c[,1]-params1$a_yr[,i])
      reef_coef[,3]<-plogis(params1$c[,3]-params1$a_yr[,i])-plogis(params1$c[,2]-params1$a_yr[,i])
      reef_coef[,4]<- 1-plogis(params1$c[,3]-params1$a_yr[,i])
      reef_coef[,5]<- 0
      reef_coef[,6]=plogis(params1$c[,1]-params1$x[,i])
      reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,i])-plogis(params1$c[,1]-params1$x[,i])
      reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,i])-plogis(params1$c[,2]-params1$x[,i])
      reef_coef[,9]<- 1-plogis(params1$c[,3]-params1$x[,i])
      reef_coef[,10]<- 0
      
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    
    if(ncol(params1$c)==4){
      reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,i])
      reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,i])-plogis(params1$c[,1]-params1$a_yr[,i])
      reef_coef[,3]<-plogis(params1$c[,3]-params1$a_yr[,i])-plogis(params1$c[,2]-params1$a_yr[,i])
      reef_coef[,4]<- plogis(params1$c[,4]-params1$a_yr[,i])-plogis(params1$c[,3]-params1$a_yr[,i])
      reef_coef[,5]<- 1- plogis(params1$c[,4]-params1$a_yr[,i])
      reef_coef[,6]=plogis(params1$c[,1]-params1$x[,i])
      reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,i])-plogis(params1$c[,1]-params1$x[,i])
      reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,i])-plogis(params1$c[,2]-params1$x[,i])
      reef_coef[,9]<- plogis(params1$c[,4]-params1$x[,i])-plogis(params1$c[,3]-params1$x[,i])
      reef_coef[,10]<- 1-plogis(params1$c[,4]-params1$x[,i])
      
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    lambda_mat[[i]]=reef_coef
  }  
  
  
  x_mat<- data.frame(median.reef=NA,l.95.rf=NA,u.95.rf=NA)
  for(i in 1:TT){
    x_mat[i,1]=median(lambda_mat[[i]]$lambda.x)
    x_mat[i,2]=quantile(lambda_mat[[i]]$lambda.x,0.05)
    x_mat[i,3]=quantile(lambda_mat[[i]]$lambda.x,0.95)
  }
  
  y_mat<- data.frame(year=seq(min(ts$year),max(ts$year)),median.reef=NA,l.95.rf=NA,u.95.rf=NA)
  
  for(i in 1:TT){
    y_mat[i,2]=median(lambda_mat[[i]]$lambda.y)
    y_mat[i,3]=quantile(lambda_mat[[i]]$lambda.y,0.05)
    y_mat[i,4]=quantile(lambda_mat[[i]]$lambda.y,0.95)
  }
  
  par(xpd=T)
  plot(y_mat$median.reef~y_mat$year,type='n',ylim=c(min(x_mat),max(c(max(y_mat[,2]),max(x_mat)))),col='darkblue',bty='l',ylab=expression('Mean count per survey'),xlab='Year',main=paste(sp,GZ,sep=', '))
  lines(x_mat[,1]~y_mat$year,lty=5,lwd=2,col='darkcyan')
  x<- c(y_mat$year, rev(y_mat$year))
  y1<- c(x_mat[,2], rev(x_mat[,3]))
  polygon(x, y1, col = adjustcolor('darkcyan', alpha = 0.1), border=NA) # Add uncertainty polygon
  
  lines(y_mat$median.reef~y_mat$year,col='dodgerblue4',lwd=2)
  points(y_mat$median.reef~y_mat$year,col='white',pch=21,bg='dodgerblue4',cex=1.5)
  text(y=rep(min(x_mat),nrow(ts)),x=ts$year,ts$n.survs,cex=0.5)
  dev.off(file.path(path,paste(sprintf("%02d",i),'_',sp,'_',GZ,'.pdf',sep='')))
  dev.off()
}


ts_reef = function(X){
  occ_by_year<- X %>% group_by(year) %>% summarize(n.occ=sum(occ),n.surv=n(),p.occ=n.occ/n.surv)
  abun_by_year<- X %>% group_by(year,geogr) %>% summarize(site_abun=mean(abund_trans),sd_abund=sd(abund_trans),n.surv=n()) %>% group_by(year) %>% summarize(mean_abund=mean(site_abun),n.survs=sum(n.surv),n.sites=n())
  total_sd<- X %>% group_by(year) %>% summarize(sd=sd(abund_trans))
  
  comb<- left_join(occ_by_year,abun_by_year)
  comb2<- left_join(comb,total_sd)
  ts_dat<- comb2
  return(ts_dat)
}

abund_tranfs<- function(probs){
  sum<- probs[1]*0+probs[2]*1+probs[3]*2+probs[4]*11+probs[5]*101
  return(sum)
}

ord_to_n<- function(x,c){
  abund_x<- numeric(length(x))
  p= matrix(ncol=5,nrow=length(x))
  if(ncol(c)==2){
    p[,1]=plogis(c[,1]-x)
    p[,2]=plogis(c[,2]-x)-plogis(c[,1]-x)
    p[,3]=1-plogis(c[,2]-x)
    p[,4]=0
    p[,5]=0
    for(i in 1:length(x)){
      abund_x[i]=abund_tranfs(p[i,])  
    }
  }
  if(ncol(c)==3){
    p[,1]=plogis(c[,1]-x)
    p[,2]=plogis(c[,2]-x)-plogis(c[,1]-x)
    p[,3]=plogis(c[,3]-x)-plogis(c[,2]-x)
    p[,4]=1-plogis(c[,3]-x)
    p[,5]=0
    for(i in 1:length(x)){
      abund_x[i]=abund_tranfs(p[i,])  
    }
  }
  if(ncol(c)==4){
    p[,1]=plogis(c[,1]-x)
    p[,2]=plogis(c[,2]-x)-plogis(c[,1]-x)
    p[,3]=plogis(c[,3]-x)-plogis(c[,2]-x)
    p[,4]=plogis(c[,4]-x[,i])-plogis(c[,3]-x[,i])
    p[,5]=1-plogis(c[,4]-x[,i])
    for(i in 1:length(x)){
      abund_x[i]=abund_tranfs(p[i,])  
    }
  }
  return(abund_x)
}

population_growth_rates<- function(params1,TT){
lambda_mat<- list()  
for(i in 1:TT){
  reef_coef<- data.frame(p_0.y=NA,p_1.y=NA,p_2.y=NA,p_11.y=NA,p_101.y=NA,p_0.x=NA,p_1.x=NA,p_2.x=NA,p_11.x=NA,p_101.x=NA,lambda.y=NA,lambda.x=NA,iter=seq(1,nrow(params$c)))
  
  if(ncol(params1$c)==2){
    reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,i])
    reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,i])-plogis(params1$c[,1]-params1$a_yr[,i])
    reef_coef[,3]<-1-plogis(params1$c[,2]-params1$a_yr[,i])
    reef_coef[,4]<- 0
    reef_coef[,5]<- 0
    reef_coef[,6]=plogis(params1$c[,1]-params1$x[,i])
    reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,i])-plogis(params1$c[,1]-params1$x[,i])
    reef_coef[,8]<-1-plogis(params1$c[,2]-params1$x[,i])
    reef_coef[,9]<- 0
    reef_coef[,10]<- 0
    
    reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
    reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
  }
  
  if(ncol(params1$c)==3){
    reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,i])
    reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,i])-plogis(params1$c[,1]-params1$a_yr[,i])
    reef_coef[,3]<-plogis(params1$c[,3]-params1$a_yr[,i])-plogis(params1$c[,2]-params1$a_yr[,i])
    reef_coef[,4]<- 1-plogis(params1$c[,3]-params1$a_yr[,i])
    reef_coef[,5]<- 0
    reef_coef[,6]=plogis(params1$c[,1]-params1$x[,i])
    reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,i])-plogis(params1$c[,1]-params1$x[,i])
    reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,i])-plogis(params1$c[,2]-params1$x[,i])
    reef_coef[,9]<- 1-plogis(params1$c[,3]-params1$x[,i])
    reef_coef[,10]<- 0
    
    reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
    reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
  }
  
  if(ncol(params1$c)==4){
    reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,i])
    reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,i])-plogis(params1$c[,1]-params1$a_yr[,i])
    reef_coef[,3]<-plogis(params1$c[,3]-params1$a_yr[,i])-plogis(params1$c[,2]-params1$a_yr[,i])
    reef_coef[,4]<- plogis(params1$c[,4]-params1$a_yr[,i])-plogis(params1$c[,3]-params1$a_yr[,i])
    reef_coef[,5]<- 1- plogis(params1$c[,4]-params1$a_yr[,i])
    reef_coef[,6]=plogis(params1$c[,1]-params1$x[,i])
    reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,i])-plogis(params1$c[,1]-params1$x[,i])
    reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,i])-plogis(params1$c[,2]-params1$x[,i])
    reef_coef[,9]<- plogis(params1$c[,4]-params1$x[,i])-plogis(params1$c[,3]-params1$x[,i])
    reef_coef[,10]<- 1-plogis(params1$c[,4]-params1$x[,i])
    
    reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
    reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
  }
  lambda_mat[[i]]=reef_coef
}  


x_mat<- matrix(nrow=nrow(lambda_mat[[1]]),ncol=TT)
for(i in 1:TT){
  x_mat[,i]=lambda_mat[[i]]$lambda.x
}
r_mat<- matrix(nrow=nrow(lambda_mat[[1]]),ncol=TT-1)
for(i in 1:nrow(r_mat)){
  for(t in 1:TT-1)
    r_mat[i,t]=log(x_mat[i,t+1])-log(x_mat[i,t])
}
return(r_mat)
}

gm_mean_r<- function(x){
  prod(exp(x))^(1/length(x))
}

population_growth_rates_multi<- function(params1,TT,n.groups){
lambda_mat<- list()
x_mat<- list()
r_mat<- list()
for(r in 1:n.groups){
  ts_estimates<- list()
  for(i in 1:TT){
    reef_coef<- data.frame(p_0.y=NA,p_1.y=NA,p_2.y=NA,p_11.y=NA,p_101.y=NA,p_0.x=NA,p_1.x=NA,p_2.x=NA,p_11.x=NA,p_101.x=NA,lambda.y=NA,lambda.x=NA,iter=seq(1,nrow(params$c)))
    
    if(ncol(params1$c)==2){
      reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,,r][,i])
      reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,,r][,i])-plogis(params1$c[,1]-params1$a_yr[,,r][,i])
      reef_coef[,3]<-1-plogis(params1$c[,2]-params1$a_yr[,,r][,i])
      reef_coef[,4]<- 0
      reef_coef[,5]<- 0
      reef_coef[,6]=plogis(params1$c[,1]-params1$x[,,r][,i])
      reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,,r][,i])-plogis(params1$c[,1]-params1$x[,,r][,i])
      reef_coef[,8]<-1-plogis(params1$c[,2]-params1$x[,,r][,i])
      reef_coef[,9]<- 0
      reef_coef[,10]<- 0
      
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    
    if(ncol(params1$c)==3){
      reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,,r][,i])
      reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,,r][,i])-plogis(params1$c[,1]-params1$a_yr[,,r][,i])
      reef_coef[,3]<-plogis(params1$c[,3]-params1$a_yr[,,r][,i])-plogis(params1$c[,2]-params1$a_yr[,,r][,i])
      reef_coef[,4]<- 1-plogis(params1$c[,3]-params1$a_yr[,,r][,i])
      reef_coef[,5]<- 0
      reef_coef[,6]=plogis(params1$c[,1]-params1$x[,,r][,i])
      reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,,r][,i])-plogis(params1$c[,1]-params1$x[,,r][,i])
      reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,,r][,i])-plogis(params1$c[,2]-params1$x[,,r][,i])
      reef_coef[,9]<- 1-plogis(params1$c[,3]-params1$x[,,r][,i])
      reef_coef[,10]<- 0
      
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    
    if(ncol(params1$c)==4){
      reef_coef[,1]<- plogis(params1$c[,1]-params1$a_yr[,,r][,i])
      reef_coef[,2]<-plogis(params1$c[,2]-params1$a_yr[,,r][,i])-plogis(params1$c[,1]-params1$a_yr[,,r][,i])
      reef_coef[,3]<-plogis(params1$c[,3]-params1$a_yr[,,r][,i])-plogis(params1$c[,2]-params1$a_yr[,,r][,i])
      reef_coef[,4]<- plogis(params1$c[,4]-params1$a_yr[,,r][,i])-plogis(params1$c[,3]-params1$a_yr[,,r][,i])
      reef_coef[,5]<- 1- plogis(params1$c[,4]-params1$a_yr[,,r][,i])
      reef_coef[,6]=plogis(params1$c[,1]-params1$x[,,r][,i])
      reef_coef[,7]<-plogis(params1$c[,2]-params1$x[,,r][,i])-plogis(params1$c[,1]-params1$x[,,r][,i])
      reef_coef[,8]<-plogis(params1$c[,3]-params1$x[,,r][,i])-plogis(params1$c[,2]-params1$x[,,r][,i])
      reef_coef[,9]<- plogis(params1$c[,4]-params1$x[,,r][,i])-plogis(params1$c[,3]-params1$x[,,r][,i])
      reef_coef[,10]<- 1-plogis(params1$c[,4]-params1$x[,,r][,i])
      
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    ts_estimates[[i]]=reef_coef
  }
  lambda_mat[[r]]=ts_estimates
} 


for(r in 1:n.groups){  
  x_mat[[r]]<- matrix(nrow=nrow(lambda_mat[[1]][[1]]),ncol=TT)
  for(i in 1:TT){
    x_mat[[r]][,i]=lambda_mat[[r]][[i]]$lambda.x
  }
  r_mat[[r]]<- matrix(nrow=nrow(lambda_mat[[1]][[1]]),ncol=TT-1)
  for(i in 1:nrow(r_mat[[r]])){
    for(t in 1:TT-1)
      r_mat[[r]][i,t]=log(x_mat[[r]][i,t+1])-log(x_mat[[r]][i,t])
  }
}
return(r_mat) 
}
