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
  
  TempDat_site2<- TempDat4 %>% group_by(geogr) %>% summarize(n=n(),n.yr=n_distinct(year)) %>% subset(n>9&n.yr>3)
  TempDat5<- subset(TempDat4, geogr %in% TempDat_site2$geogr)

  survs<- TempDat5 %>% distinct(formid,.keep_all = T) %>% group_by(geogr,day,month,year) %>% summarize(n=n()) %>% mutate(site_dmy=paste(geogr,day,month,year,sep='')) %>% subset(n>1)
  TempDat6<- TempDat5 %>% mutate(site_dmy=ifelse(paste(geogr,day,month,year,sep='') %in% survs$site_dmy,paste(geogr,day,month,year,sep=''),'baseline'))
  
  return(TempDat6)
}


scaled_mv_timeseries_plot2<- function(ts2,sp,reg,params,path,TT,yrs.rvc,yr.start,yr.end){
  
  #Extract parameters
  cuts=as.data.frame(params)[grepl('cut',colnames(params))]
  x=as.data.frame(params)[grepl('x',colnames(params))] #remove x0
  x1=x[,1:TT]
  x2=x[,c(TT+1):ncol(x)]
  y1=as.data.frame(params)[grepl('est1',colnames(params))]; 
  y2=as.data.frame(params)[grepl('a_yr2',colnames(params))]; 
  cor=as.data.frame(params)$'Cor_t[1,2]' 
  
  lambda_mat<- list()  
  for(i in 1:TT){
    reef_coef<- data.frame(p_0.y=NA,p_1.y=NA,p_2.y=NA,p_11.y=NA,p_101.y=NA,p_0.x=NA,p_1.x=NA,p_2.x=NA,p_11.x=NA,p_101.x=NA,lambda.y=NA,lambda.x=NA,iter=seq(1:nrow(x)))
    
    if(ncol(cuts)==2){
      reef_coef[,1]<- plogis(cuts[,1]-y2[,i])
      reef_coef[,2]<-plogis(cuts[,2]-y2[,i])-plogis(cuts[,1]-y2[,i])
      reef_coef[,3]<-1-plogis(cuts[,2]-y2[,i])
      reef_coef[,4]<- 0
      reef_coef[,5]<- 0
      reef_coef[,6]=plogis(cuts[,1]-(x2[,i]))
      reef_coef[,7]<-plogis(cuts[,2]-(x2[,i]))-plogis(cuts[,1]-(x2[,i]))
      reef_coef[,8]<-1-plogis(cuts[,2]-(x2[,i]))
      reef_coef[,9]<- 0
      reef_coef[,10]<- 0
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    
    if(ncol(cuts)==3){
      reef_coef[,1]<- plogis(cuts[,1]-y2[,i])
      reef_coef[,2]<-plogis(cuts[,2]-y2[,i])-plogis(cuts[,1]-y2[,i])
      reef_coef[,3]<-plogis(cuts[,3]-y2[,i])-plogis(cuts[,2]-y2[,i])
      reef_coef[,4]<- 1-plogis(cuts[,3]-y2[,i])
      reef_coef[,5]<- 0
      reef_coef[,6]=plogis(cuts[,1]-(x2[,i]))
      reef_coef[,7]<-plogis(cuts[,2]-(x2[,i]))-plogis(cuts[,1]-(x2[,i]))
      reef_coef[,8]<-plogis(cuts[,3]-(x2[,i]))-plogis(cuts[,2]-(x2[,i]))
      reef_coef[,9]<- 1-plogis(cuts[,3]-(x2[,i]))
      reef_coef[,10]<- 0
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    
    if(ncol(cuts)==4){
      reef_coef[,1]<- plogis(cuts[,1]-y2[,i])
      reef_coef[,2]<-plogis(cuts[,2]-y2[,i])-plogis(cuts[,1]-y2[,i])
      reef_coef[,3]<-plogis(cuts[,3]-y2[,i])-plogis(cuts[,2]-y2[,i])
      reef_coef[,4]<- plogis(cuts[,4]-y2[,i])-plogis(cuts[,3]-y2[,i])
      reef_coef[,5]<- 1- plogis(cuts[,4]-y2[,i])
      reef_coef[,6]<-plogis(cuts[,1]-(x2[,i]))
      reef_coef[,7]<-plogis(cuts[,2]-(x2[,i]))-plogis(cuts[,1]-(x2[,i]))
      reef_coef[,8]<-plogis(cuts[,3]-(x2[,i]))-plogis(cuts[,2]-(x2[,i]))
      reef_coef[,9]<- plogis(cuts[,4]-(x2[,i]))-plogis(cuts[,3]-(x2[,i]))
      reef_coef[,10]<- 1-plogis(cuts[,4]-(x2[,i]))
      reef_coef[,11]<- apply(reef_coef[,1:5],1,abund_tranfs)
      reef_coef[,12]<- apply(reef_coef[,6:10],1,abund_tranfs)
    }
    lambda_mat[[i]]=reef_coef
    
  }  
  median.x=median(do.call(rbind, lapply(lambda_mat, data.frame, stringsAsFactors=FALSE))$lambda.x)
  median.y=median(do.call(rbind, lapply(lambda_mat, data.frame, stringsAsFactors=FALSE))$lambda.y)
  
  x_mat<- data.frame(median.rvc=NA,l.95=NA,u.95=NA,median.reef=NA,l.95.rf=NA,u.95.rf=NA)
  for(i in 1:TT){
    x_mat[i,1]=median(x1[,i]/median(as.vector(t(x1))))
    x_mat[i,2]=quantile(x1[,i]/median(as.vector(t(x1))),0.05)
    x_mat[i,3]=quantile(x1[,i]/median(as.vector(t(x1))),0.95)
    x_mat[i,4]=median(lambda_mat[[i]]$lambda.x/median.x)
    x_mat[i,5]=quantile(lambda_mat[[i]]$lambda.x/median.x,0.05)
    x_mat[i,6]=quantile(lambda_mat[[i]]$lambda.x/median.x,0.95)
  }
  
  y_mat_rvc<- data.frame(year=yrs.rvc,median.rvc=NA,l.95.rvc=NA,u.95.rvc=NA)
  y_mat_reef<- data.frame(year=seq(yr.start,yr.end),median.reef=NA,l.95.rf=NA,u.95.rf=NA)
  for(i in 1:length(yrs.rvc)){
    y_mat_rvc[i,2]=median(y1[,i]/median(as.vector(t(y1))))
    y_mat_rvc[i,3]=quantile(y1[,i]/median(as.vector(t(y1))),0.05)
    y_mat_rvc[i,4]=quantile(y1[,i]/median(as.vector(t(y1))),0.95)
  }
  
  for(i in 1:TT){
    y_mat_reef[i,2]=median(lambda_mat[[i]]$lambda.y/median.y)
    y_mat_reef[i,3]=quantile(lambda_mat[[i]]$lambda.y/median.y,0.05)
    y_mat_reef[i,4]=quantile(lambda_mat[[i]]$lambda.y/median.y,0.95)
  }
  y_mat<- full_join(y_mat_reef,y_mat_rvc)
  
  pdf(file.path(path,paste(paste(sp,reg,sep='_'),'.pdf',sep='')),width=8,height=6)
  par(xpd=T,mar=c(5,5,2,2))
  plot(y_mat$median.rvc~y_mat$year,type='n',ylim=c(min(x_mat),max(c(max(y_mat[,2]),max(x_mat)))),col='darkblue',bty='l',ylab=expression('Scaled Relative Abundance'),xlab='Year',main=paste(sp,reg,sep=' - '))
  
  x<- c(y_mat$year, rev(y_mat$year))
  y1<- c(x_mat[,2], rev(x_mat[,3]))
  y2<-  c(x_mat[,5], rev(x_mat[,6]))
  polygon(x, y1, col = adjustcolor('dodgerblue4', alpha = 0.2), border=NA) # Add uncertainty polygon
  polygon(x, y2, col = adjustcolor('firebrick4', alpha = 0.2), border=NA) # Add uncertainty polygon
  lines(x_mat[,1]~y_mat$year,lty=5,lwd=2,col='dodgerblue3')
  lines(x_mat[,4]~y_mat$year,lty=5,lwd=2,col='firebrick')
  
  lines(y_mat$median.rvc~y_mat$year,col='white',lwd=3)
  lines(y_mat$median.rvc~y_mat$year,col='dodgerblue4',lwd=2)
  points(y_mat$median.rvc~y_mat$year,col='white',pch=21,bg='dodgerblue4',cex=1.5,lwd=1)
  
  lines(y_mat$median.reef~y_mat$year,col='white',lwd=3)
  lines(y_mat$median.reef~c(seq(yr.start,yr.end)),col='firebrick4',lwd=2)
  points(y_mat$median.reef~c(seq(yr.start,yr.end)),col='white',pch=21,bg='firebrick4',cex=1.5,lwd=1)
  
  legend(yr.start-0.5,c(max(c(max(y_mat[,2]),max(x_mat)))*1.1),c('RVC','REEF'),text.col=c('dodgerblue4','firebrick4'),bty='n',y.intersp=1.4)
  u <- par("usr")
  v <- c(
    grconvertX(u[1:2], "user", "ndc"),
    grconvertY(u[3:4], "user", "ndc")
  )
  v <- c(v[1]*1.025, (v[1]+v[2])/3, (v[3]+v[4])/1.4, v[4] )
  par( fig=v, new=TRUE, mar=c(0,0,0,0),xpd=T)
  hist(cor,xlab='',ylab='',axes=F,main='',breaks=30,col=adjustcolor('darkgray',alpha.f=0.6),border='white',xlim=c(-1,1))
  abline(v=median(cor))
  axis(side=1,at=c(-1,-0.5,0,0.5,1),labels=c('-1','','0','','1'),cex.axis=0.8,col=adjustcolor('black',alpha.f=0.6),lwd=0.8,lwd.ticks=0.8,padj=-0.5)
  mtext(expression(rho),side=3,line=-0.2)
  
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

ts_rvc = function(x){ #Takes the output from the previous function
    x1=x
    x2= x1 %>% group_by(YEAR) %>% summarize(n.occ=sum(OCC),n.surv=n(),p.occ=n.occ/n.surv,sp=unique(SPECIES_CD))
    x3= x1 %>% group_by(YEAR,PRIMARY_SAMPLE_UNIT) %>% summarize(psu_abund=mean(NUM2)) %>% group_by(YEAR) %>% summarize(mean_abund=mean(psu_abund),sd_abund=sd(psu_abund))
    x4=left_join(x2,x3) %>% complete(YEAR=seq(min(x3$YEAR),max(x3$YEAR)))
  return(x4)
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

`%notin%` <- function(x, y) {
  !(x %in% y)
}
