


## Preliminaries
rm(list=ls())

# Change working directory to where you've stored ZTRAX
path<- "P:/Peter/Hedonics/Groundwater/"
#install.packages("dplyr", repos = "http://mran.revolutionanalytics.com")
## This function will check if a package is installed, and if not, install it
pkgTest <- function(x) {
  if (!require(x, character.only = TRUE))
  {
    install.packages(x, dep = TRUE)
    if(!require(x, character.only = TRUE)) stop("Package not found")
  }
}

## These lines load the required packages
packages <- c("readxl","Hmisc","DescTools","qgam","quantreg","sphet","mgcv","McSpatial","pastecs","rdd","Matrix","psych","xtable","splines","ck37r","data.table","matrixStats","tmle","xgboost", "MatchIt","gtools","statar","foreign","multiwayvcov","lmtest","readstata13","xlsx", "data.table","doSNOW","parallel","compare","doParallel","devtools","foreach","spdep","reshape2","sm","plyr","utils","tcltk","geosphere", "matrixcalc", "dplyr","ExPosition", "randomForest","lfe", "hdm", "rdrobust", "stargazer", "ggplot2", "outliers","rpart","e1071")
lapply(packages, pkgTest)
packages <- c("readxl","Hmisc","grf","classInt","RColorBrewer","rgdal","DescTools","qgam","quantreg","sphet","mgcv","McSpatial","pastecs","rdd","Matrix","psych","xtable","splines","ck37r","data.table","matrixStats","tmle","xgboost", "MatchIt","gtools","statar","foreign","multiwayvcov","lmtest","readstata13","xlsx", "data.table","doSNOW","parallel","compare","doParallel","devtools","foreach","spdep","reshape2","sm","plyr","utils","tcltk","geosphere", "matrixcalc", "dplyr","ExPosition", "randomForest","lfe", "hdm", "rdrobust", "stargazer", "ggplot2", "outliers","rpart","e1071")
lapply(packages, pkgTest)

#library(statar)


## These lines set several options
options(scipen = 999) # Do not print scientific notation
options(stringsAsFactors = FALSE) ## Do not load strings as factors

memory.limit(10000000000000)

NPL<-readRDS(paste(path,'NPLfull.rds', sep=""), refhook = NULL)

sample<-readRDS(paste(path,'baj.rds', sep=""),refhook = NULL)


if(TRUE){
  dist<-c('10k','8k','6k','4k','2k')#,'1k','500m')
  
  #dist<-c('4k','2k')#,'1k','500m')
  
  laglead<-c("")
  llc<-laglead[1]
  treatl<-c('TATE','MUATE','WLATE')
  sitelist<-names(table(sample$closestsite))
  psite<-sitelist[1]
  samplefull<-sample
  samplefull<-samplefull[preprice>0,]
  samplefull$treatst<-ifelse(samplefull$postdel*samplefull$tgdel>0,1,0)
  buf<-1
  #samplefull$buffer<-ifelse(samplefull$date-odNPL$date[psite]+(buf*365)>0&samplefull$date-odNPL$date[psite]<0,1,0)
  samplefull$buffer<-ifelse(abs(samplefull$timetodel)-(buf*365)<0,1,0)
  samplefull<-samplefull[buffer<1,]
  samplefull$delprocess<-ifelse(samplefull$timetodel<0 &samplefull$postcc>0,1,0)
  samplefull<-samplefull[delprocess<1,]
  samplefull$timetotreat<-samplefull$timetodel
  samplefull$lsite<-samplefull$closestsite
  samplefull$demlogprice<-demeanlist(samplefull$logprice,
                                     list(as.factor(samplefull$PropertyAddressCensusTractAndBlock)))
  # sample$logprice<-sample$logprice.x
  #Total Average Treatment Effect
  
  samplefull$aftpropnpl<-ifelse(samplefull$time_to_prop>01,0)
  
  #sample$aftfinalnpl<-sample[[paste0('aftfinalnpl',fnplsite)]]
  #sample$timefinalnplfe<-sample[[paste0('timefinalnplfe',fnplsite)]]
  samplefull$treatControlsComplete<-ifelse(samplefull$time_to_cc>0,1,0)
  
  #Municipal ATE
  samplefull$treatdgwMU<- samplefull$treatst* samplefull$WaterStndCode.fMU
  samplefull$treatgroupMU<-samplefull$tgdel * samplefull$WaterStndCode.fMU
  samplefull$controlMU<-samplefull$control*samplefull$WaterStndCode.fMU
  
  samplefull$sample.MUATE<-samplefull$control+samplefull$treatgroupMU
  
  #Municipal ATE
  samplefull$treatdgwWL<- samplefull$treatst* samplefull$WaterStndCode.fWL
  samplefull$treatgroupWL<-samplefull$tgdel * samplefull$WaterStndCode.fWL
  samplefull$controlWL<-samplefull$control*samplefull$WaterStndCode.fWL
  
  samplefull$sample.WLATE<-samplefull$control+samplefull$treatgroupWL
  
  #psite<-2
  samplefull$timetotreat<-samplefull$timetotreat/365
  samplefull<-samplefull[abs(timetotreat)-10<0,]

}


dist<-c('10k','8k','6k','4k','2k')#,'1k','500m')

#dist<-c('4k','2k')#,'1k','500m')

laglead<-c("")
treatl<-c('TATE','MUATE','WLATE')
#sample<-readRDS(paste0(path,"pretreatlag.rds"), refhook = NULL)
#sample<-readRDS(paste0(path,"fullcen.rds"), refhook = NULL)
sample$treatmentgroup<-sample$tgdel
samplefull<-sample
quant<-10
qcut<-cut2(samplefull$preprice, g=quant, onlycuts = TRUE)

sumvar<-c("price","sqfeet","TotalRooms","YearBuilt","FullBath")
ts<-sample[treatmentgroup==1,c("price","sqfeet","YearBuilt","RecordingDate","FullBath","LotSizeSquareFeet")]
ts$year <- as.numeric(as.character(factor(format(as.Date(ts$RecordingDate),'%Y'))))
cs<-sample[treatmentgroup==0,c("price","sqfeet","YearBuilt","RecordingDate","FullBath","LotSizeSquareFeet")]
cs$year <- as.numeric(as.character(factor(format(as.Date(cs$RecordingDate),'%Y'))))

sum25<-rbind("25th percentile",quantile(ts$price)[2],quantile(ts$sqfeet)[2],quantile(ts$YearBuilt)[2],quantile(ts$year)[2],quantile(ts$FullBath)[2],quantile(ts$LotSizeSquareFeet)[2],
             quantile(cs$price)[2],quantile(cs$sqfeet)[2],quantile(cs$YearBuilt)[2],quantile(cs$year)[2],quantile(cs$FullBath)[2],quantile(cs$LotSizeSquareFeet)[2])
summed<-rbind("Median",median(ts$price),median(ts$sqfeet),median(ts$YearBuilt),median(ts$year),median(ts$FullBath),median(ts$LotSizeSquareFeet),
              median(cs$price),median(cs$sqfeet),median(cs$YearBuilt),median(cs$year),median(cs$FullBath),median(cs$LotSizeSquareFeet))
summean<-rbind("Mean",floor(mean(ts$price)),floor(mean(ts$sqfeet)),floor(mean(ts$YearBuilt)),floor(mean(ts$year)),floor(mean(ts$FullBath)),floor(mean(ts$LotSizeSquareFeet)),
               floor(mean(cs$price)),floor(mean(cs$sqfeet)),floor(mean(cs$YearBuilt)),floor(mean(cs$year)),floor(mean(cs$FullBath)),floor(mean(cs$LotSizeSquareFeet)))
sum75<-rbind("75th percentile",quantile(ts$price)[4],quantile(ts$sqfeet)[4],quantile(ts$YearBuilt)[4],quantile(ts$year)[4],quantile(ts$FullBath)[4],quantile(ts$LotSizeSquareFeet)[4],
             quantile(cs$price)[4],quantile(cs$sqfeet)[4],quantile(cs$YearBuilt)[4],quantile(cs$year)[4],quantile(cs$FullBath)[4],quantile(cs$LotSizeSquareFeet)[4])
tdiff<-rbind("T-Test",base::round(as.numeric(t.test(ts$price,cs$price)["statistic"]),digits=2),
             base::round(as.numeric(t.test(ts$sqfeet,cs$sqfeet)["statistic"]),digits=2),
             base::round(as.numeric(t.test(ts$year,cs$year)["statistic"]),digits=2),
             base::round(as.numeric(t.test(ts$YearBuilt,cs$YearBuilt)["statistic"]),digits=2),
             base::round(as.numeric(t.test(ts$FullBath,cs$FullBath)["statistic"]),digits=2),
             base::round(as.numeric(t.test(ts$LotSizeSquareFeet,cs$LotSizeSquareFeet)["statistic"]),digits=2),
             base::round(as.numeric(t.test(cs$price,ts$price)["statistic"]),digits=2),
             base::round(as.numeric(t.test(cs$sqfeet,ts$sqfeet)["statistic"]),digits=2),
             base::round(as.numeric(t.test(cs$year,ts$year)["statistic"]),digits=2),
             base::round(as.numeric(t.test(cs$YearBuilt,ts$YearBuilt)["statistic"]),digits=2),
             base::round(as.numeric(t.test(cs$FullBath,ts$FullBath)["statistic"]),digits=2),
             base::round(as.numeric(t.test(cs$LotSizeSquareFeet,ts$LotSizeSquareFeet)["statistic"]),digits=2))

rlab<-rbind("     ","Price","Square Feet","Year Built","Year","Full Bathroom","Lot Size","Price","Square Feet","Year Built","Year","Full Bathroom","Lot Size")
rlab2<-rbind("     ","Treatment","Group","  ","   ","  ","   ","Control","Group","  ","   ","  ","   ")
sumtab<-cbind(rlab2,rlab,sum25,summed,summean,sum75,tdiff)
rownames(sumtab)<-NULL
colnames(sumtab)<-NULL
sumtab

xtab<-xtable(sumtab)
align(xtab) <- "rl|l|rrrr|r"
print.xtable(xtab,include.rownames=FALSE,  hline.after = c(0,1,7,dim(sumtab)[1]),
             include.colnames=FALSE, sanitize.text.function = identity,
             caption = "Summary Statistics", 
             label = "tab:summary",
             type="latex", file=paste0(path,'latex/sumtab.tex'))


sumvar<-c("price","sqfeet","TotalRooms","YearBuilt","FullBath")
pts<-sample[treatmentgroup==1&timetotreat<0,c("price","sqfeet","YearBuilt","RecordingDate","FullBath","LotSizeSquareFeet")]
pts$year <- as.numeric(as.character(factor(format(as.Date(pts$RecordingDate),'%Y'))))
pcs<-sample[treatmentgroup==0&timetotreat<0,c("price","sqfeet","YearBuilt","RecordingDate","FullBath","LotSizeSquareFeet")]
pcs$year <- as.numeric(as.character(factor(format(as.Date(pcs$RecordingDate),'%Y'))))
ats<-sample[treatmentgroup==1&timetotreat>0,c("price","sqfeet","YearBuilt","RecordingDate","FullBath","LotSizeSquareFeet")]
ats$year <- as.numeric(as.character(factor(format(as.Date(ats$RecordingDate),'%Y'))))
acs<-sample[treatmentgroup==0&timetotreat>0,c("price","sqfeet","YearBuilt","RecordingDate","FullBath","LotSizeSquareFeet")]
acs$year <- as.numeric(as.character(factor(format(as.Date(acs$RecordingDate),'%Y'))))


sum25<-rbind("25th percentile",quantile(pts$price)[2],quantile(pts$sqfeet)[2],quantile(pts$YearBuilt)[2],quantile(pts$year)[2],quantile(pts$LotSizeSquareFeet)[2],
             quantile(pcs$price)[2],quantile(pcs$sqfeet)[2],quantile(pcs$YearBuilt)[2],quantile(pcs$year)[2],quantile(pcs$LotSizeSquareFeet)[2],
             quantile(ats$price)[2],quantile(ats$sqfeet)[2],quantile(ats$YearBuilt)[2],quantile(ats$year)[2],quantile(ats$LotSizeSquareFeet)[2],
             quantile(acs$price)[2],quantile(acs$sqfeet)[2],quantile(acs$YearBuilt)[2],quantile(acs$year)[2],quantile(acs$LotSizeSquareFeet)[2])
summed<-rbind("Median",median(pts$price),median(pts$sqfeet),median(pts$YearBuilt),median(pts$year),median(pts$LotSizeSquareFeet),
              median(pcs$price),median(pcs$sqfeet),median(pcs$YearBuilt),median(pcs$year),median(pcs$LotSizeSquareFeet),
              median(ats$price),median(ats$sqfeet),median(ats$YearBuilt),median(ats$year),median(ats$LotSizeSquareFeet),
              median(acs$price),median(acs$sqfeet),median(acs$YearBuilt),median(acs$year),median(acs$LotSizeSquareFeet))
summean<-rbind("Mean",mean(pts$price),mean(pts$sqfeet),mean(pts$YearBuilt),mean(pts$year),mean(pts$LotSizeSquareFeet),
               mean(pcs$price),mean(pcs$sqfeet),mean(pcs$YearBuilt),mean(pcs$year),mean(pcs$LotSizeSquareFeet),
               mean(ats$price),mean(ats$sqfeet),mean(ats$YearBuilt),mean(ats$year),mean(ats$LotSizeSquareFeet),
               mean(acs$price),mean(acs$sqfeet),mean(acs$YearBuilt),mean(acs$year),mean(acs$LotSizeSquareFeet))
sum75<-rbind("75th percentile",quantile(pts$price)[4],quantile(pts$sqfeet)[4],quantile(pts$YearBuilt)[4],quantile(pts$year)[4],quantile(pts$LotSizeSquareFeet)[4],
             quantile(pcs$price)[4],quantile(pcs$sqfeet)[4],quantile(pcs$YearBuilt)[4],quantile(pcs$year)[4],quantile(pcs$LotSizeSquareFeet)[4],
             quantile(ats$price)[4],quantile(ats$sqfeet)[4],quantile(ats$YearBuilt)[4],quantile(ats$year)[4],quantile(ats$LotSizeSquareFeet)[4],
             quantile(acs$price)[4],quantile(acs$sqfeet)[4],quantile(acs$YearBuilt)[4],quantile(acs$year)[4],quantile(acs$LotSizeSquareFeet)[4])
#rlab<-rbind("     ","Price","Square Feet","Year Built","Year","Lot Size","Price","Square Feet","Year Built","Year","Full Bathroom","Lot Size")
#rlab2<-rbind("     ","Treatment","Group","  ","   ","  ","   ","Control","Group","  ","   ","  ","   ")

sumtab1<-cbind(sum25,summed,summean,sum75)


#sample$date<-sample$date.x

#dist<-c('10k','8k','6k','5k','4k','3k','2k')#,'1k','500m')
dist<-c('8k','6k','4k','2k')#,'1k','500m')
dist<-c('10k','8k','6k','4k','2k')#,'1k','500m')

#dist<-c('4k','2k')#,'1k','500m')

laglead<-c("")
treatl<-c('TATE','MUATE','WLATE')
#di<-5
#ll<-1

#for(buf in 1:2){


specl<-c("i","a","s","y","")
#matrices

for(i in c("lm","gam","sp","match")){
  for(j in c("t","wl","mu")){
    for(k in c("did")){
      for(treat in treatl){
        for(spec in specl){
          assign(paste0('betas.',i,'.',j,'.',k,'.',treat,'.',spec),matrix(ncol = length(dist),nrow=quant))
          assign(paste0('ses.',i,'.',j,'.',k,'.',treat,'.',spec),matrix(ncol = length(dist),nrow=quant))
          assign(paste0('ps.',i,'.',j,'.',k,'.',treat,'.',spec),matrix(ncol = length(dist),nrow=quant))
          
          assign(paste0('betas.',i,'.',k),matrix(ncol = length(dist),nrow=length(treatl)))
          assign(paste0('ses.',i,'.',k),matrix(ncol = length(dist),nrow=length(treatl)))
          assign(paste0('ps.',i,'.',k),matrix(ncol = length(dist),nrow=length(treatl)))
          
        }
      }
    }
  }
}
di<-1
treat<-1
#psite<-21

for(treat in  1:length(treatl)){
  
  treatc<-treatl[[treat]]
  #for(ll in 1:length(laglead)){
  for(di in 1:length(dist)){
    for(match in c("","match")){
      #match<-""
      
      dic<-dist[[di]]
      ll<-1
      llc<-laglead[[ll]]
      sample<-samplefull[samplefull[[paste0('dist',dic)]]>0,]
      #sample<-sample[presstatusd==0,]
      
      if(match=="match"){
        #LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + sqfeet+day+prediffdate+predate
        samplem<-sample[,c("TransId","date","treatmentgroup","LotSizeSquareFeet", "sqfeet","YearBuilt", "FullBath","predate")]
        mdm<- matchit(treatmentgroup~LotSizeSquareFeet+ YearBuilt + FullBath +sqfeet+predate,data=samplem, method = "nearest", distance = "mahalanobis")
        
        mdm.treat <- match.data(mdm, group = "treat")
        mdm.control <- match.data(mdm, group = "control")
        mdm<-rbind(mdm.treat,mdm.control)
        mdm<-mdm[,c("TransId","date")]
        
        mdm.full<-merge(mdm,sample,all.x = TRUE, by = "TransId")
        sample<-mdm.full
      }
      
      #sample$treatst<-sample$treatst-sample$presstatusd
      #sample<-sample[treatst>0,]
      
      # sample$logprice<-sample$logprice.x
      if(treatc=='TATE'){
        #Total Average Treatment Effect
        
      }
      if(treatc=='MUATE'){
        #Municipal ATE
        
        sample<-subset(sample, sample.MUATE==1)
      }
      if(treatc=='WLATE'){
        #Municipal ATE
        
        sample<-subset(sample, sample.WLATE==1)
      }
      
      if(mean(sample$treatst)>0){
        
        #TATE<-as.formula(logprice ~ treatdgw+ treatmentgroup+data.matrix(treatgroupm)+ data.matrix(timefe)+ data.matrix(year[,3:25])+ as.factor(HHID))
        
        sdf<-20
        lat<-sample$PropertyAddressLatitude
        long<-sample$PropertyAddressLongitude
        splat<-bs(lat, df = sdf)
        splong<-bs(long, df = sdf)
        spint<-model.matrix(~splat:splong)
        
        spTATE<-cbind(splat,splong,spint,lat,long)
        
        #Wst<-readRDS(paste0(path,"Wmat",treatc,dic,".rds"), refhook = NULL)
        #sample$wlp<-Wst%*%sample$prelogprice
        #sample$wp<-Wst%*%sample$preprice
        
        sample$preyear<-as.numeric(format(as.Date(sample$predate),"%Y"))
        sample$premonth<-as.numeric(format(as.Date(sample$predate),"%m"))
        sample$year<-as.numeric(format(as.Date(sample$RecordingDate),"%Y"))
        sample$month<-as.numeric(format(as.Date(sample$RecordingDate),"%m"))
        
        sample$quarter<-ifelse(sample$month==1|sample$month==2|sample$month==3,1,0)
        sample$quarter<-ifelse(sample$month==4|sample$month==5|sample$month==6,2,sample$quarter)
        sample$quarter<-ifelse(sample$month==7|sample$month==8|sample$month==9,3,sample$quarter)
        sample$quarter<-ifelse(sample$month==10|sample$month==11|sample$month==12,4,sample$quarter)
        
        sample$prequarter<-ifelse(sample$premonth==1|sample$premonth==2|sample$premonth==3,1,0)
        sample$prequarter<-ifelse(sample$premonth==4|sample$premonth==5|sample$premonth==6,2,sample$prequarter)
        sample$prequarter<-ifelse(sample$premonth==7|sample$premonth==8|sample$premonth==9,3,sample$prequarter)
        sample$prequarter<-ifelse(sample$premonth==10|sample$premonth==11|sample$premonth==12,4,sample$prequarter)
        
        
        sample$prelogpricewlag<-sample$lnlagu20l0sp10c30tp10
        
        sample$tprelprice<-sample$prelogprice*ifelse(sample$treatst==1 & sample$presstatusd==0,1,0)
        sample$treatfin<-ifelse(sample$treatst==1 & sample$presstatusd==0,1,0)
        indx<- factor(as.numeric(cut2(sample$preprice,cuts=qcut,minmax=TRUE)))
        indx2<- factor(as.numeric(cut2(exp(sample$prelogpricewlag),cuts=qcut,minmax=TRUE)))
        #indx<-ifelse(as.numeric(indx)!=as.numeric(levels(indx)[1])&
        #             as.numeric(indx)!=as.numeric(levels(indx)[length(levels(indx))]),"med",indx)
        #indx<-ifelse(indx==1,"bottom",indx)
        #indx<-ifelse(indx==quant,"top",indx)
        
        treatind<-model.matrix(~treatst:indx-1,sample)
        treatind<-treatind[,1:(length(unique(indx)))]
        
        X<-model.matrix(~ treatmentgroup+indx:treatmentgroup+#:indx+#bs(timetotreat,5)+bs(day,5)+#as.factor(round(timetotreat,1))+
                          #treatmentgroup:as.factor(lsite)+indx:treatmentgroup:as.factor(lsite)+#:indx+#bs(timetotreat,5)+bs(day,5)+#as.factor(round(timetotreat,1))+
                          #poly(timetotreat,3)+poly(day,3)+
                          #treatmentgroup*timetotreat+#treatmentgroup*day+
                          LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                          timetotreat+
                          sqfeet+ prediffdate+predate+prelogprice+indx+#presstatusd+
                          #as.factor(floor(timetotreat))+
                          -1,sample)
        
        qr.X <- qr(X, tol=1e-2, LAPACK = FALSE)
        (rnkX <- qr.X$rank)  ## 4 (number of non-collinear columns)
        (keep <- qr.X$pivot[seq_len(rnkX)])
        ## 1 2 4 5 
        X <- X[,keep]
        
        
        Xlag<-model.matrix(~ treatmentgroup+indx:treatmentgroup+#:indx+#bs(timetotreat,5)+bs(day,5)+#as.factor(round(timetotreat,1))+
                             #treatmentgroup:as.factor(lsite)+indx:treatmentgroup:as.factor(lsite)+#:indx+#bs(timetotreat,5)+bs(day,5)+#as.factor(round(timetotreat,1))+
                             #poly(timetotreat,3)+poly(day,3)+
                             #treatmentgroup*timetotreat+#treatmentgroup*day+
                             timetotreat+
                             LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                             sqfeet+ prediffdate+predate+prelogprice+indx+
                             prelogpricewlag+indx2+#as.factor(floor(timetotreat))+
                             -1,sample)
        
        
        qr.Xlag <- qr(Xlag, tol=1e-2, LAPACK = FALSE)
        (rnkXlag <- qr.Xlag$rank)  ## 4 (number of non-collinear columns)
        (keeplag <- qr.Xlag$pivot[seq_len(rnkXlag)])
        ## 1 2 4 5 
        Xlag <- Xlag[,keeplag]
        #sample<-samplefull
        #sample<-mdm.full
        results.lm.t.did.i<-felm(logprice ~treatind+
                                   X|as.factor(cbl)+
                                   as.factor(year):as.factor(quarter)+
                                   as.factor(preyear):as.factor(prequarter)| 0 | lsite:year,sample) #as.factor(HHID)+
        summary(results.lm.t.did.i)
        
        results.lm.t.did.a<-felm(logprice ~treatind+
                                   X|as.factor(cbg)+
                                   as.factor(year):as.factor(quarter)+
                                   as.factor(preyear):as.factor(prequarter)| 0 | lsite:year,sample) #as.factor(HHID)+
        summary(results.lm.t.did.a)
        
        results.lm.t.did.s<-felm(logprice ~treatind+X|as.factor(ctr)+
                                   as.factor(year):as.factor(quarter)+
                                   as.factor(preyear):as.factor(prequarter)| 0 | lsite:year,sample) #as.factor(HHID)+
        summary(results.lm.t.did.s)
        
        results.lm.t.did.y<-felm(logprice ~treatind+X|as.factor(lsite)+
                                   as.factor(year):as.factor(quarter)+
                                   as.factor(preyear):as.factor(prequarter)| 0 | lsite:year,sample) #as.factor(HHID)+
        summary(results.lm.t.did.y)
        
        
        #################################
        #Spatial
        
        results.sp.t.did.i<-felm(logprice ~treatind+
                                   Xlag|as.factor(cbl)+
                                   as.factor(year):as.factor(quarter)+
                                   as.factor(preyear):as.factor(prequarter)| 0 | lsite:year,sample) #as.factor(HHID)+
        summary(results.sp.t.did.i)
        
        results.sp.t.did.a<-felm(logprice ~treatind+
                                   Xlag|as.factor(cbg)+
                                   as.factor(year):as.factor(quarter)+
                                   as.factor(preyear):as.factor(prequarter)| 0 | lsite:year,sample) #as.factor(HHID)+
        summary(results.sp.t.did.a)
        
        results.sp.t.did.s<-felm(logprice ~treatind+Xlag|as.factor(ctr)+
                                   as.factor(year):as.factor(quarter)+
                                   as.factor(preyear):as.factor(prequarter)| 0 | lsite:year,sample) #as.factor(HHID)+
        summary(results.sp.t.did.s)
        
        results.sp.t.did.y<-felm(logprice ~treatind+Xlag|as.factor(lsite)+
                                   as.factor(year):as.factor(quarter)+
                                   as.factor(preyear):as.factor(prequarter)| 0 | lsite:year,sample) #as.factor(HHID)+
        summary(results.sp.t.did.y)
        
        #results.lm.t.did.z<-felm(logprice ~treatind+X|as.factor(PropertyZip)+as.factor(year):as.factor(quarter)+
        #                           as.factor(preyear):as.factor(prequarter)| 0 | lsite:year,sample) #as.factor(HHID)+
        #summary(results.lm.t.did.z)
        
        #results.lm.t.did.c<-felm(logprice ~treatind+X|as.factor(PropertyCity)| 0 | lsite:year,sample) #as.factor(HHID)+
        #summary(results.lm.t.did.c)
        
        
        
        Xg<-model.matrix(~ treatmentgroup+#:indx+#bs(timetotreat,5)+bs(day,5)+#as.factor(round(timetotreat,1))+
                           #poly(timetotreat,3)+poly(day,3)+
                           #treatmentgroup*timetotreat+#treatmentgroup*day+
                           LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                           sqfeet+ prediffdate+predate+prelogprice+ indx+#indx+#presstatusd+
                           #as.factor(floor(timetotreat))+
                           as.factor(year):as.factor(quarter)+
                           as.factor(preyear):as.factor(prequarter)-1,sample)
        
        Xglag<-model.matrix(~ treatmentgroup+#:indx+#bs(timetotreat,5)+bs(day,5)+#as.factor(round(timetotreat,1))+
                              #poly(timetotreat,3)+poly(day,3)+
                              #treatmentgroup*timetotreat+#treatmentgroup*day+
                              LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                              sqfeet+ prediffdate+predate+prelogprice+ indx+
                              prelogpricewlag+indx2+#indx+#presstatusd+
                              #as.factor(floor(timetotreat))+
                              as.factor(year):as.factor(quarter)+
                              as.factor(preyear):as.factor(prequarter)-1,sample)
        
        
        sample$tttyear<-as.factor(ceiling(as.numeric(sample$timetotreat)))
        sample1<-subset(sample,tttyear!=1)
        ptlag<-model.matrix(~treatmentgroup:as.factor(ceiling(as.numeric(timetotreat))),sample1)
        #as.factor(cut2(prediffdate, g=quant))+
        
        resid.ptlag<-felm(logprice~ptlag[,c(-10)]+#treatmentgroup+
                            LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                            sqfeet+#timetotreat+#treatmentgroup:timetotreat+#as.factor(tttyear)+
                            prelogprice|as.factor(lsite):as.factor(year):as.factor(quarter)+
                            as.factor(lsite):as.factor(preyear):as.factor(prequarter)+
                            as.factor(ctr)|0|lsite:year,sample1)
        summary(resid.ptlag)
        
        qu<-c("Bottom","Middle","Top")
        yearstest<-9
        qun<-c(-yearstest:-2,2:yearstest)
        allModelFrame <- data.frame(Variable = qun,
                                    Coefficient = as.numeric(coef(summary(resid.ptlag))[,"Estimate"][(11-yearstest):(8+yearstest)]),
                                    SE = as.numeric(coef(summary(resid.ptlag))[,"Cluster s.e."][(11-yearstest):(8+yearstest)]),
                                    modelName = "ptlag")
        
        
        interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier
        leg<-c('10k','8k','6k','4k','2k')
        allModelFrame<-allModelFrame[allModelFrame$modelName!="2k",]
        allModelFrame$modelName<-as.factor( allModelFrame$modelName)
        allModelFrame$modelName<-factor(allModelFrame$modelName,levels(allModelFrame$modelName)[c(2:5,1)])
        
        # Plot
        zp1 <- ggplot(allModelFrame, aes(colour = modelName ))
        zp1 <- zp1 + geom_hline(yintercept = 0, colour = gray(1/2), lty = 2)
        zp1 <- zp1 + geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Coefficient - SE*interval2,
                                         ymax = Coefficient + SE*interval2,color =modelName),
                                     lwd = 1/2, position = position_dodge(width = 1/2),
                                     shape = 21, fill = "WHITE")
        #zp1 <- zp1 + coord_flip() + theme_bw()
        
        #zp1 <- zp1 + geom_line(data = allModelFrame, aes(linetype =modelName ), size = 1) +
        zp1 <- zp1 +theme(legend.position="none")+ ggtitle("Test of Difference in Treatment and Control Groups")+xlab('Year From Treatment')  
        #zp1 <- zp1 + scale_x_continuous(labels=c("1"="B","2"="M","3"="T"))
        zp1 <- zp1 + scale_x_continuous(breaks=qun)
        print(zp1)  # The trick to these is position_dodge().
        
        if(match==""){
          ggsave(file=paste(path,'latex/','pretrends',treatc,dic, 'h5.png', sep=""),height = 5,width =9)
        }
        
        resid.lm.t.did<-felm(logprice ~Xg[,-1]|as.factor(cbl),sample) #as.factor(HHID)+
        summary(resid.lm.t.did)
        
        sample$crdid<-resid.lm.t.did$residuals
        #sample$crdid<-sample$demlogprice-Xg[,-1]%*%coefdid
        #sample$crdid<-sample$crdid-mean(sample$crdid)
        sample$post<-ifelse(sample$timetotreat>0,1,0)
        treatment <- aggregate(crdid ~ treatmentgroup+as.factor(round(timetotreat,2)), data=sample, FUN=mean, na.rm=TRUE)
        names(treatment)[2]<-"timetotreat"
        treatment$timetotreat<-as.numeric(as.character(treatment$timetotreat))
        #treatment<-treatment[as.numeric(treatment$timetotreat)!=0.5,]
        
        ggplot() +
          geom_point(data=subset(treatment,treatmentgroup==1), aes(x=timetotreat, y=crdid, color= "Treatment Group")) +
          geom_point(data=subset(treatment,treatmentgroup==0), aes(x=timetotreat, y=crdid, color= "Control Group")) +
          stat_smooth(method = 'loess', formula = y ~ x  ,data=subset(sample,treatmentgroup==1&post==1),se= TRUE,  aes(x=timetotreat, y=crdid),color= "blue") +
          stat_smooth(method = 'loess', formula = y ~ x  ,data=subset(sample,treatmentgroup==1&post==0),se= TRUE,  aes(x=timetotreat, y=crdid),color= "blue") +
          
          stat_smooth(method = 'loess', formula = y ~ x  ,data=subset(sample,treatmentgroup==0&post==1),se= TRUE,  aes(x=timetotreat, y=crdid),color= "red") +
          stat_smooth(method = 'loess', formula = y ~ x  ,data=subset(sample,treatmentgroup==0&post==0),se= TRUE,  aes(x=timetotreat, y=crdid),color= "red") +
          ggtitle("Common Trends Assumption")+
          xlab('Years from Deletion') + labs(color="Legend") +  geom_vline(xintercept=0)+
          ylab('Monthly Average Residuals')
        
        if(match==""){
          ggsave(file=paste(path,'latex/','lmdidavg',treatc,dic, '.png', sep=""),height = 6,width =10)
        }
        results.lm.t.did.agg<-felm(logprice ~treatst+Xg|
                                     as.factor(PropertyAddressCensusTractAndBlock),sample) #as.factor(HHID)+
        summary(results.lm.t.did.agg)
        
        
        if(match==""){
          betas.lm.did[treat,di]<-as.numeric(coef(summary(results.lm.t.did.agg))[,"Estimate"][1])
          ses.lm.did[treat,di]<-as.numeric(coef(summary(results.lm.t.did.agg))[,"Std. Error"][1])
          ps.lm.did[treat,di]<-as.numeric(coef(summary(results.lm.t.did.agg))[,"Pr(>|t|)"][1])
          
          results.sp.t.did.agg<-felm(logprice ~treatst+Xglag|
                                       as.factor(PropertyAddressCensusTractAndBlock),sample) #as.factor(HHID)+
          summary(results.sp.t.did.agg)
          
          betas.sp.did[treat,di]<-as.numeric(coef(summary(results.sp.t.did.agg))[,"Estimate"][1])
          ses.sp.did[treat,di]<-as.numeric(coef(summary(results.sp.t.did.agg))[,"Std. Error"][1])
          ps.sp.did[treat,di]<-as.numeric(coef(summary(results.sp.t.did.agg))[,"Pr(>|t|)"][1])
          
          
        }
        if(match=="match"){
          betas.match.did[treat,di]<-as.numeric(coef(summary(results.lm.t.did.agg))[,"Estimate"][1])
          ses.match.did[treat,di]<-as.numeric(coef(summary(results.lm.t.did.agg))[,"Std. Error"][1])
          ps.match.did[treat,di]<-as.numeric(coef(summary(results.lm.t.did.agg))[,"Pr(>|t|)"][1])
        }
        if(match==""){
          results.gam<-mgcv::gam(logprice~treatst+Xg+#s(prelogprice,bs="cr")+
                                   #as.factor(year):as.factor(quarter)+
                                   #as.factor(preyear):as.factor(prequarter)+
                                   as.factor(lsite)+s(lat,long,bs="tp",m=3,k=300),data=sample)
          gam.model.t<-mgcv::summary.gam(results.gam)   
          
          betas.gam.did[treat,di]<-as.numeric(gam.model.t$p.table[,"Estimate"])[2]
          ses.gam.did[treat,di]<-as.numeric(gam.model.t$p.table[,"Std. Error"])[2]
          ps.gam.did[treat,di]<-as.numeric(gam.model.t$p.table[,"Pr(>|t|)"])[2]
        }
        
        
        
        #results.spline<-lm(logprice~treatind+X+spTATE,sample)
        #summary(results.spline)
        
        #results.semip<-semip(logprice~treatind+X[,1:8]+preprice,nonpar=~lat+long,window1 = .5, window2 = .5,
        #               kern="tcub",distance="Mahal",targetfull=NULL, print.summary=TRUE, data=sample)
        #summary(results.semip)
        #library(gam)
        library(mgcv)
        if(match==""){
          s=mgcv:::s
          intsp<-model.matrix(~as.factor(sample$lsite):lat+as.factor(sample$lsite):long-1)
          #s1<-sample[lsite!=19,]
          #X1<-X[sample$lsite!=19,]
          #treatind1<-treatind[sample$lsite!=19,]
          results.gam<-mgcv::gam(logprice~treatind+X+#s(prelogprice,bs="cr")+as.factor(lsite):
                                   #s(day,bs="cr")+s(predate,bs="cr")+
                                   as.factor(year):as.factor(quarter)+
                                   as.factor(preyear):as.factor(prequarter)+
                                   as.factor(lsite)+
                                   s(PropertyAddressLatitude,PropertyAddressLongitude,bs="tp",m=3,k=300),data=sample)
          gam.model<-mgcv::summary.gam(results.gam)  
          mgcv::summary.gam(results.gam)  
          #results.gam<-mgcv::gam(logprice~treatind+X+#s(prelogprice,bs="cr")+
          #s(predate,bs="gp")+
          #s(day,bs="gp")+
          #s(lat,long,bs="ts",m=3,k=500),data=sample)
          #mgcv::summary.gam(results.gam) 
          
        }
        if(FALSE){
          
          
          vcov_both_formula <- cluster.vcov(results.lm.t.did.t, ~ lsite + year)
          dim(results.gam$R)
          model.matrix(results.gam$R)
          vcov.HC = solve(t(X)%*%X) %*% t(X)%*%diag(ehat^2)%*%X %*% solve(t(X)%*%X)
          mg$Vp <- vcov.HC
          summary(mg)
          all.equal(as.numeric(predict(mg,se.fit=TRUE)$se.fit),se.yhat.HC)
          
          results.gam<-mgcv::gam(logprice~treatind+X+#s(prelogprice,bs="cr")+
                                   #s(day,bs="gp")+
                                   s(lsite,bs="re")+s(year,bs="re")+
                                   s(lat,long,bs="tp",m=3,k=300),data=sample)
          mgcv::summary.gam(results.gam) 
          
          
          results.gam<-mgcv::gamm(logprice~treatind+X+#s(prelogprice,bs="cr")+
                                    #s(day,bs="gp")+
                                    s(lat,long,bs="tp",m=3,k=200),
                                  correlation=corSymm(form=~1|PropertyAddressCensusTractAndBlock),data=sample)
          mgcv::summary.gam(results.gam) 
          
          qgam.fit<-qgam(logprice~treatst+nX+#s(prelogprice,bs="cr")+
                           #s(day,bs="gp")+
                           s(lat,long,bs="tp",m=3,k=300), lsig = -1,
                         data=sample,qu=.5, err = 0.05,control = list("tol" = 0.01))
          summary.gam(qgam.fit)
          
          
          
          neX<-model.matrix(~ treatmentgroup+#indx:treatmentgroup:as.factor(lsite)+#:indx+#bs(timetotreat,5)+bs(day,5)+#as.factor(round(timetotreat,1))+
                              #poly(timetotreat,3)+poly(day,3)+
                              #treatmentgroup*timetotreat+#treatmentgroup*day+
                              LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                              sqfeet+as.factor(year)+as.factor(preyear)+
                              as.factor(quarter)+as.factor(prequarter)-1,sample)
          
          nX<-neX
          qr.nX <- qr(nX, tol=1e-2, LAPACK = FALSE)
          (rnknX <- qr.nX$rank)  ## 4 (number of non-collinear columns)
          (keepnX <- qr.nX$pivot[seq_len(rnknX)])
          ## 1 2 4 5 
          nX <- nX[,keepnX]
          
          rq.o<-rq.fit.sfn(as.matrix.csr(nX),y,
                           tmpmax=floor(10000+exp(-12.1)*(dim(nX)[1]*20-1)^2.35))
          
          fit.qr<-rq(logprice~nX,tau=.5, data=sample,method="sfn",na.action = na.omit)
          fit.qr
          summary(fit.qr,se = "boot")
          
          plot(indx,sample$logprice,xlab="Treatind", ylab="Log Price")
          taus <- c(.1,.9)
          abline(rq(logprice~X2,tau=.5,data=sample),col="blue")
          abline(lm(logprice~X2,data=sample),lty = 3,col="red")
          for( i in 1:length(taus)){
            abline(rq(logprice~X2,tau=taus[i],data=sample),col="gray")
          }
          
          
          
          
          #PropertyAddressLatitude PropertyAddressLongitude
          
          
          
          Wst<-readRDS(paste0(path,"Wmat",treatc,dic,".rds"), refhook = NULL)
          summary(rowSums(Wst))
          
          if(mean(rowSums(Wst))==1){
            Wst<-mat2listw(Wst, style="W")
            
            
            #finalb.lag.2sls.robust2 <- gstslshet(logprice~treatind+X,Wst, data = sample,
            #                                    initial.value = 0.2, eps =1e-2, inverse=FALSE,sarar=FALSE)
            #summary(finalb.lag.2sls.robust2)
            #effects.finalb.lag.2sls.robust2<- impacts(finalb.lag.2sls.robust2, listw= Wst, R=100)
            #summary(effects.finalb.lag.2sls.robust2, zstats=TRUE, short=TRUE)
            
            results.stsls<-sacsarlm(logprice~treatind+X, data = sample, listw=Wst, zero.policy = NULL,
                                    na.action = na.fail)
            summary(results.stsls)
          }
          # effects.finalb.lag.2sls.robust2<- impacts(results.stsls, listw= Wst, R=100)
          #summary(effects.finalb.lag.2sls.robust2, zstats=TRUE, short=TRUE)
        }
        if(match==""){
          if(treatc=="TATE"){
            betas.lm.t.did.TATE.i[,di]<-as.numeric(coef(summary(results.lm.t.did.i))[,"Estimate"][1:quant])
            betas.lm.t.did.TATE.a[,di]<-as.numeric(coef(summary(results.lm.t.did.a))[,"Estimate"][1:quant])
            betas.lm.t.did.TATE.s[,di]<-as.numeric(coef(summary(results.lm.t.did.s))[,"Estimate"][1:quant])
            betas.lm.t.did.TATE.y[,di]<-as.numeric(coef(summary(results.lm.t.did.y))[,"Estimate"][1:quant])
            
            betas.sp.t.did.TATE.i[,di]<-as.numeric(coef(summary(results.sp.t.did.i))[,"Estimate"][1:quant])
            betas.sp.t.did.TATE.a[,di]<-as.numeric(coef(summary(results.sp.t.did.a))[,"Estimate"][1:quant])
            betas.sp.t.did.TATE.s[,di]<-as.numeric(coef(summary(results.sp.t.did.s))[,"Estimate"][1:quant])
            betas.sp.t.did.TATE.y[,di]<-as.numeric(coef(summary(results.sp.t.did.y))[,"Estimate"][1:quant])
            
            betas.gam.t.did.TATE.[,di]<-as.numeric(gam.model$p.table[,"Estimate"])[2:(quant+1)]
            #betas.lm.t.es.TATE[,di]<-as.numeric(coef(summary(results.lm.t.es))[,"Estimate"][1:quant])
            ses.lm.t.did.TATE.i[,di]<-as.numeric(coef(summary(results.lm.t.did.i))[,"Cluster s.e."][1:quant])
            ses.lm.t.did.TATE.a[,di]<-as.numeric(coef(summary(results.lm.t.did.a))[,"Cluster s.e."][1:quant])
            ses.lm.t.did.TATE.s[,di]<-as.numeric(coef(summary(results.lm.t.did.s))[,"Cluster s.e."][1:quant])
            ses.lm.t.did.TATE.y[,di]<-as.numeric(coef(summary(results.lm.t.did.y))[,"Cluster s.e."][1:quant])
            
            ses.sp.t.did.TATE.i[,di]<-as.numeric(coef(summary(results.sp.t.did.i))[,"Cluster s.e."][1:quant])
            ses.sp.t.did.TATE.a[,di]<-as.numeric(coef(summary(results.sp.t.did.a))[,"Cluster s.e."][1:quant])
            ses.sp.t.did.TATE.s[,di]<-as.numeric(coef(summary(results.sp.t.did.s))[,"Cluster s.e."][1:quant])
            ses.sp.t.did.TATE.y[,di]<-as.numeric(coef(summary(results.sp.t.did.y))[,"Cluster s.e."][1:quant])
            
            ses.gam.t.did.TATE.[,di]<-as.numeric(gam.model$p.table[,"Std. Error"])[2:(quant+1)]
            #ses.lm.t.es.TATE[,di]<-as.numeric(coef(summary(results.lm.t.es))[,"Std. Error"][1:quant])
            ps.lm.t.did.TATE.i[,di]<-as.numeric(coef(summary(results.lm.t.did.i))[,"Pr(>|t|)"][1:quant])
            ps.lm.t.did.TATE.a[,di]<-as.numeric(coef(summary(results.lm.t.did.a))[,"Pr(>|t|)"][1:quant])
            ps.lm.t.did.TATE.s[,di]<-as.numeric(coef(summary(results.lm.t.did.s))[,"Pr(>|t|)"][1:quant])
            ps.lm.t.did.TATE.y[,di]<-as.numeric(coef(summary(results.lm.t.did.y))[,"Pr(>|t|)"][1:quant])
            
            ps.sp.t.did.TATE.i[,di]<-as.numeric(coef(summary(results.sp.t.did.i))[,"Pr(>|t|)"][1:quant])
            ps.sp.t.did.TATE.a[,di]<-as.numeric(coef(summary(results.sp.t.did.a))[,"Pr(>|t|)"][1:quant])
            ps.sp.t.did.TATE.s[,di]<-as.numeric(coef(summary(results.sp.t.did.s))[,"Pr(>|t|)"][1:quant])
            ps.sp.t.did.TATE.y[,di]<-as.numeric(coef(summary(results.sp.t.did.y))[,"Pr(>|t|)"][1:quant])
            
            ps.gam.t.did.TATE.[,di]<-as.numeric(gam.model$p.table[,"Pr(>|t|)"])[2:(quant+1)]
            #ps.lm.t.es.TATE[,di]<-as.numeric(coef(summary(results.lm.t.es))[,"Pr(>|t|)"][1:quant])
          }
          
          if(treatc=="MUATE"){
            betas.lm.t.did.MUATE.i[,di]<-as.numeric(coef(summary(results.lm.t.did.i))[,"Estimate"][1:quant])
            betas.lm.t.did.MUATE.a[,di]<-as.numeric(coef(summary(results.lm.t.did.a))[,"Estimate"][1:quant])
            betas.lm.t.did.MUATE.s[,di]<-as.numeric(coef(summary(results.lm.t.did.s))[,"Estimate"][1:quant])
            betas.lm.t.did.MUATE.y[,di]<-as.numeric(coef(summary(results.lm.t.did.y))[,"Estimate"][1:quant])
            
            betas.sp.t.did.MUATE.i[,di]<-as.numeric(coef(summary(results.sp.t.did.i))[,"Estimate"][1:quant])
            betas.sp.t.did.MUATE.a[,di]<-as.numeric(coef(summary(results.sp.t.did.a))[,"Estimate"][1:quant])
            betas.sp.t.did.MUATE.s[,di]<-as.numeric(coef(summary(results.sp.t.did.s))[,"Estimate"][1:quant])
            betas.sp.t.did.MUATE.y[,di]<-as.numeric(coef(summary(results.sp.t.did.y))[,"Estimate"][1:quant])
            
            betas.gam.t.did.MUATE.[,di]<-as.numeric(gam.model$p.table[,"Estimate"])[2:(quant+1)]
            #betas.lm.t.es.TATE[,di]<-as.numeric(coef(summary(results.lm.t.es))[,"Estimate"][1:quant])
            ses.lm.t.did.MUATE.i[,di]<-as.numeric(coef(summary(results.lm.t.did.i))[,"Cluster s.e."][1:quant])
            ses.lm.t.did.MUATE.a[,di]<-as.numeric(coef(summary(results.lm.t.did.a))[,"Cluster s.e."][1:quant])
            ses.lm.t.did.MUATE.s[,di]<-as.numeric(coef(summary(results.lm.t.did.s))[,"Cluster s.e."][1:quant])
            ses.lm.t.did.MUATE.y[,di]<-as.numeric(coef(summary(results.lm.t.did.y))[,"Cluster s.e."][1:quant])
            
            ses.sp.t.did.MUATE.i[,di]<-as.numeric(coef(summary(results.sp.t.did.i))[,"Cluster s.e."][1:quant])
            ses.sp.t.did.MUATE.a[,di]<-as.numeric(coef(summary(results.sp.t.did.a))[,"Cluster s.e."][1:quant])
            ses.sp.t.did.MUATE.s[,di]<-as.numeric(coef(summary(results.sp.t.did.s))[,"Cluster s.e."][1:quant])
            ses.sp.t.did.MUATE.y[,di]<-as.numeric(coef(summary(results.sp.t.did.y))[,"Cluster s.e."][1:quant])
            
            
            ses.gam.t.did.MUATE.[,di]<-as.numeric(gam.model$p.table[,"Std. Error"])[2:(quant+1)]
            #ses.lm.t.es.TATE[,di]<-as.numeric(coef(summary(results.lm.t.es))[,"Std. Error"][1:quant])
            ps.lm.t.did.MUATE.i[,di]<-as.numeric(coef(summary(results.lm.t.did.i))[,"Pr(>|t|)"][1:quant])
            ps.lm.t.did.MUATE.a[,di]<-as.numeric(coef(summary(results.lm.t.did.a))[,"Pr(>|t|)"][1:quant])
            ps.lm.t.did.MUATE.s[,di]<-as.numeric(coef(summary(results.lm.t.did.s))[,"Pr(>|t|)"][1:quant])
            ps.lm.t.did.MUATE.y[,di]<-as.numeric(coef(summary(results.lm.t.did.y))[,"Pr(>|t|)"][1:quant])
            
            ps.sp.t.did.MUATE.i[,di]<-as.numeric(coef(summary(results.sp.t.did.i))[,"Pr(>|t|)"][1:quant])
            ps.sp.t.did.MUATE.a[,di]<-as.numeric(coef(summary(results.sp.t.did.a))[,"Pr(>|t|)"][1:quant])
            ps.sp.t.did.MUATE.s[,di]<-as.numeric(coef(summary(results.sp.t.did.s))[,"Pr(>|t|)"][1:quant])
            ps.sp.t.did.MUATE.y[,di]<-as.numeric(coef(summary(results.sp.t.did.y))[,"Pr(>|t|)"][1:quant])
            
            ps.gam.t.did.MUATE.[,di]<-as.numeric(gam.model$p.table[,"Pr(>|t|)"])[2:(quant+1)]
            #ps.lm.t.es.TATE[,di]<-as.numeric(coef(summary(results.lm.t.es))[,"Pr(>|t|)"][1:quant])
          }
          
          if(treatc=="WLATE"){
            betas.lm.t.did.WLATE.i[,di]<-as.numeric(coef(summary(results.lm.t.did.i))[,"Estimate"][1:quant])
            betas.lm.t.did.WLATE.a[,di]<-as.numeric(coef(summary(results.lm.t.did.a))[,"Estimate"][1:quant])
            betas.lm.t.did.WLATE.s[,di]<-as.numeric(coef(summary(results.lm.t.did.s))[,"Estimate"][1:quant])
            betas.lm.t.did.WLATE.y[,di]<-as.numeric(coef(summary(results.lm.t.did.y))[,"Estimate"][1:quant])
            
            betas.sp.t.did.WLATE.i[,di]<-as.numeric(coef(summary(results.sp.t.did.i))[,"Estimate"][1:quant])
            betas.sp.t.did.WLATE.a[,di]<-as.numeric(coef(summary(results.sp.t.did.a))[,"Estimate"][1:quant])
            betas.sp.t.did.WLATE.s[,di]<-as.numeric(coef(summary(results.sp.t.did.s))[,"Estimate"][1:quant])
            betas.sp.t.did.WLATE.y[,di]<-as.numeric(coef(summary(results.sp.t.did.y))[,"Estimate"][1:quant])
            
            betas.gam.t.did.WLATE.[,di]<-as.numeric(gam.model$p.table[,"Estimate"])[2:(quant+1)]
            #betas.lm.t.es.TATE[,di]<-as.numeric(coef(summary(results.lm.t.es))[,"Estimate"][1:quant])
            ses.lm.t.did.WLATE.i[,di]<-as.numeric(coef(summary(results.lm.t.did.i))[,"Cluster s.e."][1:quant])
            ses.lm.t.did.WLATE.a[,di]<-as.numeric(coef(summary(results.lm.t.did.a))[,"Cluster s.e."][1:quant])
            ses.lm.t.did.WLATE.s[,di]<-as.numeric(coef(summary(results.lm.t.did.s))[,"Cluster s.e."][1:quant])
            ses.lm.t.did.WLATE.y[,di]<-as.numeric(coef(summary(results.lm.t.did.y))[,"Cluster s.e."][1:quant])
            
            ses.sp.t.did.WLATE.i[,di]<-as.numeric(coef(summary(results.sp.t.did.i))[,"Cluster s.e."][1:quant])
            ses.sp.t.did.WLATE.a[,di]<-as.numeric(coef(summary(results.sp.t.did.a))[,"Cluster s.e."][1:quant])
            ses.sp.t.did.WLATE.s[,di]<-as.numeric(coef(summary(results.sp.t.did.s))[,"Cluster s.e."][1:quant])
            ses.sp.t.did.WLATE.y[,di]<-as.numeric(coef(summary(results.sp.t.did.y))[,"Cluster s.e."][1:quant])
            
            ses.gam.t.did.WLATE.[,di]<-as.numeric(gam.model$p.table[,"Std. Error"])[2:(quant+1)]
            #ses.lm.t.es.TATE[,di]<-as.numeric(coef(summary(results.lm.t.es))[,"Std. Error"][1:quant])
            ps.lm.t.did.WLATE.i[,di]<-as.numeric(coef(summary(results.lm.t.did.i))[,"Pr(>|t|)"][1:quant])
            ps.lm.t.did.WLATE.a[,di]<-as.numeric(coef(summary(results.lm.t.did.a))[,"Pr(>|t|)"][1:quant])
            ps.lm.t.did.WLATE.s[,di]<-as.numeric(coef(summary(results.lm.t.did.s))[,"Pr(>|t|)"][1:quant])
            ps.lm.t.did.WLATE.y[,di]<-as.numeric(coef(summary(results.lm.t.did.y))[,"Pr(>|t|)"][1:quant])
            
            ps.sp.t.did.WLATE.i[,di]<-as.numeric(coef(summary(results.sp.t.did.i))[,"Pr(>|t|)"][1:quant])
            ps.sp.t.did.WLATE.a[,di]<-as.numeric(coef(summary(results.sp.t.did.a))[,"Pr(>|t|)"][1:quant])
            ps.sp.t.did.WLATE.s[,di]<-as.numeric(coef(summary(results.sp.t.did.s))[,"Pr(>|t|)"][1:quant])
            ps.sp.t.did.WLATE.y[,di]<-as.numeric(coef(summary(results.sp.t.did.y))[,"Pr(>|t|)"][1:quant])
            
            ps.gam.t.did.WLATE.[,di]<-as.numeric(gam.model$p.table[,"Pr(>|t|)"])[2:(quant+1)]
            #ps.lm.t.es.TATE[,di]<-as.numeric(coef(summary(results.lm.t.es))[,"Pr(>|t|)"][1:quant])
          }
        }
        if(match=="match"){
          if(treatc=="TATE"){
            betas.match.t.did.TATE.i[,di]<-as.numeric(coef(summary(results.lm.t.did.i))[,"Estimate"][1:quant])
            betas.match.t.did.TATE.a[,di]<-as.numeric(coef(summary(results.lm.t.did.a))[,"Estimate"][1:quant])
            betas.match.t.did.TATE.s[,di]<-as.numeric(coef(summary(results.lm.t.did.s))[,"Estimate"][1:quant])
            betas.match.t.did.TATE.y[,di]<-as.numeric(coef(summary(results.lm.t.did.y))[,"Estimate"][1:quant])
            
            #betas.gam.t.did.TATE.[,di]<-as.numeric(gam.model$p.table[,"Estimate"])[2:(quant+1)]
            #betas.lm.t.es.TATE[,di]<-as.numeric(coef(summary(results.lm.t.es))[,"Estimate"][1:quant])
            ses.match.t.did.TATE.i[,di]<-as.numeric(coef(summary(results.lm.t.did.i))[,"Cluster s.e."][1:quant])
            ses.match.t.did.TATE.a[,di]<-as.numeric(coef(summary(results.lm.t.did.a))[,"Cluster s.e."][1:quant])
            ses.match.t.did.TATE.s[,di]<-as.numeric(coef(summary(results.lm.t.did.s))[,"Cluster s.e."][1:quant])
            ses.match.t.did.TATE.y[,di]<-as.numeric(coef(summary(results.lm.t.did.y))[,"Cluster s.e."][1:quant])
            
            #ses.gam.t.did.TATE.[,di]<-as.numeric(gam.model$p.table[,"Std. Error"])[2:(quant+1)]
            #ses.lm.t.es.TATE[,di]<-as.numeric(coef(summary(results.lm.t.es))[,"Std. Error"][1:quant])
            ps.match.t.did.TATE.i[,di]<-as.numeric(coef(summary(results.lm.t.did.i))[,"Pr(>|t|)"][1:quant])
            ps.match.t.did.TATE.a[,di]<-as.numeric(coef(summary(results.lm.t.did.a))[,"Pr(>|t|)"][1:quant])
            ps.match.t.did.TATE.s[,di]<-as.numeric(coef(summary(results.lm.t.did.s))[,"Pr(>|t|)"][1:quant])
            ps.match.t.did.TATE.y[,di]<-as.numeric(coef(summary(results.lm.t.did.y))[,"Pr(>|t|)"][1:quant])
            
            #ps.gam.t.did.TATE.[,di]<-as.numeric(gam.model$p.table[,"Pr(>|t|)"])[2:(quant+1)]
            #ps.lm.t.es.TATE[,di]<-as.numeric(coef(summary(results.lm.t.es))[,"Pr(>|t|)"][1:quant])
          }
          
          if(treatc=="MUATE"){
            betas.match.t.did.MUATE.i[,di]<-as.numeric(coef(summary(results.lm.t.did.i))[,"Estimate"][1:quant])
            betas.match.t.did.MUATE.a[,di]<-as.numeric(coef(summary(results.lm.t.did.a))[,"Estimate"][1:quant])
            betas.match.t.did.MUATE.s[,di]<-as.numeric(coef(summary(results.lm.t.did.s))[,"Estimate"][1:quant])
            betas.match.t.did.MUATE.y[,di]<-as.numeric(coef(summary(results.lm.t.did.y))[,"Estimate"][1:quant])
            
            #betas.gam.t.did.MUATE.[,di]<-as.numeric(gam.model$p.table[,"Estimate"])[2:(quant+1)]
            #betas.lm.t.es.TATE[,di]<-as.numeric(coef(summary(results.lm.t.es))[,"Estimate"][1:quant])
            ses.match.t.did.MUATE.i[,di]<-as.numeric(coef(summary(results.lm.t.did.i))[,"Cluster s.e."][1:quant])
            ses.match.t.did.MUATE.a[,di]<-as.numeric(coef(summary(results.lm.t.did.a))[,"Cluster s.e."][1:quant])
            ses.match.t.did.MUATE.s[,di]<-as.numeric(coef(summary(results.lm.t.did.s))[,"Cluster s.e."][1:quant])
            ses.match.t.did.MUATE.y[,di]<-as.numeric(coef(summary(results.lm.t.did.y))[,"Cluster s.e."][1:quant])
            
            #ses.gam.t.did.MUATE.[,di]<-as.numeric(gam.model$p.table[,"Std. Error"])[2:(quant+1)]
            #ses.lm.t.es.TATE[,di]<-as.numeric(coef(summary(results.lm.t.es))[,"Std. Error"][1:quant])
            ps.match.t.did.MUATE.i[,di]<-as.numeric(coef(summary(results.lm.t.did.i))[,"Pr(>|t|)"][1:quant])
            ps.match.t.did.MUATE.a[,di]<-as.numeric(coef(summary(results.lm.t.did.a))[,"Pr(>|t|)"][1:quant])
            ps.match.t.did.MUATE.s[,di]<-as.numeric(coef(summary(results.lm.t.did.s))[,"Pr(>|t|)"][1:quant])
            ps.match.t.did.MUATE.y[,di]<-as.numeric(coef(summary(results.lm.t.did.y))[,"Pr(>|t|)"][1:quant])
            
            #ps.gam.t.did.MUATE.[,di]<-as.numeric(gam.model$p.table[,"Pr(>|t|)"])[2:(quant+1)]
            #ps.lm.t.es.TATE[,di]<-as.numeric(coef(summary(results.lm.t.es))[,"Pr(>|t|)"][1:quant])
          }
          
          if(treatc=="WLATE"){
            betas.match.t.did.WLATE.i[,di]<-as.numeric(coef(summary(results.lm.t.did.i))[,"Estimate"][1:quant])
            betas.match.t.did.WLATE.a[,di]<-as.numeric(coef(summary(results.lm.t.did.a))[,"Estimate"][1:quant])
            betas.match.t.did.WLATE.s[,di]<-as.numeric(coef(summary(results.lm.t.did.s))[,"Estimate"][1:quant])
            betas.match.t.did.WLATE.y[,di]<-as.numeric(coef(summary(results.lm.t.did.y))[,"Estimate"][1:quant])
            
            #betas.gam.t.did.WLATE.[,di]<-as.numeric(gam.model$p.table[,"Estimate"])[2:(quant+1)]
            #betas.lm.t.es.TATE[,di]<-as.numeric(coef(summary(results.lm.t.es))[,"Estimate"][1:quant])
            ses.match.t.did.WLATE.i[,di]<-as.numeric(coef(summary(results.lm.t.did.i))[,"Cluster s.e."][1:quant])
            ses.match.t.did.WLATE.a[,di]<-as.numeric(coef(summary(results.lm.t.did.a))[,"Cluster s.e."][1:quant])
            ses.match.t.did.WLATE.s[,di]<-as.numeric(coef(summary(results.lm.t.did.s))[,"Cluster s.e."][1:quant])
            ses.match.t.did.WLATE.y[,di]<-as.numeric(coef(summary(results.lm.t.did.y))[,"Cluster s.e."][1:quant])
            
            #ses.gam.t.did.WLATE.[,di]<-as.numeric(gam.model$p.table[,"Std. Error"])[2:(quant+1)]
            #ses.lm.t.es.TATE[,di]<-as.numeric(coef(summary(results.lm.t.es))[,"Std. Error"][1:quant])
            ps.match.t.did.WLATE.i[,di]<-as.numeric(coef(summary(results.lm.t.did.i))[,"Pr(>|t|)"][1:quant])
            ps.match.t.did.WLATE.a[,di]<-as.numeric(coef(summary(results.lm.t.did.a))[,"Pr(>|t|)"][1:quant])
            ps.match.t.did.WLATE.s[,di]<-as.numeric(coef(summary(results.lm.t.did.s))[,"Pr(>|t|)"][1:quant])
            ps.match.t.did.WLATE.y[,di]<-as.numeric(coef(summary(results.lm.t.did.y))[,"Pr(>|t|)"][1:quant])
            
            #ps.gam.t.did.WLATE.[,di]<-as.numeric(gam.model$p.table[,"Pr(>|t|)"])[2:(quant+1)]
            #ps.lm.t.es.TATE[,di]<-as.numeric(coef(summary(results.lm.t.es))[,"Pr(>|t|)"][1:quant])
          }
        }
        if(FALSE){
          sdf<-10
          lat<-sample$PropertyAddressLatitude
          long<-sample$PropertyAddressLongitude
          splat<-bs(lat, df = sdf)
          splong<-bs(long, df = sdf)
          spint<-model.matrix(~splat:splong)
          
          xcTATE<-cbind(splat,splong,spint,lat,long,poly(sample$day,5),bs(sample$day, df = 10))
          year<-dplyr::select(sample, starts_with('year'))
          
          feTATE<-model.matrix(~ treatind+indx+
                                 #prelogprice+prediffdate+predate+
                                 data.matrix(year[,4:25]),sample)#+#timefedControlsComplete+timefed +
          #aftfinalnpl+timefinalnplfe+
          #data.matrix(exdum)+
          #data.matrix(year[,25]),sample)
          feTATE<-as.matrix(feTATE[,SD(feTATE)>0])
          feTATE<-as.matrix(feTATE[,!duplicated(cor(feTATE))])
          qr.X <- qr(feTATE, tol=1e-3, LAPACK = FALSE)
          (rnkX <- qr.X$rank)  ## 4 (number of non-collinear columns)
          (keep <- qr.X$pivot[seq_len(rnkX)])
          ## 1 2 4 5 
          feTATE <- feTATE[,keep]
          
          W<-cbind(sample$day,lat,long,X)
          colnames(W)[1]<-"day"
          colnames(W)[2]<-"lat"
          colnames(W)[3]<-"long"
          A<-sample$treatmentgroup
          V<-X
          Time<-sample$day
          
          results.tmle.t.did <- tmleMSM(Y = sample$logprice, A = A, W = W, V = V, #T= Time,
                                        MSM = "A + V",family="gaussian", 
                                        Q.SL.library = SL.library2$as.list(),
                                        g.SL.library = SL.library2$as.list(),
                                        #Qform = Y ~ A+V+W,
                                        #gform = A~1,
                                        #hAVform = A~ 1,
                                        ub = 20,
                                        V_SL =5,
                                        alpha = 0.90,
                                        inference = TRUE,
                                        verbose=TRUE)
          print(results.tmle.t.did)
          summary(results.tmle.t.did)
          
          
          sdf<-5
          lat<-sample$PropertyAddressLatitude
          long<-sample$PropertyAddressLongitude
          splat<-bs(lat, df = sdf)
          splong<-bs(long, df = sdf)
          spint<-model.matrix(~splat:splong)
          
          xcTATE<-cbind(splat,splong,spint,lat,long,poly(sample$day,5),bs(sample$day, df = 10))
          feTATEes<-model.matrix(~ treatexCC+#aftfinalnpl+
                                   #as.matrix(exdum)
                                   -1,tsample)
          feTATEes<-as.matrix(feTATEes[,SD(feTATEes)>0])
          if(dim(feTATEes)[2]>1){
            feTATEes<-feTATEes[,!duplicated(cor(feTATEes))]
          }
          qr.X <- qr(feTATEes, tol=1e-3, LAPACK = FALSE)
          (rnkX <- qr.X$rank)  ## 4 (number of non-collinear columns)
          (keep <- qr.X$pivot[seq_len(rnkX)])
          ## 1 2 4 5 
          feTATEes <- feTATEes[,keep]
          if(nocc==0&dim(as.matrix(feTATEes))[2]>0){
            if(dim(data.matrix(tsample))[1]>40){
              W<-cbind(xcTATE[sample$control==0,],xTATE[sample$control==0,],feTATEes)
              
              A<-sample[control==0,treatst]
              V<-feTATEes
              Time<-sample[control==0,day]
              
              results.tmle.t.es <- tmleMSM(Y = tsample$logprice, A = A, W = W, V = V, #T= Time,
                                           MSM = "A + V",family="gaussian", 
                                           Q.SL.library = SL.library2$as.list(),
                                           g.SL.library = PS.library2$as.list(),
                                           #Qform = Y ~ A+V+W,
                                           #gform = A~1,
                                           #hAVform = A~ 1,
                                           ub = 20,
                                           V_SL =5,
                                           alpha = 0.90,
                                           inference = TRUE,
                                           verbose=TRUE)
              print(results.tmle.t.es)
              
              betas.tmle.t.es[di,treat]<-results.tmle.t.es$psi["A"]
              ses.tmle.t.es[di,treat]<-results.tmle.t.es$se["A"]
              ps.tmle.t.es[di,treat]<-results.tmle.t.es$pvalue["A"]
              
              cc.betas.tmle.t.es[di,treat]<-results.tmle.t.es$psi["V"]
              cc.ses.tmle.t.es[di,treat]<-results.tmle.t.es$se["V"]
              cc.ps.tmle.t.es[di,treat]<-results.tmle.t.es$pvalue["V"]
              
              
            }}
          
          
          
          betas.tmle.t.did[di,treat]<-results.tmle.t.did$psi["A"]
          ses.tmle.t.did[di,treat]<-results.tmle.t.did$se["A"]
          ps.tmle.t.did[di,treat]<-results.tmle.t.did$pvalue["A"]
          
          if(nocc==0 &notg==0){
            cc.betas.tmle.t.did[di,treat]<-results.tmle.t.did$psi["VtreatexCC"]
            cc.ses.tmle.t.did[di,treat]<-results.tmle.t.did$se["VtreatexCC"]
            cc.ps.tmle.t.did[di,treat]<-results.tmle.t.did$pvalue["VtreatexCC"]
          }
        }
        
        print(paste0('distance = ',dic))
        print(paste0('treat = ',treatc))
        
      }
    }
  }
}

for(statchange in c('')){
  for(meth in c('lm','gam','match','sp')){
    for(inf in c('did')){
      for(treat in  treatl){
        for(spec in specl){
          #treat<-"TATE"
          
          p<-get(paste0(statchange,'ps.',meth,'.t.',inf,'.',treat,'.',spec))
          mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ", ifelse(p < .1, "^\\bullet  ", " "))))
          if(!is.na(p[1,1])){
            #pb<-exp(get(paste0(statchange,'betas.',meth,'.t.',inf,'.',treat)))-1
            #rpb<-round(pb,3)
            #se<-round(exp(get(paste0(statchange,'ses.',meth,'.t.',inf,'.',treat)))-1,3)
            
            pb<-get(paste0(statchange,'betas.',meth,'.t.',inf,'.',treat,'.',spec))
            rpb<-round(pb,3)
            se<-round(get(paste0(statchange,'ses.',meth,'.t.',inf,'.',treat,'.',spec)),3)
            
            srpb <- matrix(paste(rpb, mystars, sep=""), ncol=dim(pb)[2] )
            nsrpb<-rbind(c("",laglead),cbind(dist,srpb))
            
            #colnames(srpb)<-laglead
            #rownames(srpb)<-dist
            
            results.mat<-matrix(nrow= 2*dim(srpb)[1],ncol= dim(srpb)[2])
            
            for(i in 1:dim(results.mat)[1]){
              if(i %% 2 != 0){
                results.mat[i,]<-srpb[ceiling(i/2),]
                #    rownames(ols.mat)[i]<-rownames(srpb)[ceiling(i/2)]
              }
              if(i %% 2 == 0){
                results.mat[i,]<-paste0('(',se[ceiling(i/2),],')')
              }
            }
            
            results.mat<-rbind(c('10k','8k','6k','4k','2k'),results.mat)
            if(meth=='lm'){
              results.mat<-rbind(c('OLS','','','',''),results.mat)
            }
            if(meth=='gam'){
              results.mat<-rbind(c('GAM','','','',''),results.mat)
            }
            if(meth=='match'){
              results.mat<-rbind(c('Matching','','','',''),results.mat)
            }
            if(meth=='sp'){
              results.mat<-rbind(c('Spatial Lag','','','',''),results.mat)
            }
            # rn<-c(paste0("(",qcut[1],","),paste0(qcut[2],"]"),paste0("(",qcut[2],","),paste0(qcut[3],"]"),
            #      paste0("(",qcut[3],","),paste0(qcut[4],"]"),
            #     paste0("(",qcut[4],","),paste0(qcut[5],"]"),paste0("(",qcut[5],","),
            #    paste0(qcut[6],"]"),paste0("(",qcut[6],","),paste0(qcut[7],"]"),paste0("(",qcut[7],","),
            #   paste0(qcut[8],"]"),paste0("(",qcut[8],","),paste0(qcut[9],"]"),paste0("(",qcut[9],","),
            #  paste0(qcut[10],"]"),
            # paste0("(",qcut[10],","),paste0(qcut[11],"]"))
            rn<-c(paste0("(",qcut[1],",",qcut[2],"]"),"Bottom Ten Percentile" ,
                  paste0("(",qcut[2],",",qcut[quant],"]"),"Middle 80 Percentile" ,
                  paste0("(",qcut[quant],",",qcut[quant+1],"]"), "Top Ten Percentile")
            rn2<-c("                                      ","                      ",paste0("(",qcut[1],", ",qcut[2],"]")," ",paste0("(",qcut[2],",  ",qcut[3],"]"),"  ",
                   paste0("(",qcut[3],", ",qcut[4],"]"),"   ",
                   paste0("(",qcut[4],",",qcut[5],"]"),"     ",paste0("(",qcut[5],",",qcut[6],"]"),"      ",
                   paste0("(",qcut[6],",",qcut[7],"]"),"       ",paste0("(",qcut[7],",",qcut[8],"]"),"        ",
                   paste0("(",qcut[8],",",qcut[9],"]"),"           ", paste0("(",qcut[9],",",qcut[10],"]"),"           ",
                   paste0("(",qcut[10],",",qcut[11],"]"),"                       ")
            #rownames(results.mat)<-rn2
            results.mat<-cbind(rn2,results.mat)
            xtab<-xtable(results.mat)
            align(xtab) <- "rl|rrrrr"
            print.xtable(xtab,include.rownames=FALSE, hline.after = c(0,1,2,dim(results.mat)[1]),
                         include.colnames=FALSE, sanitize.text.function = identity,
                         #caption = "example", 
                         label = paste0("tab:",meth,inf,statchange,treat,spec),
                         type="latex", file=paste0(path,'latex/',meth,inf,statchange,treat,spec,".tex"))
            
            
            
            
            qu<-c("Bottom","Middle","Top")
            qun<-c(1:10)
            allModelFrame <- data.frame(Variable = qun,
                                        Coefficient = pb[,1],
                                        SE = se[, 1],
                                        modelName = "10k")
            for(i in 2:length(dist)){
              di<-dist[i]
              modelFrame <- data.frame(Variable =  qun,
                                       Coefficient = pb[,i],
                                       SE = se[, i],
                                       modelName = di)
              allModelFrame <- data.frame(rbind(allModelFrame,modelFrame))
              
            }
            
            interval2 <- -qnorm((1-0.90)/2)  # 95% multiplier
            leg<-c('10k','8k','6k','4k','2k')
            allModelFrame<-allModelFrame[allModelFrame$modelName!="2k",]
            allModelFrame$modelName<-as.factor( allModelFrame$modelName)
            allModelFrame$modelName<-factor(allModelFrame$modelName,levels(allModelFrame$modelName)[c(2:5,1)])
            
            # Plot
            zp1 <- ggplot(allModelFrame, aes(colour = modelName ))
            zp1 <- zp1 + geom_hline(yintercept = 0, colour = gray(1/2), lty = 2)
            zp1 <- zp1 + geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Coefficient - SE*interval2,
                                             ymax = Coefficient + SE*interval2,color =modelName),
                                         lwd = 1/2, position = position_dodge(width = 1/2),
                                         shape = 21, fill = "WHITE")
            #zp1 <- zp1 + coord_flip() + theme_bw()
            
            #zp1 <- zp1 + geom_line(data = allModelFrame, aes(linetype =modelName ), size = 1) +
            zp1 <- zp1 + ggtitle("Comparing distance cut-offs")+xlab('Quantile')  
            #zp1 <- zp1 + scale_x_continuous(labels=c("1"="B","2"="M","3"="T"))
            zp1 <- zp1 + scale_x_continuous(breaks=qun)
            print(zp1)  # The trick to these is position_dodge().
            
            ggsave(file=paste(path,'latex/','coeff',meth,treat,spec, '.png', sep=""),height = 7,width =9)
          }
        }
        
      }
    }
  }
}

for(meth in c("lm","gam","match",'sp')){
  p<-get(paste0('ps.',meth,'.did'))
  mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ", ifelse(p < .1, "^\\bullet  ", " "))))
  
  #pb<-exp(get(paste0(statchange,'betas.',meth,'.t.',inf,'.',treat)))-1
  #rpb<-round(pb,3)
  #se<-round(exp(get(paste0(statchange,'ses.',meth,'.t.',inf,'.',treat)))-1,3)
  
  pb<-get(paste0('betas.',meth,'.did'))
  rpb<-round(pb,3)
  se<-round(get(paste0('ses.',meth,'.did')),3)
  
  srpb <- matrix(paste(rpb, mystars, sep=""), ncol=dim(pb)[2] )
  nsrpb<-rbind(c("",laglead),cbind(dist,srpb))
  
  #colnames(srpb)<-laglead
  #rownames(srpb)<-dist
  
  results.mat<-matrix(nrow= 2*dim(srpb)[1],ncol= dim(srpb)[2])
  
  for(i in 1:dim(results.mat)[1]){
    if(i %% 2 != 0){
      results.mat[i,]<-srpb[ceiling(i/2),]
      #    rownames(ols.mat)[i]<-rownames(srpb)[ceiling(i/2)]
    }
    if(i %% 2 == 0){
      results.mat[i,]<-paste0('(',se[ceiling(i/2),],')')
    }
  }
  
  results.mat<-rbind(c('10k','8k','6k','4k','2k'),results.mat)
  if(meth=='lm'){
    results.mat<-rbind(c('OLS','','','',''),results.mat)
  }
  if(meth=='gam'){
    results.mat<-rbind(c('GAM','','','',''),results.mat)
  }
  if(meth=='match'){
    results.mat<-rbind(c('Matching','','','',''),results.mat)
  }
  if(meth=='sp'){
    results.mat<-rbind(c('Spatial Lag','','','',''),results.mat)
  }
  # rn<-c(paste0("(",qcut[1],","),paste0(qcut[2],"]"),paste0("(",qcut[2],","),paste0(qcut[3],"]"),
  #      paste0("(",qcut[3],","),paste0(qcut[4],"]"),
  #     paste0("(",qcut[4],","),paste0(qcut[5],"]"),paste0("(",qcut[5],","),
  #    paste0(qcut[6],"]"),paste0("(",qcut[6],","),paste0(qcut[7],"]"),paste0("(",qcut[7],","),
  #   paste0(qcut[8],"]"),paste0("(",qcut[8],","),paste0(qcut[9],"]"),paste0("(",qcut[9],","),
  #  paste0(qcut[10],"]"),
  # paste0("(",qcut[10],","),paste0(qcut[11],"]"))
  rn<-c(" ","  ", "Total","  ","Municipal Water","   ","Well Water","    ")
  
  #rownames(results.mat)<-rn
  results.mat<-cbind(rn,results.mat)
  xtab<-xtable(results.mat)
  align(xtab) <- "rl|rrrrr"
  print.xtable(xtab,include.rownames=FALSE, hline.after = c(0,1,2,dim(results.mat)[1]),
               include.colnames=FALSE, sanitize.text.function = identity,
               #caption = "example", 
               label = paste0("tab:",meth,"water"),
               type="latex", file=paste0(path,'latex/ATEcomp',meth,'.tex'))
  
  qu<-c("Bottom","Middle","Top")
  if(FALSE){
    allModelFrame <- data.frame(Variable = qu,
                                Coefficient = pb[,1],
                                SE = se[, 1],
                                modelName = "10k")
    for(i in 2:length(dist)){
      di<-dist[i]
      modelFrame <- data.frame(Variable =  qu,
                               Coefficient = pb[,i],
                               SE = se[, i],
                               modelName = di)
      allModelFrame <- data.frame(rbind(allModelFrame,modelFrame))
      
    }
    
    
    interval2 <- -qnorm((1-0.90)/2)  # 95% multiplier
    leg<-c('10k','8k','6k','4k','2k')
    allModelFrame<-allModelFrame[allModelFrame$modelName!="2k",]
    allModelFrame$modelName<-as.factor( allModelFrame$modelName)
    allModelFrame$modelName<-factor(allModelFrame$modelName,levels(allModelFrame$modelName)[c(2:5,1)])
    
    # Plot
    zp1 <- ggplot(allModelFrame, aes(colour = modelName ))
    zp1 <- zp1 + geom_hline(yintercept = 0, colour = gray(1/2), lty = 2)
    zp1 <- zp1 + geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Coefficient - SE*interval2,
                                     ymax = Coefficient + SE*interval2,color =modelName),
                                 lwd = 1/2, position = position_dodge(width = 1/2),
                                 shape = 21, fill = "WHITE")
    #zp1 <- zp1 + coord_flip() + theme_bw()
    
    #zp1 <- zp1 + geom_line(data = allModelFrame, aes(linetype =modelName ), size = 1) +
    zp1 <- zp1 + ggtitle("Comparing distance cut-offs")+xlab('Quantile')  
    #zp1 <- zp1 + scale_x_discrete(breaks=qu)
    print(zp1)  # The trick to these is position_dodge().
    
    ggsave(file=paste(path,'latex/','coeff',meth,treat, '.png', sep=""),height = 7,width =9)
  }
}

for(statchange in c('')){
  for(meth in c('lm',"match",'sp')){
    for(inf in c('did')){
      for(treat in  treatl){
        for(di in dist){
          #treat<-"TATE"
          #di<-1
          #meth<-'lm'
          
          p1<-get(paste0(statchange,'ps.',meth,'.t.',inf,'.',treat,'.',specl[1]))
          p2<-get(paste0(statchange,'ps.',meth,'.t.',inf,'.',treat,'.',specl[2]))
          p3<-get(paste0(statchange,'ps.',meth,'.t.',inf,'.',treat,'.',specl[3]))
          p4<-get(paste0(statchange,'ps.',meth,'.t.',inf,'.',treat,'.',specl[4]))
          p<-rbind(p1,p2,p3,p4)
          mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ", ifelse(p < .1, "^\\bullet  ", " "))))
          if(!is.na(p[1,1])){
            #pb<-exp(get(paste0(statchange,'betas.',meth,'.t.',inf,'.',treat)))-1
            #rpb<-round(pb,3)
            #se<-round(exp(get(paste0(statchange,'ses.',meth,'.t.',inf,'.',treat)))-1,3)
            
            pb1<-get(paste0(statchange,'betas.',meth,'.t.',inf,'.',treat,'.',specl[1]))
            pb2<-get(paste0(statchange,'betas.',meth,'.t.',inf,'.',treat,'.',specl[2]))
            pb3<-get(paste0(statchange,'betas.',meth,'.t.',inf,'.',treat,'.',specl[3]))
            pb4<-get(paste0(statchange,'betas.',meth,'.t.',inf,'.',treat,'.',specl[4]))
            pb<-rbind(pb1,pb2,pb3,pb4)
            
            pb10k<-cbind(pb[1:10,1],pb[11:20,1],pb[21:30,1],pb[31:40,1])
            pb8k<-cbind(pb[1:10,2],pb[11:20,2],pb[21:30,2],pb[31:40,2])
            pb6k<-cbind(pb[1:10,3],pb[11:20,3],pb[21:30,3],pb[31:40,3])
            pb4k<-cbind(pb[1:10,4],pb[11:20,4],pb[21:30,4],pb[31:40,4])
            
            rpb<-round(pb,3)
            se1<-round(get(paste0(statchange,'ses.',meth,'.t.',inf,'.',treat,'.',specl[1])),3)
            se2<-round(get(paste0(statchange,'ses.',meth,'.t.',inf,'.',treat,'.',specl[2])),3)
            se3<-round(get(paste0(statchange,'ses.',meth,'.t.',inf,'.',treat,'.',specl[3])),3)
            se4<-round(get(paste0(statchange,'ses.',meth,'.t.',inf,'.',treat,'.',specl[4])),3)
            se<-rbind(se1,se2,se3,se4)
            
            se10k<-cbind(se[1:10,1],se[11:20,1],se[21:30,1],se[31:40,1])
            se8k<-cbind(se[1:10,2],se[11:20,2],se[21:30,2],se[31:40,2])
            se6k<-cbind(se[1:10,3],se[11:20,3],se[21:30,3],se[31:40,3])
            se4k<-cbind(se[1:10,4],se[11:20,4],se[21:30,4],se[31:40,4])
            
            srpb <- matrix(paste(rpb, mystars, sep=""), ncol=dim(pb)[2] )
            nsrpb<-rbind(c("",laglead),cbind(dist,srpb))
            
            #colnames(srpb)<-laglead
            #rownames(srpb)<-dist
            #10k
            
            
            results.mat<-matrix(nrow= 2*dim(srpb)[1],ncol= dim(srpb)[2])
            
            for(i in 1:dim(results.mat)[1]){
              if(i %% 2 != 0){
                results.mat[i,]<-srpb[ceiling(i/2),]
                #    rownames(ols.mat)[i]<-rownames(srpb)[ceiling(i/2)]
              }
              if(i %% 2 == 0){
                results.mat[i,]<-paste0('(',se[ceiling(i/2),],')')
              }
            }
            
            colnames(results.mat)<-c('10k','8k','6k','4k','2k')
            # rn<-c(paste0("(",qcut[1],","),paste0(qcut[2],"]"),paste0("(",qcut[2],","),paste0(qcut[3],"]"),
            #      paste0("(",qcut[3],","),paste0(qcut[4],"]"),
            #     paste0("(",qcut[4],","),paste0(qcut[5],"]"),paste0("(",qcut[5],","),
            #    paste0(qcut[6],"]"),paste0("(",qcut[6],","),paste0(qcut[7],"]"),paste0("(",qcut[7],","),
            #   paste0(qcut[8],"]"),paste0("(",qcut[8],","),paste0(qcut[9],"]"),paste0("(",qcut[9],","),
            #  paste0(qcut[10],"]"),
            # paste0("(",qcut[10],","),paste0(qcut[11],"]"))
            
            results10k<-cbind(results.mat[1:20,1],results.mat[21:40,1],results.mat[41:60,1],results.mat[61:80,1])
            results8k<-cbind(results.mat[1:20,2],results.mat[21:40,2],results.mat[41:60,2],results.mat[61:80,2])
            results6k<-cbind(results.mat[1:20,3],results.mat[21:40,3],results.mat[41:60,3],results.mat[61:80,3])
            results4k<-cbind(results.mat[1:20,4],results.mat[21:40,4],results.mat[41:60,4],results.mat[61:80,4])
            FEmatrix<- c("Census Tract","Census Tract","Superfund Site","Superfund Site")
            Cluster<-c("Tract by Year","Site by Year","Tract by Year","Site by Year")
            rn2<-c(paste0("(",qcut[1],",",qcut[2],"]"),"  ",paste0("(",qcut[2],",",qcut[3],"]"),"    ",
                   paste0("(",qcut[3],",",qcut[4],"]"),"      ",
                   paste0("(",qcut[4],",",qcut[5],"]"),"       ",paste0("(",qcut[5],",",qcut[6],"]"),"        ",
                   paste0("(",qcut[6],",",qcut[7],"]"),"           ",paste0("(",qcut[7],",",qcut[8],"]"),"            ",
                   paste0("(",qcut[8],",",qcut[9],"]"),"             ", paste0("(",qcut[9],",",qcut[10],"]"),"          ",
                   paste0("(",qcut[10],",",qcut[11],"]"),"                  ", "Fixed Effects", "Cluster")
            
            
            results10k<-rbind(results10k,FEmatrix,Cluster)
            results8k<-rbind(results8k,FEmatrix,Cluster)
            results6k<-rbind(results6k,FEmatrix,Cluster)
            results4k<-rbind(results4k,FEmatrix,Cluster)
            
            rownames(results10k)<-rn2
            rownames(results8k)<-rn2
            rownames(results6k)<-rn2
            rownames(results4k)<-rn2
            
            rn<-c(paste0("(",qcut[1],",",qcut[2],"]"),"Bottom Ten Percentile" ,
                  paste0("(",qcut[2],",",qcut[quant],"]"),"Middle 80 Percentile" ,
                  paste0("(",qcut[quant],",",qcut[quant+1],"]"), "Top Ten Percentile")
            
            xtable(results10k)
            print.xtable(xtable(results10k),include.rownames=TRUE, 
                         include.colnames=FALSE, sanitize.text.function = identity,
                         type="latex", file=paste0(path,'latex/',meth,inf,statchange,treat,"10k.tex"))
            xtable(results8k)
            print.xtable(xtable(results8k),include.rownames=TRUE, 
                         include.colnames=FALSE, sanitize.text.function = identity,
                         type="latex", file=paste0(path,'latex/',meth,inf,statchange,treat,"8k.tex"))
            xtable(results6k)
            print.xtable(xtable(results6k),include.rownames=TRUE, 
                         include.colnames=FALSE, sanitize.text.function = identity,
                         type="latex", file=paste0(path,'latex/',meth,inf,statchange,treat,"6k.tex"))
            xtable(results4k)
            print.xtable(xtable(results4k),include.rownames=TRUE, 
                         include.colnames=FALSE, sanitize.text.function = identity,
                         type="latex", file=paste0(path,'latex/',meth,inf,statchange,treat,"4k.tex"))
            
            
            
            qu<-c("Bottom","Middle","Top")
            leg<-c('Block','Block Group','Tract','Site')
            qun<-c(1:10)
            for(j in c("4k","6k","8k","10k")){
              pb<-get(paste0("pb",j))
              se<-get(paste0("se",j))
              allModelFrame <- data.frame(Variable = qun,
                                          Coefficient = pb[,1],
                                          SE = se[, 1],
                                          modelName = 'Block')
              for(i in 2:length(leg)){
                le<-leg[i]
                modelFrame <- data.frame(Variable =  qun,
                                         Coefficient = pb[,i],
                                         SE = se[, i],
                                         modelName = le)
                allModelFrame <- data.frame(rbind(allModelFrame,modelFrame))
                
              }
              
              interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier
              allModelFrame<-allModelFrame[allModelFrame$modelName!="2k",]
              allModelFrame$modelName<-factor( allModelFrame$modelName,ordered = TRUE)
              allModelFrame$modelName<-factor(allModelFrame$modelName,levels(allModelFrame$modelName)[c(1,2,4,3)],ordered = TRUE)
              
              
              # Plot
              zp1 <- ggplot(allModelFrame, aes(colour = modelName ))
              
              zp1 <- zp1 + geom_hline(yintercept = 0, colour = gray(1/2), lty = 2)
              zp1 <- zp1 + geom_pointrange(aes(x = Variable, y = Coefficient, ymin = Coefficient - SE*interval2,
                                               ymax = Coefficient + SE*interval2,color =modelName),
                                           lwd = 1/2, position = position_dodge(width = 1/2),
                                           shape = 21, fill = "WHITE")
              #zp1 <- zp1 + coord_flip() + theme_bw()
              
              #zp1 <- zp1 + geom_line(data = allModelFrame, aes(linetype =modelName ), size = 1) +
              zp1 <- zp1 + ggtitle("Comparing distance cut-offs")+xlab('Quantile')  
              #zp1 <- zp1 + scale_x_continuous(labels=c("1"="B","2"="M","3"="T"))
              #scale_color_discrete(breaks=c("1","3","10")
              
              zp1 <- zp1 + scale_x_continuous(breaks=qun)
              zp1 <- zp1  +  labs(color="Fixed Effects")
              print(zp1)  # The trick to these is position_dodge().
              
              ggsave(file=paste(path,'latex/','coeff',meth,treat,j, 'ols.png', sep=""),height = 7,width =9)
            }
          }
        }
        
      }
    }
  }
}
