
## Preliminaries
rm(list=ls())

# Change working directory to where you've stored ZTRAX
path<- "P:/Peter/Hedonics/Groundwater/"
path2<- "P:/Peter/Hedonics/"

## This function will check if a package is installed, and if not, install it
pkgTest <- function(x) {
  if (!require(x, character.only = TRUE))
  {
    install.packages(x, dep = TRUE)
    if(!require(x, character.only = TRUE)) stop("Package not found")
  }
}

## These lines load the required packages
packages <- c("readxl","matrixStats","matrixStats","MatchIt","Matching","foreign","multiwayvcov","lmtest","readstata13","xlsx", "data.table","doSNOW","parallel","compare","doParallel","devtools","foreach","spdep","reshape2","sm","plyr","utils","tcltk","geosphere", "matrixcalc", "dplyr","ExPosition", "randomForest","lfe", "hdm", "rdrobust", "stargazer", "ggplot2", "outliers","rpart","e1071")
lapply(packages, pkgTest)


## These lines set several options
options(scipen = 999) # Do not print scientific notation
options(stringsAsFactors = FALSE) ## Do not load strings as factors
options(max.print=1000000)
memory.limit(10000000000000)

## read data
sample.va<- readRDS("P:/Peter/Hedonics/Groundwater/VA_generate_result_fix.rds", refhook = NULL)
#sample.pa<- readRDS("P:/Peter/Hedonics/Groundwater/PA_generate_result_fix.rds", refhook = NULL)
sample.ny<- readRDS("P:/Peter/Hedonics/Groundwater/NY_generate_result_fix.rds", refhook = NULL)
#sample.nj<- readRDS("P:/Peter/Hedonics/Groundwater/NJ_generate_result_fix.rds", refhook = NULL)
#sample.nc<- readRDS("P:/Peter/Hedonics/Groundwater/NC_generate_result_fix.rds", refhook = NULL)
#sample.ma<- readRDS("P:/Peter/Hedonics/Groundwater/MA_generate_result_fix.rds", refhook = NULL)
#sample.mi<- readRDS("P:/Peter/Hedonics/Groundwater/MIHedonics_fix.rds", refhook = NULL)
#sample.de<- readRDS("P:/Peter/Hedonics/Groundwater/DE_Hedonics.rds", refhook = NULL)
#sample.ct<-readRDS(paste0(path2,"Hedonics/CTHedonics_withTract.rds", refhook = NULL))
#sample.de<-readRDS(paste0(path2,"Hedonics/DEHedonics_withTract.rds", refhook = NULL))
sample.ga<-readRDS(paste0(path2,"Hedonics/GAHedonics_withTract.rds", refhook = NULL))
#sample.il<-readRDS(paste0(path2,"Hedonics/ILHedonics_withTract.rds", refhook = NULL))
#sample.in<-readRDS(paste0(path2,"Hedonics/INHedonics_withTract.rds", refhook = NULL))
#sample.md<-readRDS(paste0(path2,"Hedonics/MDHedonics_withTract.rds", refhook = NULL))
#sample.me<-readRDS(paste0(path2,"Hedonics/MEHedonics_withTract.rds", refhook = NULL))
sample.mi<-readRDS(paste0(path2,"Hedonics/MIHedonics_withTract.rds", refhook = NULL))
#sample.nh<-readRDS(paste0(path2,"Hedonics/NHHedonics_withTract.rds", refhook = NULL))
#sample.nj<-readRDS(paste0(path2,"Hedonics/NJHedonics_withTract.rds", refhook = NULL))
sample.oh<-readRDS(paste0(path2,"Hedonics/OHHedonics_withTract.rds", refhook = NULL))
sample.pa<-readRDS(paste0(path2,"Hedonics/PAHedonics_withTract.rds", refhook = NULL))
#sample.ri<-readRDS(paste0(path2,"Hedonics/RIHedonics_withTract.rds", refhook = NULL))
sample.sc<-readRDS(paste0(path2,"Hedonics/SCHedonics_withTract.rds", refhook = NULL))
#sample.vt<-readRDS(paste0(path2,"Hedonics/VTHedonics_withTract.rds", refhook = NULL))
#sample.ct<-readRDS(paste0(path2,"Hedonics/CTHedonics_withTract.rds", refhook = NULL))


table(sample.pa$WaterStndCode)
table(sample.ny$WaterStndCode)
table(sample.nc$WaterStndCode)
table(sample.ma$WaterStndCode)
table(sample.nj$WaterStndCode)
table(sample.va$WaterStndCode)
table(sample.mi$WaterStndCode)
table(sample.ct$WaterStndCode)
table(sample.de$WaterStndCode)
table(sample.ga$WaterStndCode)
table(sample.il$WaterStndCode)
table(sample.in$WaterStndCode)
table(sample.md$WaterStndCode)
table(sample.me$WaterStndCode)
table(sample.nh$WaterStndCode)
table(sample.oh$WaterStndCode)
table(sample.pa$WaterStndCode)
table(sample.ri$WaterStndCode)
table(sample.sc$WaterStndCode)
table(sample.vt$WaterStndCode)

#table(sample.de$WaterStndCode)

sample<-sample.ny
library(gtools)
sample<-smartbind(sample.va,sample.ny)#sample.ct,
#sample.de,
#sample.ga,#sample.il,
#sample.in,sample.md,sample.me,sample.mi,sample.nh,
#sample.oh#sample.mi#sample.pa,sample.ri,
#sample.sc#,sample.vt
#)

saveRDS(sample, file = paste(path2,'Tdrop', sep=""), ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)
#rm(list=ls())
gc()


for(s in c("SC","PA","OH","MI","GA")){
  #s<-"NY"
  sample<-readRDS(paste0(path2,'Hedonics/',s,'Hedonics_withTract.rds', refhook = NULL))
  
  #sample<-readRDS(paste(path2,'tdrop', sep=""), refhook = NULL)
  
  
  ## Prepare Data_________________________________________________________________________________________ 
  
  
  #Manipulating sale price data_________________________________________________________________________
  sample$logprice <- as.numeric(log(sample$SalesPriceAmount)) #calculate log of price
  
  #Exclude houses <$10,000 & >$10,000,000 (Follows Currie et al. (2015))
  sample$price <- as.numeric(sample$SalesPriceAmount)    #
  sample<-sample[sample$price>10000 & sample$price<10000000,]
  
  #sample <- subset(sample, price>10000)  #Exclude houses with sales prices <$10,000
  #sample <- subset(sample, price<10000000)  #Exclude houses with sales prices >$10 mil
  
  #exclude years before 1994
  sample$date <- as.Date(sample$RecordingDate, format="%Y-%m-%d") 
  
  sample<-sample[sample$date>="1994-01-01" & sample$date<="2018-01-01",]
  #sample <- subset(sample, date>="1994-01-01")
  #sample <- subset(sample, date<="2018-01-01")
  
  
  #Subset Housing Price Data_________________________________________________________________________
  #Exclude intrafamily transfers & tax exempt transactions
  sample <- subset(sample, IntraFamilyTransferFlag!="Y")	  
  sample <- subset(sample, TransferTaxExemptFlag!="Y")	
  
  #Keeping only land use type RR101 (most common by far) COME BACK TO WHEN KNOW CODES
  sample <- subset(sample, PropertyLandUseStndCode=="RR101")	   
  
  sample<-sample[!is.na(sample$sqfeet),]
  sample<-sample[!is.na(sample$YearBuilt),]
  sample<-sample[!is.na(sample$LotSizeSquareFeet),]
  sample<-sample[!is.na(sample$PropertyAddressLatitude),]
  
  MU<-subset(sample, WaterStndCode == "MU")
  WL<-subset(sample, WaterStndCode == "WL")
  
  sample<-rbind(MU,WL)
  
  ## save
  saveRDS(sample, file = paste(path,'tdrop', sep=""), ascii = FALSE, version = NULL,
          compress = TRUE, refhook = NULL)
  
  #####################################################################################################################################
  #rm(list=ls())
  gc()
  sample<-readRDS(paste(path,'tdrop', sep=""), refhook = NULL)
  
  # Retain only Repeat Sales
  sample <- sample[(duplicated(sample$ImportParcelID) | duplicated(sample$ImportParcelID, fromLast = TRUE)), ]
  
  #### Houses near NPL sites
  
  #current status of NPL sites keeping only the location
  pNPL=read.csv(paste(path2,'Superfund/Data/ProposedNPLSites.csv', sep=""))
  #pNPL <- subset(pNPL, St=="VA")	
  myvars <- c("Site.ID", "Latitude", "Longitude")
  pNPL <- pNPL[myvars]
  fNPL=read.csv(paste(path2,'Superfund/Data/FinalNPLSites.csv', sep=""))
  #fNPL <- subset(fNPL, St=="VA")
  fNPL <- fNPL[myvars]
  dNPL=read.csv(paste(path2,'Superfund/Data/DeletedNPLSites.csv', sep=""))
  #dNPL <- subset(dNPL, St=="VA")
  dNPL <- dNPL[myvars]
  compNPL<-rbind(pNPL,fNPL,dNPL)
  
  #Sequence of activies for NPL 
  NPL=read.csv(paste(path2,'Superfund/Data/npl1.csv', sep=""))
  #NPL <- subset(NPL1, rstate_code=="VA")
  
  #only look at changes to NPL status
  keep<-c("PROPOSAL TO NPL","FINAL LISTING ON NPL","DELETION FROM NPL","PARTIAL NPL DELETION")
  #keep<-c("FINAL LISTING ON NPL","DELETION FROM NPL","PARTIAL NPL DELETION")
  NPL<-subset(NPL, rat_name %in% keep)
  
  #only take first listings
  NPL <- subset(NPL, rat_act_code!="NP002")
  NPL <- subset(NPL, rat_act_code!="NF002")
  
  #add latitude and longitude ot NPL set
  ll=matrix(nrow=dim(NPL)[1],ncol=3)
  ll[,1]=NPL$site_id
  
  for(i in 1:dim(NPL)[1]){
    for(j in 1:dim(compNPL)[1]){
      if(ll[i,1]==compNPL$Site.ID[j]){
        ll[i,2]<-compNPL$Latitude[j]
        ll[i,3]<-compNPL$Longitude[j]
      }
    }
  }
  
  NPL$lat<-ll[,2]
  NPL$long<-ll[,3]
  
  #vector of contamination at each superfund site
  cNPL=read.csv(paste(path2,'Superfund/Data/ContaminantsinSites1.csv', sep=""))
  #cNPL <- subset(cNPL, rstate_code=="VA")
  
  contamcat=matrix(nrow=dim(cNPL)[1],ncol=3)
  for( i in 1:dim(cNPL)[1]){
    if(cNPL$rmedia_desc[i]=="Surface Water"){
      contamcat[i,1]<-1
    }  
    if(cNPL$rmedia_desc[i]!="Surface Water"){
      contamcat[i,1]<-0
    }
    if(cNPL$rmedia_desc[i]=="Groundwater"){
      contamcat[i,2]<-1
    } 
    if(cNPL$rmedia_desc[i]!="Groundwater"){
      contamcat[i,2]<-0
    }  
    if(cNPL$rmedia_desc[i]=="Air"){
      contamcat[i,3]<-1
    }  
    if(cNPL$rmedia_desc[i]!="Air"){
      contamcat[i,3]<-0
    }
    
  }
  
  cNPL$sw<-contamcat[,1]
  cNPL$gw<-contamcat[,2]
  cNPL$a<-contamcat[,3]
  
  aggsw<-aggregate(cNPL$sw, by=list(Category=cNPL$site_id), FUN=sum)
  sw<-aggsw
  sw[sw > 0] = 1
  
  agggw<-aggregate(cNPL$gw, by=list(Category=cNPL$site_id), FUN=sum)
  gw<-agggw
  gw[gw > 0] = 1
  
  agga<-aggregate(cNPL$a, by=list(Category=cNPL$site_id), FUN=sum)
  a<-agga
  a[a > 0] = 1
  
  agg<-cbind(aggsw[1],sw[2],gw[2],a[2])
  
  cont=matrix(nrow=dim(NPL)[1],ncol=5)
  cont[,4]=NPL$site_id
  cont[,5]=NPL$act_actl_cmpltn_date
  for(i in 1:dim(NPL)[1]){
    for(j in 1:dim(agg)[1])
      if(NPL$site_id[i]-agg[j,1]==0){
        cont[i,1]=agg[j,2]
        cont[i,2]=agg[j,3]
        cont[i,3]=agg[j,4]
      }
  }
  
  NPL$sw=cont[,1]
  NPL$gw=cont[,2]
  NPL$a=cont[,3]
  
  #exclude years before 1994
  cont.df<-as.data.frame(cont)
  cont.df$date <- as.Date(cont.df$V5, format="%m/%d/%Y %H:%M:%S") 
  cont.df <- subset(cont.df, date>="1994-01-01")
  cont.df <- subset(cont.df, date<="2018-01-01")
  
  cont.df<-cont.df[!is.na(cont.df$V1),]
  
  #exclude years before 1994
  NPL$date <- as.Date(NPL$act_actl_cmpltn_date, format="%m/%d/%Y %H:%M:%S") 
  NPL <- subset(NPL, date>="1994-01-01")
  NPL<- subset(NPL, date<="2018-01-01")
  
  NPL<-NPL[!is.na(NPL$gw),]
  ## save
  saveRDS(NPL, file = paste(path,'npl.rds', sep=""), ascii = FALSE, version = NULL,
          compress = TRUE, refhook = NULL)
  
  sitedist.mat<-distm(cbind(sample$PropertyAddressLongitude, sample$PropertyAddressLatitude),cbind(NPL$long,NPL$lat), fun = distHaversine)
  #sample<-cbind(sample,data.frame(sitedist.mat))
  sample$min <- rowMins(sitedist.mat,na.rm = TRUE)
  
  #keep only houses within 5 km of superfund site
  dcut<-10000
  sample <- subset(sample, min<dcut)
  ## save
  saveRDS(sample, file = paste(path,'repeatpre3',s,'.rds', sep=""), ascii = FALSE, version = NULL,
          compress = TRUE, refhook = NULL)
}
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################

sample.SC<-readRDS(paste(path,'repeatpre3SC.rds', sep=""), refhook = NULL)
sample.PA<-readRDS(paste(path,'repeatpre3PA.rds', sep=""), refhook = NULL)
sample.OH<-readRDS(paste(path,'repeatpre3OH.rds', sep=""), refhook = NULL)
sample.MI<-readRDS(paste(path,'repeatpre3MI.rds', sep=""), refhook = NULL)
sample.GA<-readRDS(paste(path,'repeatpre3GA.rds', sep=""), refhook = NULL)
sample.NY<-readRDS(paste(path,'repeatpre3NY.rds', sep=""), refhook = NULL)
NPL<-readRDS(paste(path,'npl.rds', sep=""), refhook = NULL)

library(gtools)
sample<-smartbind(sample.SC,sample.PA,sample.OH,sample.MI,sample.GA,sample.NY
                  #sample.de,
                  #sample.ga,#sample.il,
                  #sample.in,sample.md,sample.me,sample.mi,sample.nh,
                  #sample.oh#sample.mi#sample.pa,sample.ri,
                  #sample.sc#,sample.vt
)
dcut<-10000
#change date format
sample$date <- as.Date(sample$RecordingDate, format="%Y-%m-%d")                           # transform into "R date"
#sample$datepos <- as.POSIXlt(sample$date)
#sample$month <- as.factor(sample$datepos$mon+1)
#sample$year <- as.factor(sample$datepos$year+1900)
#sample$day <- as.factor(sample$datepos$yday+1)

sample$year1994<-ifelse(sample$date<"1995-01-01" & sample$date>="1994-01-01",1,0)
sample$year1995<-ifelse(sample$date<"1996-01-01" & sample$date>="1995-01-01",1,0)
sample$year1996<-ifelse(sample$date<"1997-01-01" & sample$date>="1996-01-01",1,0)
sample$year1997<-ifelse(sample$date<"1998-01-01" & sample$date>="1997-01-01",1,0)
sample$year1998<-ifelse(sample$date<"1999-01-01" & sample$date>="1998-01-01",1,0)
sample$year1999<-ifelse(sample$date<"2000-01-01" & sample$date>="1999-01-01",1,0)
sample$year2000<-ifelse(sample$date<"2001-01-01" & sample$date>="2000-01-01",1,0)
sample$year2001<-ifelse(sample$date<"2002-01-01" & sample$date>="2001-01-01",1,0)
sample$year2002<-ifelse(sample$date<"2003-01-01" & sample$date>="2002-01-01",1,0)
sample$year2003<-ifelse(sample$date<"2004-01-01" & sample$date>="2003-01-01",1,0)
sample$year2004<-ifelse(sample$date<"2005-01-01" & sample$date>="2004-01-01",1,0)
sample$year2005<-ifelse(sample$date<"2006-01-01" & sample$date>="2005-01-01",1,0)
sample$year2006<-ifelse(sample$date<"2007-01-01" & sample$date>="2006-01-01",1,0)
sample$year2007<-ifelse(sample$date<"2008-01-01" & sample$date>="2007-01-01",1,0)
sample$year2008<-ifelse(sample$date<"2009-01-01" & sample$date>="2008-01-01",1,0)
sample$year2009<-ifelse(sample$date<"2010-01-01" & sample$date>="2009-01-01",1,0)
sample$year2010<-ifelse(sample$date<"2011-01-01" & sample$date>="2010-01-01",1,0)
sample$year2011<-ifelse(sample$date<"2012-01-01" & sample$date>="2011-01-01",1,0)
sample$year2012<-ifelse(sample$date<"2013-01-01" & sample$date>="2012-01-01",1,0)
sample$year2013<-ifelse(sample$date<"2014-01-01" & sample$date>="2013-01-01",1,0)
sample$year2014<-ifelse(sample$date<"2015-01-01" & sample$date>="2014-01-01",1,0)
sample$year2015<-ifelse(sample$date<"2016-01-01" & sample$date>="2015-01-01",1,0)
sample$year2016<-ifelse(sample$date<"2017-01-01" & sample$date>="2016-01-01",1,0)

#sample<-sample[, -grep("X", colnames(sample))]
#sample<-sample[, -grep(".1", colnames(sample))]

sample$HHID <- as.factor(sample$ImportParcelID)

#dummy variables indicating in the neighborhood of site contaminating through each vector
sitedist.mat<-distm(cbind(sample$PropertyAddressLongitude, sample$PropertyAddressLatitude),cbind(NPL$long,NPL$lat), fun = distHaversine)

d.sitedist.mat<-sitedist.mat
d.sitedist.mat[d.sitedist.mat<=dcut]=1
d.sitedist.mat[d.sitedist.mat>dcut]=0

gw.cont<-NPL$gw
gw.cont.m<-as.numeric(gw.cont)
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
gw.cont.mat<-rep.row(gw.cont.m,dim(d.sitedist.mat)[1])
d.gw.mat<-d.sitedist.mat*gw.cont.mat

sample<-cbind(sample,data.frame(d.gw.mat))
sumsdgw<-rowSums(d.gw.mat)
dumsdgw<-sumsdgw
dumsdgw[dumsdgw>0]=1

a.cont<-NPL$a
a.cont.m<-as.numeric(a.cont)
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
a.cont.mat<-rep.row(a.cont.m,dim(d.sitedist.mat)[1])
d.a.mat<-d.sitedist.mat*a.cont.mat
sample<-cbind(sample,data.frame(d.a.mat))
sumsda<-rowSums(d.a.mat)
dumsda<-sumsda
dumsda[dumsda>0]=1

#Sequence of activies for NPL 
NPLcontrol1<-read.csv(paste(path2,'Superfund/Data/NPLcontrol2.csv', sep=""))
NPLcontrol1$date <- as.Date(NPLcontrol1$ACT_ACTL_CMPLTN_DATE, format="%d-%b-%y") 
NPLcontrol1<-subset(NPLcontrol1, Engineering.Controls!="No Action")
NPLcontrol1<-subset(NPLcontrol1, Engineering.Controls!="No Further Action")
NPLcontrol<-NPLcontrol1

d <- aggregate(NPLcontrol1$date,by=list(NPLcontrol$SITE_EPA_ID),max)
NPLcontrol <-merge(NPLcontrol,d,by.x=1,by.y=1)
names(NPLcontrol)[names(NPLcontrol) == 'x'] <- 'ControlsComplete'

NPLcontrol2<-subset(NPLcontrol1, Expr1=="Engineering Control")
decc <- aggregate(NPLcontrol2$date,by=list(NPLcontrol2$SITE_EPA_ID),max)
NPLcontrol <-merge(NPLcontrol,decc,by.x=1,by.y=1,all.x=TRUE)
names(NPLcontrol)[names(NPLcontrol) == 'x'] <- 'EngControlsComplete'

NPLcontrol2<-subset(NPLcontrol1, RMEDIA_DESC=="Groundwater")
dgc <- aggregate(NPLcontrol2$date,by=list(NPLcontrol2$SITE_EPA_ID),max)
NPLcontrol <-merge(NPLcontrol,dgc,by.x=1,by.y=1,all.x=TRUE)
names(NPLcontrol)[names(NPLcontrol) == 'x'] <- 'GroundwaterComplete'

NPLcontrol2<-subset(NPLcontrol1, Expr1=="Institutional Control")
dicc <- aggregate(NPLcontrol2$date,by=list(NPLcontrol2$SITE_EPA_ID),max)
NPLcontrol <-merge(NPLcontrol,dicc,by.x=1,by.y=1,all.x=TRUE)
names(NPLcontrol)[names(NPLcontrol) == 'x'] <- 'InstControlsComplete'

NPLcontrol2<-subset(NPLcontrol1, RMEDIA_DESC=="Groundwater")
NPLcontrol2<-subset(NPLcontrol2, Expr1=="Engineering Control")
dgec <- aggregate(NPLcontrol2$date,by=list(NPLcontrol2$SITE_EPA_ID),max)
NPLcontrol <-merge(NPLcontrol,dgec,by.x=1,by.y=1,all.x=TRUE)
names(NPLcontrol)[names(NPLcontrol) == 'x'] <- 'GroundwaterEngComplete'

NPLcontrol2<-subset(NPLcontrol1, Engineering.Controls=="Drinking Water Advisory")
ddwa <- aggregate(NPLcontrol2$date,by=list(NPLcontrol2$SITE_EPA_ID),max)
NPLcontrol <-merge(NPLcontrol,ddwa,by.x=1,by.y=1,all.x=TRUE)
names(NPLcontrol)[names(NPLcontrol) == 'x'] <- 'DrinkingWaterAdvisory'

NPLcontrol2<-subset(NPLcontrol1, Engineering.Controls=="Water Supply Use Restriction")
dwsu <- aggregate(NPLcontrol2$date,by=list(NPLcontrol2$SITE_EPA_ID),max)
NPLcontrol <-merge(NPLcontrol,dwsu,by.x=1,by.y=1,all.x=TRUE)
names(NPLcontrol)[names(NPLcontrol) == 'x'] <- 'WaterSupplyUseRestriction'

NPLcontrol2<-subset(NPLcontrol1, Engineering.Controls=="Deed Notices")
ddn <- aggregate(NPLcontrol2$date,by=list(NPLcontrol2$SITE_EPA_ID),max)
NPLcontrol <-merge(NPLcontrol,ddn,by.x=1,by.y=1,all.x=TRUE)
names(NPLcontrol)[names(NPLcontrol) == 'x'] <- 'DeedNotices'

NPLcontrol2<-subset(NPLcontrol1, Engineering.Controls=="Deed Restriction")
ddn <- aggregate(NPLcontrol2$date,by=list(NPLcontrol2$SITE_EPA_ID),max)
NPLcontrol <-merge(NPLcontrol,ddn,by.x=1,by.y=1,all.x=TRUE)
names(NPLcontrol)[names(NPLcontrol) == 'x'] <- 'DeedRestriction'

NPLcontrol2<-subset(NPLcontrol1, Engineering.Controls=="Groundwater Use/Well Drilling Regulation")
dgur <- aggregate(NPLcontrol2$date,by=list(NPLcontrol2$SITE_EPA_ID),max)
NPLcontrol <-merge(NPLcontrol,dgur,by.x=1,by.y=1,all.x=TRUE)
names(NPLcontrol)[names(NPLcontrol) == 'x'] <- 'GroundwaterUseRegulation'

NPLcontrol2<-subset(NPLcontrol1, Engineering.Controls=="Access Restriction")
dgur <- aggregate(NPLcontrol2$date,by=list(NPLcontrol2$SITE_EPA_ID),max)
NPLcontrol <-merge(NPLcontrol,dgur,by.x=1,by.y=1,all.x=TRUE)
names(NPLcontrol)[names(NPLcontrol) == 'x'] <- 'AccessRestriction'

NPLcontrol2<-subset(NPLcontrol1, Engineering.Controls=="Health Advisory")
dha <- aggregate(NPLcontrol2$date,by=list(NPLcontrol2$SITE_EPA_ID),max)
NPLcontrol <-merge(NPLcontrol,dha,by.x=1,by.y=1,all.x=TRUE)
names(NPLcontrol)[names(NPLcontrol) == 'x'] <- 'HealthAdvisory'


NPLcontrol2<-subset(NPLcontrol1, Engineering.Controls=="Notices to State Regulators Before Changes in Land Ownership")
dha <- aggregate(NPLcontrol2$date,by=list(NPLcontrol2$SITE_EPA_ID),max)
NPLcontrol <-merge(NPLcontrol,dha,by.x=1,by.y=1,all.x=TRUE)
names(NPLcontrol)[names(NPLcontrol) == 'x'] <- 'NoticestoStateRegulators'

myvars<-c("SITE_EPA_ID","GroundwaterEngComplete","ControlsComplete",
          "InstControlsComplete","GroundwaterComplete","EngControlsComplete",
          "HealthAdvisory","GroundwaterUseRegulation","DeedNotices","DeedRestriction",
          #"NoticestoStateRegulators",
          "DrinkingWaterAdvisory","WaterSupplyUseRestriction" #,"AccessRestriction"
)
NPLcont<-NPLcontrol[myvars]
NPLcont<-unique(NPLcont)
names(NPLcont)[names(NPLcont) == 'SITE_EPA_ID'] <- 'site_epa_id'
NPL <- merge(NPL,NPLcont, by="site_epa_id")
myvars <- names(NPL) %in% names(NPL)[grepl( ".x" , names( NPL ) )]
NPL<-NPL[,!myvars]
myvars <- names(NPL) %in% names(NPL)[grepl( ".y.1" , names( NPL ) )]
NPL<-NPL[,!myvars]
NPL<-NPL %>%  rename_all(.funs = funs(sub("\\..*", "", names(NPL)))) 

saveRDS(NPL, file = paste(path,'NPLfull.rds', sep=""), ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

pNPL<-subset(NPL,rat_name=="PROPOSAL TO NPL")
fNPL<-subset(NPL,rat_name=="FINAL LISTING ON NPL")
dNPL<-subset(NPL,rat_name=="DELETION FROM NPL"|rat_name=="PARTIAL NPL DELETION")
pdNPL<-subset(NPL,rat_name=="PARTIAL NPL DELETION")
#ccNPL<-subset(NPL,rat_name=="CONTRUCTION COMPLETE")


opNPL<- pNPL[order(pNPL$date),] 
ofNPL<- fNPL[order(fNPL$date),] 
odNPL<- dNPL[order(dNPL$date),] 
opdNPL<- dNPL[order(pdNPL$date),] 

of.sitedist.mat<-distm(cbind(sample$PropertyAddressLongitude, sample$PropertyAddressLatitude),cbind(ofNPL$long,ofNPL$lat), fun = distHaversine)
colnames(of.sitedist.mat)<-colnames(of.sitedist.mat, do.NULL = FALSE, prefix = "fNPL")
od.sitedist.mat<-distm(cbind(sample$PropertyAddressLongitude, sample$PropertyAddressLatitude),cbind(odNPL$long,odNPL$lat), fun = distHaversine)
colnames(od.sitedist.mat)<-colnames(od.sitedist.mat, do.NULL = FALSE, prefix = "dNPL")
opd.sitedist.mat<-distm(cbind(sample$PropertyAddressLongitude, sample$PropertyAddressLatitude),cbind(opdNPL$long,opdNPL$lat), fun = distHaversine)
colnames(opd.sitedist.mat)<-colnames(opd.sitedist.mat, do.NULL = FALSE, prefix = "pdNPL")


sample<-cbind(sample,opd.sitedist.mat,of.sitedist.mat,od.sitedist.mat)
attach(sample)


for(i in 1:dim(ofNPL)[1]){
  sample[[paste0('aftfinalnpl',i)]]<-ifelse(sample$date - ofNPL$date[i]>0 & get(paste('fNPL',i,sep=""))<dcut,1,0)
  #print(i)
}
#sample<-subset(sample,  select=-c(z,u))


final<-select(sample, starts_with('aftfinalnpl'))

sample <- subset(sample, rowSums(final,na.rm = TRUE)==1)


attach(sample)

deletedCompletion<-c("GroundwaterEngComplete","ControlsComplete",
                     "InstControlsComplete","GroundwaterComplete","EngControlsComplete")
finalInstitutional<-c("HealthAdvisory","GroundwaterUseRegulation","DeedNotices","DeedRestriction",
                      "NoticestoStateRegulators",
                      "DrinkingWaterAdvisory","WaterSupplyUseRestriction","AccessRestriction")
for(i in deletedCompletion){
  for(j in 1:dim(odNPL)[1]){
    treat<-ifelse(sample$date-odNPL[[paste0(i)]][j]>0 &get(paste('dNPL',j,sep=""))<dcut,1,0 )
    treat[is.na(treat)]<-0
    sample[[paste0('treat',i,j)]]<-treat
  }
}
for(i in deletedCompletion){
  for(j in 1:dim(opdNPL)[1]){
    treat<-ifelse(sample$date-opdNPL[[paste0(i)]][j]>0 &get(paste('pdNPL',j,sep=""))<dcut,1,0 )
    treat[is.na(treat)]<-0
    sample[[paste0('treat',i,'p',j)]]<-treat
  }
}

for(i in deletedCompletion){
  treatd<-select(sample, starts_with(paste0('treat',i)))
  treat<-rowSums(treatd)
  treat[treat>0]<-1
  treat[is.na(treat)]<-0
  sample[[paste0('treat',i)]]<-treat
}

for(i in finalInstitutional){
  for(j in 1:dim(ofNPL)[1]){
    treat<-ifelse(sample$date-ofNPL[[paste0(i)]][j]>0 &get(paste('fNPL',j,sep=""))<dcut,1,0 )
    treat[is.na(treat)]<-0
    sample[[paste0('treat',i,j)]]<-treat
  }
}
for(i in finalInstitutional){
  treatd<-select(sample, starts_with(paste0('treat',i)))
  treat<-rowSums(treatd)
  treat[treat>0]<-1
  treat[is.na(treat)]<-0
  sample[[paste0('treat',i)]]<-treat
}

for(i in deletedCompletion){
  for(j in 1:dim(odNPL)[1]){
    treat<-ifelse(sample$date-odNPL[[paste0(i)]][j]<0 &get(paste('dNPL',j,sep=""))<dcut,1,0 )
    treat[is.na(treat)]<-0
    sample[[paste0('pretreat',i,j)]]<-treat
  }
}
for(i in deletedCompletion){
  for(j in 1:dim(opdNPL)[1]){
    treat<-ifelse(sample$date-opdNPL[[paste0(i)]][j]<0 &get(paste('pdNPL',j,sep=""))<dcut,1,0 )
    treat[is.na(treat)]<-0
    sample[[paste0('pretreat',i,'p',j)]]<-treat
  }
}

for(i in deletedCompletion){
  treatd<-select(sample, starts_with(paste0('pretreat',i)))
  treat<-rowSums(treatd)
  treat[treat>0]<-1
  treat[is.na(treat)]<-0
  sample[[paste0('pretreat',i)]]<-treat
}

for(i in finalInstitutional){
  for(j in 1:dim(ofNPL)[1]){
    treat<-ifelse(sample$date-ofNPL[[paste0(i)]][j]<0 &get(paste('fNPL',j,sep=""))<dcut,1,0 )
    treat[is.na(treat)]<-0
    sample[[paste0('pretreat',i,j)]]<-treat
  }
}
for(i in finalInstitutional){
  treatd<-select(sample, starts_with(paste0('pretreat',i)))
  treat<-rowSums(treatd)
  treat[treat>0]<-1
  treat[is.na(treat)]<-0
  sample[[paste0('pretreat',i)]]<-treat
}

for(i in 1:dim(odNPL)[1]){
  sample[[paste0('treatd',i)]]<-ifelse(sample$date - odNPL$date[i]>0 & get(paste('dNPL',i,sep=""))<dcut,1,0)
}

for(i in 1:dim(opdNPL)[1]){
  sample[[paste0('treatpd',i)]]<-ifelse(sample$date - opdNPL$date[i]>0 & get(paste('pdNPL',i,sep=""))<dcut,1,0)
}

treatd<-select(sample, starts_with('treatd'))
treatpd<-select(sample, starts_with('treatpd'))

treatd<-rowSums(cbind(treatd,treatpd))
treatd[treatd>0]<-1

for(i in 1:dim(odNPL)[1]){
  sample[[paste0('treatgroupd',i)]]<-ifelse(get(paste('dNPL',i,sep=""))<dcut,1,0)
}

treatgroupd<-select(sample, starts_with('treatgroupd'))

for(i in 1:dim(opdNPL)[1]){
  sample[[paste0('treatgrouppd',i)]]<-ifelse(get(paste('pdNPL',i,sep=""))<dcut,1,0)
}
treatgrouppd<-select(sample, starts_with('treatgrouppd'))
treatgroup<-rowSums(cbind(treatgroupd,treatgrouppd))

treat<-select(sample, starts_with('treat'))
pretreat<-select(sample, starts_with('pretreat'))

treat[is.na(treat)]<-0
pretreat[is.na(pretreat)]<-0

treatgroupd<-rowSums(cbind(treat,pretreat))
treatgroupd[treatgroupd>0]<-1

control<-ifelse(treatgroupd==1, 0,1)
for(i in 1:dim(odNPL)[1]){
  sample$control<-ifelse(sample$control==1 & sample[[paste0('dNPL',i)]]>20000,1,0)
  
}

sample<-cbind(sample,control)

for(i in 1:dim(ofNPL)[1]){
  sample[[paste0('control',i)]]<-ifelse(sample$control==1 & get(paste('fNPL',i,sep=""))<dcut & ofNPL$rnpl_status_desc[i] == "Currently on the Final NPL",1,0)
}

for(i in 1:dim(odNPL)[1]){
  sample[[paste0('treatdgw',i)]]<-ifelse(sample$date - odNPL$date[i]>0 & odNPL$gw[i]==1 & get(paste('dNPL',i,sep=""))<dcut,1,0)
}

treatdgw<-rowSums(select(sample, starts_with('treatdgw')))
sample<-cbind(sample,treatdgw)

for(i in 1:dim(odNPL)[1]){
  sample[[paste0('treatda',i)]]<-ifelse(sample$date - odNPL$date[i]>0 & odNPL$a[i]==1 & get(paste('dNPL',i,sep=""))<dcut,1,0)
}
treatda<-rowSums(select(sample, starts_with('treatda')))
sample<-cbind(sample,treatda)

for(i in 1:dim(odNPL)[1]){
  sample[[paste0('treatgroupdgw',i)]]<-ifelse(odNPL$gw[i]==1 & get(paste('dNPL',i,sep=""))<dcut,1,0)
}

treatgroupdgw<-rowSums(select(sample, starts_with('treatgroupdgw')))
sample<-cbind(sample,treatgroupdgw)

for(i in 1:dim(odNPL)[1]){
  sample[[paste0('treatgroupda',i)]]<-ifelse(odNPL$a[i]==1 & get(paste('dNPL',i,sep=""))<dcut,1,0)
}

treatgroupda<-rowSums(select(sample, starts_with('treatgroupda')))
sample<-cbind(sample,treatgroupda)

for(i in 1:dim(odNPL)[1]){
  sample[[paste0('treatdn',i)]]<-ifelse(odNPL$gw[i]==0 & sample$date - odNPL$date[i]>0 & get(paste('dNPL',i,sep=""))<dcut,1,0)
}

treatdn<-rowSums(select(sample, starts_with('treatdn')))
sample<-cbind(sample,treatdn)

for(i in 1:dim(odNPL)[1]){
  sample[[paste0('treatgroupdn',i)]]<-ifelse(odNPL$a[i]==1 & get(paste('dNPL',i,sep=""))<dcut,1,0)
}

treatgroupdn<-rowSums(select(sample, starts_with('treatgroupdn')))
sample<-cbind(sample,treatgroupdn)

MU<-subset(sample, WaterStndCode == "MU")
WL<-subset(sample, WaterStndCode == "WL")

sample<-rbind(MU,WL)

WaterStndCode.f <-factor(sample$WaterStndCode)
WaterStndCodeD <- model.matrix(~WaterStndCode.f-1)

sample<-cbind(sample,WaterStndCodeD)

sample$gwcontrol<- sample$control * sample$WaterStndCode.fMU
sample$gwtreat<- sample$treatdgw * sample$WaterStndCode.fWL
sample$gwtreatgroup<-sample$treatgroupdgw * sample$WaterStndCode.fWL

sample$sample<-sample$gwcontrol+sample$gwtreatgroup

for(i in 1:dim(odNPL)[1]){
  sample[[paste0('timefed',i)]]<-ifelse(sample$date - odNPL$date[i]>0,1,0)
}
for(i in 1:dim(opdNPL)[1]){
  sample[[paste0('timefepd',i)]]<-ifelse(sample$date - opdNPL$date[i]>0,1,0)
}



## save
saveRDS(sample, file = paste(path,'repeatpre3.rds', sep=""), ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
sample<-readRDS(paste(path,'repeatpre3.rds', sep=""), refhook = NULL)
NPL<-readRDS(paste(path,'NPLfull.rds', sep=""), refhook = NULL)

pNPL<-subset(NPL,rat_name=="PROPOSAL TO NPL")
fNPL<-subset(NPL,rat_name=="FINAL LISTING ON NPL")
dNPL<-subset(NPL,rat_name=="DELETION FROM NPL"|rat_name=="PARTIAL NPL DELETION")
pdNPL<-subset(NPL,rat_name=="PARTIAL NPL DELETION")

opNPL<- pNPL[order(pNPL$date),] 
ofNPL<- fNPL[order(fNPL$date),] 
odNPL<- dNPL[order(dNPL$date),] 
opdNPL<- dNPL[order(pdNPL$date),] 

deletedCompletion<-c("GroundwaterEngComplete","ControlsComplete",
                     "InstControlsComplete","GroundwaterComplete","EngControlsComplete")
finalInstitutional<-c("HealthAdvisory","GroundwaterUseRegulation","DeedNotices","DeedRestriction",
                      "NoticestoStateRegulators",
                      "DrinkingWaterAdvisory","WaterSupplyUseRestriction","AccessRestriction")
deletedContaminant<-c("gw","a","n")


sample<-sample[!is.na(sample$sqfeet),]
sample<-sample[!is.na(sample$YearBuilt),]
sample<-sample[!is.na(sample$LotSizeSquareFeet),]
sample<-sample[!is.na(sample$PropertyAddressLatitude),]    
sample<-sample[!is.na(sample$NoOfStories),] 
#sample<-sample[!is.na(sample$TotalRooms),] 
#sample<-subset(sample, TotalRooms>0)
sample<-sample[!is.na(sample$TotalBedrooms),] 


d.sample.data<-sample

myvars <- names(d.sample.data) %in% c("Intercept") 
d.sample.data <- d.sample.data[!myvars]

for(i in 1:40){
  d.sample.data[[paste0("binsqfeet",i)]]<-ifelse(d.sample.data$sqfeet-(i*250)<=0 & d.sample.data$sqfeet-((i-1)*250)>0,1,0 )
}

for(i in 1:80){
  d.sample.data[[paste0("binyearbuilt",i)]]<-ifelse(d.sample.data$YearBuilt-(2017-(i*5))>=0 & d.sample.data$YearBuilt-(2017-((i-1)*5))<0 ,1,0)
}

for(i in 1:7){
  d.sample.data[[paste0("binBedroom",i)]]<-ifelse(d.sample.data$TotalBedrooms-((i))==0 ,1,0)
}

for(i in 1:3){
  d.sample.data[[paste0("binStories",i)]]<-ifelse(d.sample.data$NoOfStories-((i))==0 ,1,0)
}

for(i in 1:17){
  d.sample.data[[paste0("binRoom",i)]]<-ifelse(d.sample.data$TotalRooms-((i))==0 ,1,0)
}

for(i in 1:17){
  d.sample.data[[paste0("binLotsize",i)]]<-ifelse(d.sample.data$LotSizeSquareFeet-(i*500)<=0 & d.sample.data$LotSizeSquareFeet-((i-1)*500)>0 ,1,0)
}


PropertyAddressCensusTractAndBlock.f. = factor(d.sample.data$PropertyAddressCensusTractAndBlock)
PropertyAddressCensusTractAndBlock.dummies = model.matrix(~PropertyAddressCensusTractAndBlock.f.)
d.sample.data<- cbind(d.sample.data,PropertyAddressCensusTractAndBlock.dummies)

saveRDS(d.sample.data, file = paste(path,'repeat10.rds', sep=""), ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

#repeat sales Bajari

d.sample.data<-d.sample.data[!duplicated(d.sample.data[,c("date","HHID")]),]

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
A<-rep.row(as.numeric(d.sample.data$HHID),nrow(d.sample.data))
At<-t(A)
D<-At-A
D[D>0]<-2
D[D<0]<-2
D[D==0]<-1
D[D==2]<-0
sameHouse<-D

A<-rep.row(as.numeric(d.sample.data$TransId),nrow(d.sample.data))
At<-t(A)
D<-At-A
D[D>0]<-2
D[D<0]<-2
D[D==0]<-1
D[D==2]<-0
sameSale<-D

otherSales<-sameHouse-sameSale

A<-rep.row(d.sample.data$date,nrow(d.sample.data))
At<-t(A)
D<-At-A
D[D<0]<-0
diffDates<-D

library(matrixStats)
dumDiffDates<-diffDates*otherSales
dumDiffDates[dumDiffDates==0]<- 10000000000000000
#dumDiffDates[dumDiffDates-rowMins(dumDiffDates)>0]<- -1
dumDiffDates[dumDiffDates-rowMins(dumDiffDates,na.rm = TRUE)==0]<-1
dumDiffDates[dumDiffDates<0]<-0
dumDiffDates[dumDiffDates>1]<-0
dumDiffDates[dumDiffDates==10000000000000000]<-0
dumDiffDates[rowSums(dumDiffDates)-dim(dumDiffDates)[1]==0]<-0

d.sample.data$preprice<-dumDiffDates%*%d.sample.data$price
d.sample.data$prelogprice<-dumDiffDates%*%d.sample.data$logprice

for(i in 1:60){
  d.sample.data[[paste0("binPrePrice",i)]]<-ifelse(d.sample.data$preprice-((i-1)*10000)>=0 & d.sample.data$preprice-(i*10000)<0 ,1,0)
}

for(i in 1:33){
  d.sample.data[[paste0("binPreLogPrice",i)]]<-ifelse(d.sample.data$prelogprice-((i-1))>=0 & d.sample.data$prelogprice-(i)<0 ,1,0)
}

d.sample.data$predate<-dumDiffDates%*%as.numeric(d.sample.data$date)
d.sample.data$prediffdate<-as.numeric(d.sample.data$date)-d.sample.data$predate

for(i in 1:44){
  d.sample.data[[paste0("binPreDate",i)]]<-ifelse(d.sample.data$predate-((i-1)*365)>=0 & d.sample.data$predate-(i*365)<0 ,1,0)
}

for(i in 1:15){
  d.sample.data[[paste0("binPreDiffDate",i)]]<-ifelse(d.sample.data$prediffdate-((i-1)*365)>=0 & d.sample.data$prediffdate-(i*365)<0 ,1,0)
}

HHID.f. = factor(d.sample.data$HHID)
HHID.dummies = model.matrix(~HHID.f.)
d.sample.data<- cbind(d.sample.data,HHID.dummies)


saveRDS(d.sample.data, file = paste(path,'repeat10.rds', sep=""), ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

sample<-readRDS(paste(path,'repeat5.rds', sep=""), refhook = NULL)
sampleout<-sample

sample$day<-as.numeric(sample$date)
fmmatdf<-as.formula(~TransId+ date +lagu10l0sp1d1c30tp1 + treatmentgroup+ YearBuilt + NoOfStories + sqfeet)
fmmatch<-as.formula(treatmentgroup ~ date+lagu10l0sp1d1c30tp1 + YearBuilt + NoOfStories+  sqfeet)
for(j in deletedContaminant){ 
  for(i in 1:dim(odNPL)[1]){
    sample$treatmentgroup<-sampleout[[paste0('treatgroupd',j,i)]]  
    #sample[[paste0('treatmentgroup',j)]][is.na(sample[[paste0('treatmentgroup',j)]])]<-0
    sample$control<-ifelse(sampleout$control==1 & sampleout[[paste0('dNPL',i)]]<200000 & sampleout[[paste0('dNPL',i)]]>20000,1,0)
    if(mean(sample$treatmentgroup)>0 & mean(sample$control)>0){
      sample1<-subset(sample, treatmentgroup==1 | control ==1)
      #sample1$treatmentgroup<-sample1[[paste0('treatmentgroup',j,i)]]
      sample1[[paste0('timeFE',j,i)]]<-ifelse(sample1$date-odNPL$date[i]>0 & odNPL$gw[i]==1,1,0)
      sample.df<-data.frame(model.matrix(fmmatdf ,sample1))
      
      saveRDS(sample1, file = paste(path,'fulldeletion',j,i,'.rds', sep=""), ascii = FALSE, version = NULL,
              compress = TRUE, refhook = NULL)
      
      sample.m.df<-subset(sample.df, as.numeric(date-min(sample.df$date[sample.df$treatmentgroup==1]))>0)
      
      mdm<- matchit(fmmatch,data = sample.m.df, method = "nearest", distance = "mahalanobis")
      
      mdm.treat <- match.data(mdm, group = "treat")
      mdm.control <- match.data(mdm, group = "control")
      mdm<-rbind(mdm.treat,mdm.control)
      mdm<-mdm[,c("TransId","lagu10l0sp1d1c30tp1")]
      mdm.full<-merge(mdm,sample1,all.x = TRUE, by = "TransId")
      
      saveRDS(mdm.full , file = paste(path,'mdmdeletion',j,i,'.rds', sep=""), ascii = FALSE, version = NULL,
              compress = TRUE, refhook = NULL)
      
      lm<-lm(fmmatch,sample.m.df)
      
      psm<- matchit(treatmentgroup~ lm$fitted.values,data = sample.m.df, method = "nearest", distance = "mahalanobis")
      
      psm.treat <- match.data(psm, group = "treat")
      psm.control <- match.data(psm, group = "control")
      psm<-rbind(psm.treat,psm.control)
      psm<-psm[,c("TransId","lagu10l0sp1d1c30tp1")]
      psm.full<-merge(psm,sample1,all.x = TRUE, by = "TransId")
      
      saveRDS(psm.full , file = paste(path,'psmdeletion',j,i,'.rds', sep=""), ascii = FALSE, version = NULL,
              compress = TRUE, refhook = NULL)
    }
  }
}
for(j in 1:dim(odNPL)[1]){
  sample$treatmentgroup<-sampleout[[paste0('treatgroupdgw',j)]]  
  #sample[[paste0('treatmentgroup',j)]][is.na(sample[[paste0('treatmentgroup',j)]])]<-0
  sample$control<-ifelse(sampleout$control==1 & sampleout[[paste0('dNPL',j)]]<150000 & sampleout[[paste0('dNPL',j)]]>20000,1,0)
  if(mean(sample$treatmentgroup)>0 & mean(sample$control)>0){
    sample1<-subset(sample, treatmentgroup==1 | control ==1)
    #sample1$treatmentgroup<-sample1[[paste0('treatmentgroup',j)]]
    sample1[[paste0('timeFE',j)]]<-ifelse(sample1$date-odNPL$date[j]>0 & odNPL$gw[j]==1,1,0)
    sample.df<-data.frame(model.matrix(fmmatdf,sample1))
    
    saveRDS(sample1, file = paste(path,'fulldeletion',j,'.rds', sep=""), ascii = FALSE, version = NULL,
            compress = TRUE, refhook = NULL)
    
    sample.m.df<-subset(sample.df, as.numeric(date-min(sample.df$date[sample.df$treatmentgroup==1]))>0)
    
    mdm<- matchit(fmmatch,data = sample.m.df, method = "nearest", distance = "mahalanobis")
    
    mdm.treat <- match.data(mdm, group = "treat")
    mdm.control <- match.data(mdm, group = "control")
    mdm<-rbind(mdm.treat,mdm.control)
    mdm<-mdm[,c("TransId","lagu10l0sp1d1c30tp1")]
    mdm.full<-merge(mdm,sample1,all.x = TRUE, by = "TransId")
    
    saveRDS(mdm.full , file = paste(path,'mdmdeletion',j,'.rds', sep=""), ascii = FALSE, version = NULL,
            compress = TRUE, refhook = NULL)
    
    lm<-lm(fmmatch,sample.m.df)
    
    psm<- matchit(treatmentgroup~ lm$fitted.values,data = sample.m.df, method = "nearest", distance = "mahalanobis")
    
    psm.treat <- match.data(psm, group = "treat")
    psm.control <- match.data(psm, group = "control")
    psm<-rbind(psm.treat,psm.control)
    psm<-psm[,c("TransId","lagu10l0sp1d1c30tp1")]
    psm.full<-merge(psm,sample1,all.x = TRUE, by = "TransId")
    
    saveRDS(psm.full , file = paste(path,'psmdeletion',j,'.rds', sep=""), ascii = FALSE, version = NULL,
            compress = TRUE, refhook = NULL)
  }
}

for(j in deletedCompletion){
  for(i in 1:dim(odNPL)[1]){
    sample$treatmentgroup<-sampleout[[paste0('treat',j,i)]]+sampleout[[paste0('pretreat',j,i)]]  
    sample$control<-ifelse(sampleout$control==1 & sampleout[[paste0('dNPL',i)]]<150000 & sampleout[[paste0('dNPL',i)]]>20000,1,0)
    if(mean(sample$treatmentgroup)>0 & mean(sample$control)>0){
      sample1<-subset(sample, treatmentgroup==1 | control ==1)
      #sample1$treatmentgroup<-sample1[[paste0('treatmentgroup',j,i)]]
      sample1[[paste0('timeFE',j,i)]]<-ifelse(sample1$date-odNPL$date[i]>0 & odNPL$gw[i]==1,1,0)
      
      sample.df<-data.frame(model.matrix(fmmatdf,sample1))
      
      saveRDS(sample1, file = paste(path,'full',j,i,'.rds', sep=""), ascii = FALSE, version = NULL,
              compress = TRUE, refhook = NULL)
      
      sample.m.df<-subset(sample.df, as.numeric(date-min(sample.df$date[sample.df$treatmentgroup==1]))>0)
      
      mdm<- matchit(fmmatch,  data = sample.m.df, method = "nearest", distance = "mahalanobis")
      
      mdm.treat <- match.data(mdm, group = "treat")
      mdm.control <- match.data(mdm, group = "control")
      mdm<-rbind(mdm.treat,mdm.control)
      mdm<-mdm[,c("TransId","lagu10l0sp1d1c30tp1")]
      
      mdm.full<-merge(mdm,sample1,all.x = TRUE, by = "TransId")
      
      saveRDS(mdm.full , file = paste(path,'mdm',j,i,'.rds', sep=""), ascii = FALSE, version = NULL,
              compress = TRUE, refhook = NULL)
      lm<-lm(fmmatch,sample.m.df)
      
      psm<- matchit(treatmentgroup~ lm$fitted.values,data = sample.m.df, method = "nearest", distance = "mahalanobis")
      
      psm.treat <- match.data(psm, group = "treat")
      psm.control <- match.data(psm, group = "control")
      psm<-rbind(psm.treat,psm.control)
      psm<-psm[,c("TransId","lagu10l0sp1d1c30tp1")]
      psm.full<-merge(psm,sample1,all.x = TRUE, by = "TransId")
      
      saveRDS(psm.full , file = paste(path,'psm',j,i,'.rds', sep=""), ascii = FALSE, version = NULL,
              compress = TRUE, refhook = NULL)
    }
  }
}
for(j in finalInstitutional){
  for(i in 1:dim(ofNPL)[1]){
    sample$treatmentgroup<-sampleout[[paste0('treat',j,i)]]+sampleout[[paste0('pretreat',j,i)]]  
    sample$control<-ifelse(sampleout$control==1 & sampleout[[paste0('fNPL',i)]]<150000 & sampleout[[paste0('fNPL',i)]]>20000,1,0)
    if(mean(sample$treatmentgroup)>0 & mean(sample$control)>0){
      sample1<-subset(sample, treatmentgroup==1 | control ==1)
      #sample1$treatmentgroup<-sample1[[paste0('treatmentgroup',j,i)]]
      sample1[[paste0('timeFE',j,i)]]<-ifelse(sample1$date-ofNPL$date[i]>0 & ofNPL$gw[i]==1,1,0)
      
      sample.df<-data.frame(model.matrix(fmmatdf,sample1))
      
      saveRDS(sample1, file = paste(path,'full',j,i,'.rds', sep=""), ascii = FALSE, version = NULL,
              compress = TRUE, refhook = NULL)
      
      sample.m.df<-subset(sample.df, as.numeric(date-min(sample.df$date[sample.df$treatmentgroup==1]))>0)
      
      mdm<- matchit(fmmatch,data = sample.m.df, method = "nearest", distance = "mahalanobis")
      
      mdm.treat <- match.data(mdm, group = "treat")
      mdm.control <- match.data(mdm, group = "control")
      mdm<-rbind(mdm.treat,mdm.control)
      mdm<-mdm[,c("TransId","lagu10l0sp1d1c30tp1")]
      
      mdm.full<-merge(mdm,sample1,all.x = TRUE, by = "TransId")
      
      saveRDS(mdm.full , file = paste(path,'mdm',j,i,'.rds', sep=""), ascii = FALSE, version = NULL,
              compress = TRUE, refhook = NULL)
      
      lm<-lm(fmmatch,sample.m.df)
      
      psm<- matchit(treatmentgroup~ lm$fitted.values,data = sample.m.df, method = "nearest", distance = "mahalanobis")
      
      psm.treat <- match.data(psm, group = "treat")
      psm.control <- match.data(psm, group = "control")
      psm<-rbind(psm.treat,psm.control)
      psm<-psm[,c("TransId","lagu10l0sp1d1c30tp1")]
      psm.full<-merge(psm,sample1,all.x = TRUE, by = "TransId")
      
      saveRDS(psm.full , file = paste(path,'psm',j,i,'.rds', sep=""), ascii = FALSE, version = NULL,
              compress = TRUE, refhook = NULL)
    }
  }
}



