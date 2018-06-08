
## Preliminaries
rm(list=ls())
gc()

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
sample<- readRDS("P:/Peter/Hedonics/Groundwater/NY_generate_result_fix.rds", refhook = NULL)
#sample<- readRDS("P:/Peter/Hedonics/Groundwater/PA_generate_result_fix.rds", refhook = NULL)
#cosample<- readRDS(paste0(path2,"Hedonics/COHedonics_withTract.rds"), refhook = NULL)

table(sample$WaterStndCode)

s<-"NY"

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
#sample <- sample[(duplicated(sample$ImportParcelID) | duplicated(sample$ImportParcelID, fromLast = TRUE)), ]
#NPL<-readRDS(paste(path2,'NPLfull.rds', sep=""), refhook = NULL)
ccNPL<-read.csv(paste(path2,'Superfund/Data/nplconcomp.csv', sep=""))
ccNPL<-ccNPL[ c(FALSE,TRUE), ]
ccNPL$ConstructionComplete<-as.Date(ccNPL$Construction.1, format="%m/%d/%Y")
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
keepl<-c("DISCOVERY","PROPOSAL","FINAL","DELETION")
keepcodel<-c("DS001","NP001","NF001","ND001")
#keep<-c("FINAL LISTING ON NPL","DELETION FROM NPL","PARTIAL NPL DELETION")
unpl<-subset(NPL , rat_act_code=="NP001")
#unpl<-NPL[unique(NPL$site_name),]
NPL$Date<-as.Date(strptime(NPL$act_actl_cmpltn_date, format="%m/%d/%Y %H:%M:%S"))
mydate = strptime('16/Oct/2005:07:51:00',format='%d/%b/%Y:%H:%M:%S')

disnpl<-NPL[NPL$rat_act_code=="DS001",c("site_epa_id","site_id","rsite_section_name","site_name","rsitinc_desc","rat_act_code","Date")]


for(i in names(table(NPL$ralt_desc))[-1]){
  #i<-"EPA Fund-Financed"
  disnpl[[paste0(i)]]<-0
  for(j in names(table(NPL$site_epa_id))){
    #j<-names(table(NPL$site_epa_id))[1]
    if(i %in% names(table(NPL[NPL$site_epa_id==j,"ralt_desc"]))){
      disnpl[disnpl$site_epa_id==j,paste0(i)]<-1
    }
  }
}


pronpl<-NPL[NPL$rat_act_code=="NP001",c("site_epa_id","Date")]
finnpl<-NPL[NPL$rat_act_code=="NF001",c("site_epa_id","Date")]
delnpl<-NPL[NPL$rat_act_code=="ND001",c("site_epa_id","Date")]

npldata<-merge(disnpl,pronpl,by="site_epa_id",all.x =TRUE)
names(npldata)[names(npldata)=="Date.x"]<-"discoveryDate"
names(npldata)[names(npldata)=="Date.y"]<-"proposalDate"
npldata<-merge(npldata,finnpl,by="site_epa_id",all.x =TRUE)
names(npldata)[names(npldata)=="Date"]<-"finalDate"
npldata<-merge(npldata,delnpl,by="site_epa_id",all.x =TRUE)
names(npldata)[names(npldata)=="Date"]<-"deleteDate"
ccNPLn<-ccNPL[,c("Site.EPA.ID","ConstructionComplete")]
npldata<-merge(npldata,ccNPLn,by.x="site_epa_id",by.y="Site.EPA.ID",all.x =TRUE)
names(npldata)[names(npldata)=="ConstructionComplete"]<-"ccDate"


npldata$dtoprop[!is.na(npldata$proposalDate)]<-npldata$proposalDate[!is.na(npldata$proposalDate)]-
  npldata$discoveryDate[!is.na(npldata$proposalDate)]
npldata$dtofin[!is.na(npldata$finalDate)]<-npldata$finalDate[!is.na(npldata$finalDate)]-
  npldata$proposalDate[!is.na(npldata$finalDate)]
npldata$dtocc[!is.na(npldata$ccDate)]<-npldata$ccDate[!is.na(npldata$ccDate)]-
  npldata$finalDate[!is.na(npldata$ccDate)]
npldata$dtodel[!is.na(npldata$deleteDate)]<-npldata$deleteDate[!is.na(npldata$deleteDate)]-
  npldata$ccDate[!is.na(npldata$deleteDate)]

npldata$ytoprop<-npldata$dtoprop/365
npldata$ytofin<-npldata$dtofin/365
npldata$ytocc<-npldata$dtocc/365
npldata$ytodel<-npldata$dtodel/365

summary(npldata$ytoprop)
summary(npldata$ytofin)
summary(npldata$ytocc)
summary(npldata$ytodel)
npldata<-npldata[npldata$ytodel>0 | is.na(npldata$ytodel),]
npldata<-npldata[npldata$ytoprop>0| is.na(npldata$ytoprop),]
#table(ifelse(npldata$ytodel<0,1,0))
NPL<-npldata

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
#cont[,5]=NPL$act_actl_cmpltn_date
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

#NPL<-NPL[!is.na(NPL$gw),]
## save
saveRDS(NPL, file = paste(path,'npl.rds', sep=""), ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)
NPL<-NPL[NPL$rsitinc_desc=="LANDFILL",]
sitedist.mat<-distm(cbind(sample$PropertyAddressLongitude, sample$PropertyAddressLatitude),cbind(NPL$long,NPL$lat), fun = distHaversine)
#sample<-cbind(sample,data.frame(sitedist.mat))
sample$min <- rowMins(sitedist.mat,na.rm = TRUE)

#keep only houses within 5 km of superfund site
#dcut<-15000
#sample <- subset(sample, min<dcut)
## save
saveRDS(sample, file = paste(path,'repeatpre3',s,'.rds', sep=""), ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

###########################################################################################################################
###########################################################################################################################
###########################################################################################################################

#sample<-readRDS(paste(path,'repeatpre3PA.rds', sep=""), refhook = NULL)
sample<-readRDS(paste(path,'repeatpre3NY.rds', sep=""), refhook = NULL)
NPL<-readRDS(paste(path,'npl.rds', sep=""), refhook = NULL)
NPL<-NPL[NPL$rsitinc_desc=="LANDFILL",]
NPL<-NPL[substr(NPL$site_epa_id, 1, 2)=="NY",]

#NPLct<-NPL[NPL$rstate_code=="CT",]
#NPLvt<-NPL[NPL$rstate_code=="VT",]
#NPLma<-NPL[NPL$rstate_code=="MA",]
#NPLpa<-NPL[NPL$rstate_code=="PA",]

#NPL<-rbind(NPLny,NPLct,NPLvt,NPLma,NPLpa)
dcut<-15000
sample <- subset(sample, min<dcut)
sample$dist14k<-ifelse(sample$min-14000<0,1,0)
sample$dist12k<-ifelse(sample$min-12000<0,1,0)
sample$dist10k<-ifelse(sample$min-10000<0,1,0)
sample$dist8k<-ifelse(sample$min-8000<0,1,0)
sample$dist6k<-ifelse(sample$min-6000<0,1,0)
sample$dist5k<-ifelse(sample$min-5000<0,1,0)
sample$dist4k<-ifelse(sample$min-4000<0,1,0)
sample$dist3k<-ifelse(sample$min-3000<0,1,0)
sample$dist2k<-ifelse(sample$min-2000<0,1,0)
sample$dist1k<-ifelse(sample$min-1000<0,1,0)
sample$dist500m<-ifelse(sample$min-500<0,1,0)
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
colnames(sitedist.mat) <- 1:(dim(sitedist.mat)[2])
colnames(sitedist.mat) <- paste("dist_to_site", colnames(sitedist.mat), sep = "_")
sample<-cbind(sample,data.frame(sitedist.mat))
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

#sample<-cbind(sample,data.frame(d.gw.mat))
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
#sample<-cbind(sample,data.frame(d.a.mat))
sumsda<-rowSums(d.a.mat)
dumsda<-sumsda
dumsda[dumsda>0]=1

#Sequence of activies for NPL 
NPLcontrol1<-read.csv(paste(path2,'Superfund/Data/NPLcontrol2.csv', sep=""))
NPLcontrol1$date <- as.Date(NPLcontrol1$ACT_ACTL_CMPLTN_DATE, format="%d-%b-%y") 
#NPLcontrol1<-subset(NPLcontrol1, Engineering.Controls!="No Action")
#NPLcontrol1<-subset(NPLcontrol1, Engineering.Controls!="No Further Action")
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

saveRDS(NPL, file = paste(path,'NPLfullny.rds', sep=""), ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

sample<-sample[!is.na(sample$sqfeet),]
sample<-sample[!is.na(sample$YearBuilt),]
sample<-sample[!is.na(sample$LotSizeSquareFeet),]
sample<-sample[!is.na(sample$PropertyAddressLatitude),]    
sample<-sample[!is.na(sample$NoOfStories),] 
sample<-sample[!is.na(sample$TotalRooms),] 
#sample<-subset(sample, TotalRooms>0)
sample<-sample[!is.na(sample$TotalBedrooms),] 



attach(sample)

tNPL<-NPL[!is.na(NPL$deleteDate),]
tsites<-row.names(tNPL)
tsitesdist<-paste('dist_to_site',tsites,sep = "_")
sample<-data.frame(sample)
dist_to_site_mint<-rowMins(as.matrix(sample[, tsitesdist]))

closestt<-NA
for(i in tsites){
  closestt<-ifelse(dist_to_site_mint-sample[[paste0('dist_to_site_',i)]]==0,i,closestt)
}
sample$closesttsite<-closestt

#tNPL<-NPL[!is.na(NPL$deleteDate),]
sites<-row.names(NPL)
sitesdist<-paste('dist_to_site',sites,sep = "_")
sample<-data.frame(sample)
dist_to_site_min<-rowMins(as.matrix(sample[, sitesdist]))

closest<-NA
for(i in sites){
  closest<-ifelse(dist_to_site_min-sample[[paste0('dist_to_site_',i)]]==0,i,closest)
}
sample$closestsite<-closest
dcut<-10000
tgdel<-rep(0,dim(sample)[1])
postdel<-rep(0,dim(sample)[1])
timetodel<-rep(0,dim(sample)[1])
cdel<-rep(0,dim(sample)[1])
#tgprop<-NA
for(i in as.numeric(tsites)){
  postdel<-ifelse(sample$date - NPL$deleteDate[i]>0 & as.numeric(sample$closesttsite)-i==0,1,postdel)
  timetodel<-ifelse( as.numeric(sample$closesttsite)-i==0,sample$date - NPL$deleteDate[i],timetodel)
  tgdel<-ifelse(sample[[paste0('dist_to_site_',i)]]<dcut& as.numeric(sample$closesttsite)-i==0,1,tgdel)
  cdel<-ifelse(sample[[paste0('dist_to_site_',i)]]>30000& as.numeric(sample$closesttsite)-i==0,1,cdel)
  
}
sample$postdel<-postdel
sample$timetodel<-timetodel
sample$tgdel<-tgdel
sample$cdel<-cdel
sample<-sample[sample$tgdel>0|sample$cdel>0,]
#sample<-sample[sample$tgprop<2,]
postprop<-rep(0,dim(sample)[1])
postfin<-rep(0,dim(sample)[1])
postcc<-rep(0,dim(sample)[1])

for(i in 1:dim(NPL)[1]){
  postprop<-ifelse(sample$date - NPL$proposalDate[i]>0 & as.numeric(sample$closesttsite)-i==0,1,postprop)
  postfin<-ifelse(sample$date - NPL$finalDate[i]>0 & as.numeric(sample$closesttsite)-i==0,1,postfin)
  postcc<-ifelse(sample$date - NPL$ccDate[i]>0 & as.numeric(sample$closesttsite)-i==0,1,postcc)
  
}
sample$postprop<-postprop
sample$postfin<-postfin
sample$postcc<-postcc


WaterStndCode.f <-factor(sample$WaterStndCode)
WaterStndCodeD <- model.matrix(~WaterStndCode.f-1)

sample<-cbind(sample,WaterStndCodeD)

sample$gwcdel<- sample$cdel* sample$WaterStndCode.fWL
sample$gwtgdel<- sample$tgdel * sample$WaterStndCode.fWL

sample$mucdel<- sample$cdel* sample$WaterStndCode.fMU
sample$mutgdel<- sample$tgdel * sample$WaterStndCode.fMU

#sample1<-sample[dist_to_site_min<10000,]
#sample1<-sample1[]
## save
saveRDS(sample, file = paste(path,'repeatpr5.rds', sep=""), ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

###########################################################################################################################
###########################################################################################################################
###########################################################################################################################

## Preliminaries
rm(list=ls())
gc()

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

packages <- c("readxl","Hmisc","DescTools","qgam","quantreg","sphet","mgcv","McSpatial","pastecs","rdd","Matrix","psych","xtable","splines","ck37r","data.table","matrixStats","tmle","xgboost", "MatchIt","gtools","statar","foreign","multiwayvcov","lmtest","readstata13","xlsx", "data.table","doSNOW","parallel","compare","doParallel","devtools","foreach","spdep","reshape2","sm","plyr","utils","tcltk","geosphere", "matrixcalc", "dplyr","ExPosition", "randomForest","lfe", "hdm", "rdrobust", "stargazer", "ggplot2", "outliers","rpart","e1071")
lapply(packages, pkgTest)
packages <- c("readxl","Hmisc","grf","classInt","RColorBrewer","rgdal","DescTools","qgam","quantreg","sphet","mgcv","McSpatial","pastecs","rdd","Matrix","psych","xtable","splines","ck37r","data.table","matrixStats","tmle","xgboost", "MatchIt","gtools","statar","foreign","multiwayvcov","lmtest","readstata13","xlsx", "data.table","doSNOW","parallel","compare","doParallel","devtools","foreach","spdep","reshape2","sm","plyr","utils","tcltk","geosphere", "matrixcalc", "dplyr","ExPosition", "randomForest","lfe", "hdm", "rdrobust", "stargazer", "ggplot2", "outliers","rpart","e1071")
lapply(packages, pkgTest)


#Shape file
cblock <- readOGR("P:/Peter/Hedonics/Groundwater/census","tl_2010_36_tabblock10")
cblock <- spTransform(cblock,"+proj=longlat +datum=WGS84")
cbg <- readOGR("P:/Peter/Hedonics/Groundwater/census","tl_2010_36_bg10")
cbg <- spTransform(cbg,"+proj=longlat +datum=WGS84")
ctract <- readOGR("P:/Peter/Hedonics/Groundwater/census","tl_2010_36_tract10")
ctract <- spTransform(ctract,"+proj=longlat +datum=WGS84")


NPL<-readRDS(paste(path,'NPLfullny.rds', sep=""), refhook = NULL)

sample<-readRDS(paste(path,'repeatpr5.rds', sep=""),refhook = NULL)

samplefull<-sample

samplefull$la<-samplefull$PropertyAddressLatitude
samplefull$lo<-samplefull$PropertyAddressLongitude
cblock$block <-as.numeric(as.character(cblock$GEOID10))
cbg$bg<-as.numeric(as.character(cbg$GEOID10))
ctract$tract<-as.numeric(as.character(ctract$GEOID10))
lmat <- samplefull[,c("lo","la")]
spdata <- SpatialPointsDataFrame(lmat,lmat)
proj4string(spdata) <- "+proj=longlat +datum=WGS84"
samplefull$ctr <- over(spdata,ctract)$tract
samplefull$cbg <- over(spdata,cbg)$bg
samplefull$cbl <- over(spdata,cblock)$block

saveRDS(samplefull, file = paste(path,'fullcen.rds', sep=""), ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)
###############################################################################################


sample<-readRDS(paste(path,'fullcen.rds', sep=""),refhook = NULL)
NPL<-readRDS(paste(path,'NPLfullny.rds', sep=""), refhook = NULL)

#i<-tracts[2]
tracts<-names(table(sample$ctr))
sample2<-data.table(sample)
for(i in tracts){
gc()
sample1<-sample[as.numeric(sample$ctr)-as.numeric(i)==0,]
ud<-!duplicated(sample1[,c("date","HHID")])
sample1<-sample1[ud,]
sample1<-sample1[!is.na(sample1$ctr),]
d.sample.data<-sample1[,c("date","HHID","TransId","price","logprice")]
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
D<-rep.row(as.numeric(d.sample.data$HHID),nrow(d.sample.data))
D<-t(D)-D
D[D>0]<-2
D[D<0]<-2
D[D==0]<-1
D[D==2]<-0
sameHouse<-D

D<-rep.row(as.numeric(d.sample.data$TransId),nrow(d.sample.data))

D<-t(D)-D
D[D>0]<-2
D[D<0]<-2
D[D==0]<-1
D[D==2]<-0
sameSale<-D

otherSales<-sameHouse-sameSale
rm(sameHouse,sameSale)
gc()
D<-rep.row(d.sample.data$date,nrow(d.sample.data))
D<-t(D)-D
D[D<0]<-0
diffDates<-D
rm(D)
gc()

library(matrixStats)
dumDiffDates<-diffDates*otherSales
dumDiffDates[dumDiffDates==0]<- 10000000000000000
#dumDiffDates[dumDiffDates-rowMins(dumDiffDates)>0]<- -1
dumDiffDates[dumDiffDates-rowMins(dumDiffDates,na.rm = TRUE)==0]<-1
dumDiffDates[dumDiffDates<0]<-0
dumDiffDates[dumDiffDates>1]<-0
dumDiffDates[dumDiffDates==10000000000000000]<-0
dumDiffDates[rowSums(dumDiffDates)-dim(dumDiffDates)[1]==0]<-0
rm(diffDates,otherSales)
gc()


d.sample.data$preprice<-dumDiffDates%*%d.sample.data$price
d.sample.data$prelogprice<-dumDiffDates%*%d.sample.data$logprice
d.sample.data$predate<-dumDiffDates%*%as.numeric(d.sample.data$date)
d.sample.data$prediffdate<-as.numeric(d.sample.data$date)-d.sample.data$predate


d.sample.data<-data.table(d.sample.data)
if(i ==tracts[1]){
  samplepre<-d.sample.data
}
if(i !=tracts[1]){
  samplepre<-rbind(samplepre, d.sample.data, fill = TRUE)
}
#sample2<-sample


print(i)
}
sample2<-merge(sample2,samplepre[,c("date","HHID","TransId",
                                        "preprice","prelogprice","predate","prediffdate")],
               by = c("TransId","date","HHID"),all.x=TRUE)       
sample2<-sample2[preprice>0,]

saveRDS(sample2, file = paste0(path,'baj.rds'), ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)
###############################################


NPL<-readRDS(paste(path,'NPLfullny.rds', sep=""), refhook = NULL)

samplefull<-readRDS(paste(path,'baj.rds', sep=""),refhook = NULL)

sitelist<-names(table(samplefull$closestsite))

for(i in sitelist){
  
 sample<-samplefull[samplefull$closestsite==i,]
 #sample$treatst<-sample$treatst-sample$presstatusd
  #sample<-sample[treatst>0,]
  #sample<-sample[presstatusd==0,]
  upper.spatial.range<-c(20)
  lower.spatial.range<-c(0)
  spatial.power.range<-c(10)
  temporal.cut.range<-c(30)
  temporal.power.range<-c(10)
  
  urange<-upper.spatial.range
  lrange<-lower.spatial.range
  prange<-spatial.power.range
  crange<-temporal.cut.range
  qrange<-temporal.power.range
  
  d.sample.data<-sample
  W.trend.lag.variables<-function(urange,lrange,prange,crange,qrange,path){
    
    
    denom<-0 
    for(c in crange){
      for(q in qrange){
        denom<-denom+1
      }
    }   
    for(u in urange){
      for(l in lrange){
        for(p in prange){
          if(u>l){
            denom<-denom+1
          }
          
        }
      }
    }
    
    
    dist.mat<-distm (cbind(d.sample.data$PropertyAddressLongitude, d.sample.data$PropertyAddressLatitude), fun = distHaversine)
    
    Wtime<- function(c,q){  
      t<-d.sample.data$predate
      rep.row<-function(x,n){
        matrix(rep(x,each=n),nrow=n)
      }
      t<-d.sample.data$predate
      T1<-rep.row(t,nrow(d.sample.data))
      T2<-T1
      T1<-t(T1)
      Tdp<-t(t(T1-T2))
      Tdp[Tdp >= c*365] = -2
      Tdp[Tdp<=-365] = -2
      Wtp<-1/(1+abs(Tdp))
      Wtp[Wtp<0]<-0
      #Wtp[Wtp == 'Inf'] = 0
      Wt<-Wtp^q
      Wt<-return(Wt)
      print(Wt)
      rm(t,T1,T2,Tdp,Wtp)
      gc()
    }
    
    
    num<-0
    for(c in crange){
      for(q in qrange){
        
        assign(paste('wt',c,'m',q,sep=""),Wtime(c,q))
        num<-num+1
        print(paste(num, 'of', denom,sep=" "))
      }
    }
    
    
    Wspat<- function(u,l,p){  
      
      Sp<-dist.mat
      Sp[Sp>= u*500] = -2
      Sp[Sp <= l*500] = -2
      Wsp<-1/(1+Sp)
      Wsp[Wsp<0]<-0
      #Wsp[Wsp == 'Inf'] = 0
      Ws<-Wsp^p
      Ws<-return(Ws)
      print(Ws)
      rm(Sp,Wsp)
      gc()
    }
    
    
    
    for(u in urange){
      
      for(l in lrange){
        for(p in prange){
          assign(paste('ws',u,'m',l, 'm',p,sep=""),Wspat(u,l,p))
          if(u>l){
            num<-num+1
            print(paste(num,'of',denom,sep=" "))
          }
          else{num<-num}
        }
      }
    }
    
    #Dummies
    
    rep.row<-function(x,n){
      matrix(rep(x,each=n),nrow=n)
    }
    if(FALSE){
      A<-rep.row(d.sample.data$YearBuilt,nrow(d.sample.data))
      At<-t(A)
      D<-At-A
      D[D<10000000000]<-1
      D1<-D
    }
    
    Wst<-hadamard.prod(get(paste('ws',u,'m',l, 'm',p,sep="")), get(paste('wt',c,'m',q,sep="")))
    #Wst<-Wst*get(paste('D',d, sep=""))
    weight<-rowSums(Wst)
    for(i in 1:nrow(d.sample.data)){
      if(weight[i]==0){
        weight[i]<-1
      }
    }
    Wst<-Wst*(1/weight)
    
    #diag(Wst)<-0
    #Wst<-mat2listw(Wst, style="W")
    sample[[paste0('lagu',u,'l',l, 'sp',p,'c',c, 'tp',q,sep="")]]<-Wst %*% sample$price
    sample[[paste0('lnlagu',u,'l',l, 'sp',p,'c',c, 'tp',q,sep="")]]<-Wst %*%sample$logprice
    
    
    saveRDS(sample, file = path, ascii = FALSE, version = NULL,
            compress = TRUE, refhook = NULL)
    #stop cluster
    
  }
  start.time <- Sys.time()
  
  W.trend.lag.variables(upper.spatial.range,lower.spatial.range,spatial.power.range,
                        temporal.cut.range,temporal.power.range,paste0(path,"pretreatlag",i".rds"))
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken
  
  
}

      
      
      
      
      
      
      
      
      
      
      