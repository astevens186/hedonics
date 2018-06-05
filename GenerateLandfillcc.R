
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
sample<-sample[sample$tgdel|sample$cdel,]
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

sample1<-sample[dist_to_site_min<10000,]
sample1<-sample1[]
## save
saveRDS(sample, file = paste(path,'repeatpr5.rds', sep=""), ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
sample<-readRDS(paste(path,'repeatpr5.rds', sep=""), refhook = NULL)
NPL<-readRDS(paste(path,'NPLfull.rds', sep=""), refhook = NULL)

NPL<-NPL[NPL$rstate_code=="NY",]
#NPLct<-NPL[NPL$rstate_code=="CT",]
#NPLvt<-NPL[NPL$rstate_code=="VT",]
#NPLma<-NPL[NPL$rstate_code=="MA",]
#NPLpa<-NPL[NPL$rstate_code=="PA",]

#NPL<-rbind(NPLny,NPLct,NPLvt,NPLma,NPLpa)
pNPL<-subset(NPL,rat_name=="PROPOSAL TO NPL")
fNPL<-subset(NPL,rat_name=="FINAL LISTING ON NPL")
dNPL<-subset(NPL,rat_name=="DELETION FROM NPL"|rat_name=="PARTIAL NPL DELETION")
pdNPL<-subset(NPL,rat_name=="PARTIAL NPL DELETION")

opNPL<- pNPL#[order(pNPL$date),] 
ofNPL<- fNPL#[order(fNPL$date),] 
odNPL<- dNPL#[order(dNPL$date),] 
opdNPL<- dNPL#[order(pdNPL$date),] 

deletedCompletion<-c("ControlsComplete")#,
#"GroundwaterEngComplete","InstControlsComplete","GroundwaterComplete","EngControlsComplete")

finalInstitutional<-c("HealthAdvisory","GroundwaterUseRegulation","DeedNotices","DeedRestriction",
                      "NoticestoStateRegulators",
                      "DrinkingWaterAdvisory","WaterSupplyUseRestriction","AccessRestriction")
deletedContaminant<-c("gw","a","n")

d.sample.data<-sample

#myvars <- names(d.sample.data) %in% c("Intercept") 
#d.sample.data <- d.sample.data[!myvars]

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


#PropertyAddressCensusTractAndBlock.f. = factor(d.sample.data$PropertyAddressCensusTractAndBlock)
#PropertyAddressCensusTractAndBlock.dummies = model.matrix(~PropertyAddressCensusTractAndBlock.f.)
#d.sample.data<- cbind(d.sample.data,PropertyAddressCensusTractAndBlock.dummies)

saveRDS(d.sample.data, file = paste(path,'repeat5NYwells.rds', sep=""), ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

#colMeans(sample[,2667:2725])
#odNPL[odNPL$date>"2013-01-01",]
sample<-readRDS(paste(path,'repeat5NYwells.rds', sep=""), refhook = NULL)
sampleout<-sample

sample$day<-as.numeric(sample$date)
fmmatdf<-as.formula(~TransId+ date + treatmentgroup+ YearBuilt + NoOfStories + sqfeet)
fmmatch<-as.formula(treatmentgroup ~ date+ YearBuilt + NoOfStories+  sqfeet)
#psite<-c("209","210","217","228")
#for(j in deletedContaminant){ 
j<-"gw"
for(i in 1:dim(odNPL)[1]){
  #i<-9
  sample$treatmentgroup<-sampleout[[paste0('treatgroupd',j,i)]] 
  sample$treatdgw<-sampleout[[paste0('treatd',j,i)]] 
  #sample[[paste0('treatmentgroup',j)]][is.na(sample[[paste0('treatmentgroup',j)]])]<-0
  sample$control<-ifelse(sampleout[[paste0('control',i)]] ==1 & sampleout[[paste0('dNPL',i)]]<100000 & sampleout[[paste0('dNPL',i)]]>30000,1,0)
  if(mean(sample$treatmentgroup)>0 & mean(sample$control)>0 &mean(sample[[paste0('treatd',j,i)]])>0){
    sample1<-subset(sample, treatmentgroup==1 | control ==1)
    
    
    saveRDS(sample1, file = paste(path,'fulldeletionbaj',j,i,'.rds', sep=""), ascii = FALSE, version = NULL,
            compress = TRUE, refhook = NULL)
    
    
  }
}
#}






#############################################################################
#potential sites


expandingList <- function(capacity = 10) {
  buffer <- vector('list', capacity)
  length <- 0
  
  methods <- list()
  
  methods$double.size <- function() {
    buffer <<- c(buffer, vector('list', capacity))
    capacity <<- capacity * 2
  }
  
  methods$add <- function(val) {
    if(length == capacity) {
      methods$double.size()
    }
    
    length <<- length + 1
    buffer[[length]] <<- val
  }
  
  methods$as.list <- function() {
    b <- buffer[0:length]
    return(b)
  }
  
  methods
}

psites<-expandingList()
sample<-NULL
for(k in c("full")){
  for(i in 1:dim(odNPL)[1]){
    #k<-"full"
    #i<-73
    if (file.exists(paste(path,k,'deletionbajgw',i,'.rds', sep=""))){
      sample.1<-readRDS(paste(path,k,'deletionbajgw',i,'.rds', sep=""), refhook = NULL)
      print(paste0(k,"and",i,"sample size = ", dim(sample.1)))
      print(paste0(k,"and",i,"treatment group = ", dim(sample.1[sample.1$treatmentgroup>0,])[1]))
      print(paste0(k,"and",i,"treated = ", dim(sample.1[sample.1$treatdgw>0,])[1]))
      print(paste0(k,"and",i,"control = ", dim(sample.1[sample.1$control>0,])[1]))
      sample.1$treatgwWL<- sample.1$treatdgw * sample.1$WaterStndCode.fWL
      print(paste0(k,"and",i,"treat well = ", dim(sample.1[sample.1$treatgwWL>0,])[1]))
      sample.1$treatgwMU<- sample.1$treatdgw * sample.1$WaterStndCode.fMU
      print(paste0(k,"and",i,"treat public = ", dim(sample.1[sample.1$treatgwMU>0,])[1]))
      cut<-500
      if(dim(sample.1[sample.1$treatgwWL>0,])[1]-cut>0 & dim(sample.1[sample.1$treatgwMU>0,])[1]-cut>0 &
         dim(sample.1[sample.1$control>0,])[1]-cut>0 & 
         dim(sample.1)[1]-dim(sample.1[sample.1$control>0,])[1]-dim(sample.1[sample.1$treatdgw>0,])[1]-cut>0){
        psites$add(paste0("site",i))
        #sample<-rbind(sample,sample.1)
      }
      
      
    }
  }
}
psites$as.list()

saveRDS(psites$as.list(), file = paste(path,'repeat5wellslist.rds', sep=""), ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

psitelist<-readRDS(paste(path,'repeat5wellslist.rds', sep=""), refhook = NULL)
###############################################################################################






for(j in 1:dim(odNPL)[1]){
  sample$treatmentgroup<-sampleout[[paste0('treatgroupdgw',j)]]  
  #sample[[paste0('treatmentgroup',j)]][is.na(sample[[paste0('treatmentgroup',j)]])]<-0
  sample$control<-ifelse(sampleout$control==1 & sampleout[[paste0('dNPL',j)]]<100000 & sampleout[[paste0('dNPL',j)]]>20000,1,0)
  if(mean(sample$treatmentgroup)>0 & mean(sample$control)>0){
    sample1<-subset(sample, treatmentgroup==1 | control ==1)
    #sample1$treatmentgroup<-sample1[[paste0('treatmentgroup',j)]]
    
    sample1[[paste0('timeFE',j)]]<-ifelse(sample1$date-odNPL$date[j]>0 & odNPL$gw[j]==1,1,0)
    sample.df<-data.frame(model.matrix(fmmatdf,sample1))
    
    saveRDS(sample1, file = paste(path,'fulldeletion',j,'.rds', sep=""), ascii = FALSE, version = NULL,
            compress = TRUE, refhook = NULL)
    
    sample.m.df<-subset(sample.df, as.numeric(date-min(sample.df$date[sample.df$treatmentgroup==1]))>0)
    if(FALSE){
      mdm<- matchit(fmmatch,data = sample.m.df, method = "nearest", distance = "mahalanobis")
      
      mdm.treat <- match.data(mdm, group = "treat")
      mdm.control <- match.data(mdm, group = "control")
      mdm<-rbind(mdm.treat,mdm.control)
      mdm<-mdm[,c("TransId","date")]
      mdm.full<-merge(mdm,sample1,all.x = TRUE, by = "TransId")
      
      saveRDS(mdm.full , file = paste(path,'mdmdeletion',j,'.rds', sep=""), ascii = FALSE, version = NULL,
              compress = TRUE, refhook = NULL)
      
      lm<-lm(fmmatch,sample.m.df)
      
      psm<- matchit(treatmentgroup~ lm$fitted.values,data = sample.m.df, method = "nearest", distance = "mahalanobis")
      
      psm.treat <- match.data(psm, group = "treat")
      psm.control <- match.data(psm, group = "control")
      psm<-rbind(psm.treat,psm.control)
      psm<-psm[,c("TransId","date")]
      psm.full<-merge(psm,sample1,all.x = TRUE, by = "TransId")
      
      saveRDS(psm.full , file = paste(path,'psmdeletion',j,'.rds', sep=""), ascii = FALSE, version = NULL,
              compress = TRUE, refhook = NULL)
    }
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
      if(FALSE){
        mdm<- matchit(fmmatch,  data = sample.m.df, method = "nearest", distance = "mahalanobis")
        
        mdm.treat <- match.data(mdm, group = "treat")
        mdm.control <- match.data(mdm, group = "control")
        mdm<-rbind(mdm.treat,mdm.control)
        mdm<-mdm[,c("TransId","date")]
        
        mdm.full<-merge(mdm,sample1,all.x = TRUE, by = "TransId")
        
        saveRDS(mdm.full , file = paste(path,'mdm',j,i,'.rds', sep=""), ascii = FALSE, version = NULL,
                compress = TRUE, refhook = NULL)
        lm<-lm(fmmatch,sample.m.df)
        
        psm<- matchit(treatmentgroup~ lm$fitted.values,data = sample.m.df, method = "nearest", distance = "mahalanobis")
        
        psm.treat <- match.data(psm, group = "treat")
        psm.control <- match.data(psm, group = "control")
        psm<-rbind(psm.treat,psm.control)
        psm<-psm[,c("TransId","date")]
        psm.full<-merge(psm,sample1,all.x = TRUE, by = "TransId")
        
        saveRDS(psm.full , file = paste(path,'psm',j,i,'.rds', sep=""), ascii = FALSE, version = NULL,
                compress = TRUE, refhook = NULL)
      }
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
      if(FALSE){
        mdm<- matchit(fmmatch,data = sample.m.df, method = "nearest", distance = "mahalanobis")
        
        mdm.treat <- match.data(mdm, group = "treat")
        mdm.control <- match.data(mdm, group = "control")
        mdm<-rbind(mdm.treat,mdm.control)
        mdm<-mdm[,c("TransId","date")]
        
        mdm.full<-merge(mdm,sample1,all.x = TRUE, by = "TransId")
        
        saveRDS(mdm.full , file = paste(path,'mdm',j,i,'.rds', sep=""), ascii = FALSE, version = NULL,
                compress = TRUE, refhook = NULL)
        
        lm<-lm(fmmatch,sample.m.df)
        
        psm<- matchit(treatmentgroup~ lm$fitted.values,data = sample.m.df, method = "nearest", distance = "mahalanobis")
        
        psm.treat <- match.data(psm, group = "treat")
        psm.control <- match.data(psm, group = "control")
        psm<-rbind(psm.treat,psm.control)
        psm<-psm[,c("TransId","date")]
        psm.full<-merge(psm,sample1,all.x = TRUE, by = "TransId")
        
        saveRDS(psm.full , file = paste(path,'psm',j,i,'.rds', sep=""), ascii = FALSE, version = NULL,
                compress = TRUE, refhook = NULL)
      }
    }
  }
}

psitel<-c(2,4,11,12,15,16,19,20,21)
psitel<-c(2,12,15,16)

psitel<-c(2,15,16)

#psite<-2
for(psite in psitel){
  data<-readRDS(paste(path,'fulldeletionbajgw',psite,'.rds', sep=""), refhook = NULL)
  
  #repeat sales Bajari
  #d.sample.data<-sample1
  #rm(sample1)
  
  gc()
  data<-data[data$price>0,]
  data<-data[!duplicated(data[,c("date","HHID")]),]
  sample.data<-data[,c("date","HHID","TransId","price","logprice")]
  quants<-20
  sample.data$indx<- factor(as.numeric(cut2(as.numeric(sample.data$HHID), g=quants)))
  sample<-NULL
  for(i in 1:quants){
    #i<-1
    d.sample.data<-sample.data[indx==i,]
    
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
    d.sample.data$presstatusd<-ifelse(d.sample.data$predate-as.numeric(odNPL$date[i])>0,1,0)
    d.sample.data$presstatuscc<-ifelse(d.sample.data$predate-as.numeric(odNPL$ControlsComplete[i])>0,1,0)
    
    
    d.sample.data<-d.sample.data[d.sample.data$predate>0,]
    #sample1<-sample1[sample1$presstatus<1 ,]
    #sample1<-sample1[sample1$treatdgw<1 ,]
    
    
    d.sample.data$difflogprice<-d.sample.data$logprice-d.sample.data$prelogprice
    
    sample<-rbind(sample,d.sample.data)
  }
  library(data.table)
  data.dt<-data.table(data)
  sample.dt<-data.table(sample)
  sample.new<-merge(data.dt,sample.dt,all.x=TRUE,by="TransId")
  
  saveRDS(sample.new, file = paste(path,'fullbaj',psite,'.rds', sep=""), ascii = FALSE, version = NULL,
          compress = TRUE, refhook = NULL)
}

