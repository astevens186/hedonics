
## Preliminaries
rm(list=ls())
gc()

# Change working directory to where you've stored ZTRAX
path<- "P:/Peter/Hedonics/Time/"
path2<- "P:/Peter/Hedonics/"
path<- "P:/ansteve2/Peter/Hedonics/Time/"
path2<- "P:/ansteve2/Peter/Hedonics/"


## This function will check if a package is installed, and if not, install it
pkgTest <- function(x) {
  if (!require(x, character.only = TRUE))
  {
    install.packages(x, dep = TRUE)
    if(!require(x, character.only = TRUE)) stop("Package not found")
  }
}

## These lines load the required packages
packages <- c("readxl","rgdal","matrixStats","matrixStats","MatchIt","Matching","foreign","multiwayvcov","lmtest","readstata13","xlsx", "data.table","doSNOW","parallel","compare","doParallel","devtools","foreach","spdep","reshape2","sm","plyr","utils","tcltk","geosphere", "matrixcalc", "dplyr","ExPosition", "randomForest","lfe", "hdm", "rdrobust", "stargazer", "ggplot2", "outliers","rpart","e1071")
lapply(packages, pkgTest)


## These lines set several options
options(scipen = 999) # Do not print scientific notation
options(stringsAsFactors = FALSE) ## Do not load strings as factors
options(max.print=1000000)
memory.limit(10000000000000)

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

#only take first listings
#NPL <- subset(NPL, rat_act_code!="NP002")
#NPL <- subset(NPL, rat_act_code!="NF002")

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



#NPL<-NPL[!is.na(NPL$gw),]
## save
saveRDS(NPL, file = paste(path,'npl1.rds', sep=""), ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
NPL<-readRDS(paste(path,'npl1.rds', sep=""), refhook = NULL)

setwd("//141.142.208.117/share")
##ADD CENSUS DATA
census <- read.csv("projects/HUD_Discrimination/stores/Census_RSEI/IPF_race_ethnic_acs2012_geo2010.csv")
census$blockgroup<-census$blkgrp
census$blkgrp <- as.character(sprintf('%012.0f', census$blkgrp)) 
census$tot2012 <- as.numeric(census$w2012+census$b2012+census$a2012+census$hisp2012+census$oth2012)
census$w2012pc <- as.numeric(census$w2012/(census$tot2012))
census$b2012pc <- as.numeric(census$b2012/(census$tot2012))
census$a2012pc <- as.numeric(census$a2012/(census$tot2012))
census$hisp2012pc <- as.numeric(census$hisp2012/(census$tot2012))
census$oth2012pc <- as.numeric(census$oth2012/(census$tot2012))
acs <- read.csv("projects/HUD_Discrimination/stores/Census_RSEI/nhgis0022_ds191_20125_2012_blck_grp.csv")
acs$share40 <- as.numeric((acs$QU0E002 + acs$QU0E003 + acs$QU0E004 + acs$QU0E005 + acs$QU0E006 + acs$QU0E007 + acs$QU0E008)/acs$QU0E001)
acs$STATEa <- sprintf('%02d', acs$STATEA)
acs$COUNTYa <- sprintf('%03d', acs$COUNTYA)
acs$TRACTa <- sprintf('%06d', acs$TRACTA)
acs$blkgrp <- paste0(acs$STATEa, acs$COUNTYa, acs$TRACTa, acs$BLKGRPA)
acsmerge <- acs[c("blkgrp", "share40")]
mergeacs <- merge(census,acsmerge, by = c("blkgrp"))
acs2 <- read.csv("projects/HUD_Discrimination/stores/nhgis0001_csv/nhgis0001_ds191_20125_2012_blck_grp.csv")
acs2$workcity <- as.numeric((acs2$QS9E004 + acs2$QS9E007 + acs2$QS9E010 + acs2$QS9E015 + acs2$QS9E018 + acs2$QS9E021)/acs2$QS9E001)
acs2$pubtrans <- as.numeric((acs2$QTFE010 + acs2$QTFE011 + acs2$QTFE012 + acs2$QTFE013 + acs2$QTFE014 + acs2$QTFE015)/acs2$QTFE001)
acs2$walkbike <- as.numeric((acs2$QTFE018 + acs2$QTFE019)/acs2$QTFE001)
acs2$commuteless20 <- as.numeric((acs2$QTHE002 + acs2$QTHE003 + acs2$QTHE004 + acs2$QTHE005)/acs2$QTHE001)
acs2$commute <- as.numeric((acs2$QTHE002 + acs2$QTHE003 + acs2$QTHE004 + acs2$QTHE005 + acs2$QTHE006 + acs2$QTHE007 + acs2$QTHE008 + acs2$QTHE009)/acs2$QTHE001)
acs2$commute60 <- as.numeric((acs2$QTHE012 + acs2$QTHE013)/acs2$QTHE001)
acs2$college <- as.numeric((acs2$QUSE021 + acs2$QUSE022 + acs2$QUSE023 + acs2$QUSE024 + acs2$QUSE025)/acs2$QUSE001)
acs2$unemployment <- as.numeric((acs2$QXSE005)/acs2$QXSE003)
acs2$skill <- as.numeric((acs2$QXTE003+acs2$QXTE039)/acs2$QXTE001)
acs2$prod <- as.numeric((acs2$QXTE030+acs2$QXTE034+acs2$QXTE066+acs2$QXTE070)/acs2$QXTE001)
acs2$skillfemale <- as.numeric((acs2$QXTE039)/acs2$QXTE038)
acs2$vacancy <- as.numeric((acs2$QX7E003)/acs2$QX7E001)
acs2$ownocc <- as.numeric((acs2$QX8E002)/acs2$QX8E001)
acs2$povrate <- as.numeric((acs2$QUVE002 + acs2$QUVE003)/acs2$QUVE001)
acs2$STATEa <- sprintf('%02d', acs2$STATEA)
acs2$COUNTYa <- sprintf('%03d', acs2$COUNTYA)
acs2$TRACTa <- sprintf('%06d', acs2$TRACTA)
acs2$blkgrp <- paste0(acs2$STATEa, acs2$COUNTYa, acs2$TRACTa, acs2$BLKGRPA)
acs2merge <- acs2[c("blkgrp", "workcity", "pubtrans", "commute", "college", "unemployment", "vacancy", "ownocc", "commute60", "skill", "walkbike", "commuteless20", "povrate","skillfemale","prod")]
mergeacs2 <- merge(mergeacs,acs2merge, by = c("blkgrp"))
acs3 <- read.csv("projects/HUD_Discrimination/stores/nhgis0003_csv/nhgis0003_ds195_20095_2009_blck_grp.csv")
acs4 <- read.csv("projects/HUD_Discrimination/stores/nhgis0003_csv/nhgis0003_ds225_20165_2016_blck_grp.csv")
acs3$blkgrp <- paste0(sprintf('%02d', acs3$STATEA), sprintf('%03d', acs3$COUNTYA), sprintf('%06d', acs3$TRACTA), acs3$BLKGRPA)
acs4$blkgrp <- paste0(sprintf('%02d', acs4$STATEA), sprintf('%03d', acs4$COUNTYA), sprintf('%06d', acs4$TRACTA), acs4$BLKGRPA)
medinc <- merge(acs3,acs4, by = c("blkgrp"))
d <- 1.10789
medinc$deltainc <- as.numeric((medinc$AF49E001/d) - medinc$RNHE001)
RACEINC <- read.csv("projects/HUD_Discrimination/stores/Census_RSEI/snipped_master_dataset_2012_geo2010.csv")
RACEINC$blkgrp <- RACEINC$fixed_blockg
RACEINC$blkgrp <- as.character(sprintf('%012.0f', RACEINC$blkgrp)) 
RSEI <- read.csv("projects/HUD_Discrimination/stores/Census_RSEI/toxic_2012_2010_blockgroup_consistent2000_2012.csv")
RSEI$blkgrp <- RSEI$block
RSEI$RSEI <- RSEI$concentration
RSEI$blkgrp <- as.character(sprintf('%012.0f', RSEI$blkgrp)) 
RACEINC_RSEI <- merge(RACEINC,RSEI, by = c("blkgrp"))
census <- merge(mergeacs2,RACEINC_RSEI, by = c("blkgrp"))

#cblock <- readOGR("projects/HUD_Discrimination/stores/GriddedRSEI","con_us_810m_1")
#cblock <- spTransform(cblock,"+proj=longlat +datum=WGS84")

cs<-c("01","02","04","05","06","08","09",10:13,15:42,44:51,53:56,72)
#cs<-c("01","02","04")
i<-"02"
NPL<-NPL[!is.na(NPL$lat),]
cbtid<-NULL
cbtid$site_id<-NPL$site_id
cbtid<-data.frame(cbtid)
cbtid$mcbg<-NA
for(i in cs){
cbg <- readOGR(paste0(path2,"Superfund/Data/gz_2010_",i,"_150_00_500k"),paste0("gz_2010_",i,"_150_00_500k"))
cbg <- spTransform(cbg,"+proj=longlat +datum=WGS84")
cbg$bg <-as.numeric(substr(cbg$GEO_ID, 10, 22))

lmat <- NPL[,c("long","lat")]
spdata <- SpatialPointsDataFrame(lmat,lmat)
proj4string(spdata) <- "+proj=longlat +datum=WGS84"
matchedcbg<- over(spdata,cbg)$bg
cbtid$mcbg<-ifelse(!is.na(cbtid$mcbg),cbtid$mcbg,ifelse(!is.na(matchedcbg) ,matchedcbg,NA))
#mcbg<-matchedcbg[!is.na(matchedcbg)]
#cbtid<-merge(cbtid,mcbg, by="")
}
NPL1<-merge(NPL,data.frame(cbtid),by="site_id")
census$nblkgrp<-as.numeric(census$blkgrp)
NPL<-merge(NPL1,census,by.x="mcbg",by.y="nblkgrp")

# Preliminary results
"year"
NPL$hiinc<-NPL$ir_p_asian_hiincome+NPL$ir_p_hispanic_hiincome+NPL$ir_p_other_hiincome+NPL$ir_p_white_hiincome+NPL$ir_p_black_hiincome
NPL$midinc<-NPL$ir_p_asian_midincome+NPL$ir_p_hispanic_midincome+NPL$ir_p_other_midincome+NPL$ir_p_white_midincome+NPL$ir_p_black_midincome
NPL$lowinc<-NPL$ir_p_asian_lowincome+NPL$ir_p_hispanic_lowincome+NPL$ir_p_other_lowincome+NPL$ir_p_white_lowincome+NPL$ir_p_black_lowincome

NPL$noincdata<-ifelse(NPL$hiinc+NPL$midinc+NPL$lowinc<.99,1,0)
var<-c(#"ir_p_asian_midincome","ir_p_asian_hiincome","ir_p_asian_lowincome","ir_p_black_hiincome","ir_p_black_lowincome",
       #"ir_p_black_midincome","ir_p_hispanic_hiincome","ir_p_hispanic_lowincome",
       #"ir_p_hispanic_midincome","ir_p_other_hiincome","ir_p_other_lowincome",
       #"ir_p_other_midincome","ir_p_white_hiincome","ir_p_white_lowincome","ir_p_white_midincome",
       #"w2012","b2012","a2012","hisp2012","oth2012","tot2012",
       "hiinc","lowinc","noincdata",
       #"w2012pc",
       "b2012pc","a2012pc","hisp2012pc","oth2012pc",
       "workcity" ,"pubtrans","commute","college","unemployment","vacancy","ownocc","commute60","skill",
       "walkbike" ,"commuteless20","povrate" ,"skillfemale","prod","sw","gw","a","Coast Guard",
       "Community Organization","EPA Fund-Financed","EPA In-House","Federal Enforcement","Federal Facilities",
       "Mixed Funding Federal/RP","Other","Prospective Purchaser","PRP Lead Under State",
       "PRP Response Under State","Responsible Party","Special Account Financed Action - EPA" ,
       "Special Account Financed Action - State", "State Deferral","State Enforcement","State, Fund Financed",
       "State, No Fund Money","State, with EPA Concurrence","State, without EPA Concurrence","Tribal Lead, Fund-Financed")

varchar<-c(#"ir_p_asian_midincome","ir_p_asian_hiincome","ir_p_asian_lowincome","ir_p_black_hiincome","ir_p_black_lowincome",
  #"ir_p_black_midincome","ir_p_hispanic_hiincome","ir_p_hispanic_lowincome",
  #"ir_p_hispanic_midincome","ir_p_other_hiincome","ir_p_other_lowincome",
  #"ir_p_other_midincome","ir_p_white_hiincome","ir_p_white_lowincome","ir_p_white_midincome",
  #"w2012","b2012","a2012","hisp2012","oth2012","tot2012",
  "hiinc","lowinc","noincdata",
  #"w2012pc",
  "b2012pc","a2012pc","hisp2012pc","oth2012pc",
  "workcity" ,"pubtrans","commute","college","unemployment","vacancy","ownocc","commute60","skill",
  "walkbike" ,"commuteless20","povrate" ,"skillfemale","prod")
  

as.factor(year(NPL$proposalDate))
#+as.factor(year(NPL$finalDate))
fitcc<-lm(NPL$ytocc~as.matrix(NPL[,var]))
summary(fitcc)
sumtabcc<-summary(fitcc)
sigcc<-sumtabcc$coefficients[sumtabcc$coefficients[,4]<0.1,]
sigcc[1:3,]
p<-sigcc[,4]
mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ", ifelse(p < .1, "^\\bullet  ", " "))))

#pb<-exp(get(paste0(statchange,'betas.',meth,'.t.',inf,'.',treat)))-1
#rpb<-round(pb,3)
#se<-round(exp(get(paste0(statchange,'ses.',meth,'.t.',inf,'.',treat)))-1,3)

pb<-sigcc[,1]
rpb<-as.matrix(round(pb,3))
se<-as.matrix(round(sigcc[,2],3))

srpb <- matrix(paste(rpb, mystars, sep=""), ncol=dim(rpb)[2] )
rownames(srpb)<-rownames(rpb)
#nsrpb<-rbind(c("",laglead),cbind(dist,srpb))


results.mat<-matrix(nrow= 2*dim(srpb)[1],ncol= dim(srpb)[2])
rownames(results.mat)<-rep
rerownames(results.mat, do.NULL = FALSE, prefix = " ")

for(i in 1:dim(results.mat)[1]){
  if(i %% 2 != 0){
    results.mat[i,]<-srpb[ceiling(i/2),]
    rownames(results.mat)[i]<-rownames(srpb)[ceiling(i/2)]
  }
  if(i %% 2 == 0){
    results.mat[i,]<-paste0('(',se[ceiling(i/2),],')')
  }
}
rownames(results.mat)<-c(rownames(sigcc)," ")
collapse(paste0(rownames(sigcc)," "))
results.mat<-rbind(c(rownames(sigcc),""),results.mat)

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

fitdel<-lm(NPL$ytodel~as.matrix(NPL[,var])+as.factor(year(NPL$ccDate)))
summary(fitdel)
sumtabd<-summary(fitdel)
sigd<-sumtabd$coefficients[sumtabd$coefficients[,4]<0.1,]
sigd[1:4,]

fitcc<-lm(NPL$ytocc[year(NPL$finalDate)>1994]~as.matrix(NPL[year(NPL$finalDate)>1994,var]))
summary(fitcc)
sumtabcc<-summary(fitcc)
sigcc<-sumtabcc$coefficients[sumtabcc$coefficients[,4]<0.1,]
sigcc[1:3,]

fitcc<-lm(NPL$ytodel[year(NPL$ccDate)>1994]~as.matrix(NPL[year(NPL$ccDate)>1994,var]))
summary(fitcc)
sumtabcc<-summary(fitcc)
sigcc<-sumtabcc$coefficients[sumtabcc$coefficients[,4]<0.1,]
sigcc[1:3,]



ggplot() +
  geom_point(data=NPL, aes(x=finalDate, y=ytocc, color= "Final")) +
  stat_smooth(method = 'loess', formula = y ~ x  ,data=NPL,se= TRUE,  aes(x=finalDate, y=ytocc, color= "Final")) +
  #ggtitle("Common Trends Assumption")+
  xlab('Final Listing') + #labs(color="Legend") +  geom_vline(xintercept=0)+
  ylab('Years from Final Listing to Construction Completion')

  ggsave(file=paste(path,'latex/cctime.png', sep=""),height = 6,width =10)

  ggplot() +
    geom_point(data=NPL, aes(x=ccDate, y=ytodel, color= "Final")) +
    stat_smooth(method = 'loess', formula = y ~ x  ,data=NPL,se= TRUE,  aes(x=ccDate, y=ytodel, color= "Final")) +
    #ggtitle("Common Trends Assumption")+
    xlab('Construction Completion') + #labs(color="Legend") +  geom_vline(xintercept=0)+
    ylab('Years from Construction Completion to Deletion')
  
  ggsave(file=paste(path,'latex/deltime.png', sep=""),height = 6,width =10)
  

  ggplot(data=NPL, aes(year(NPL$finalDate))) + geom_histogram()+xlab("Final NPL")+ylab("Count")
  ggsave(file=paste(path,'latex/countfinal.png', sep=""),height = 6,width =8)
  
  ggplot(data=NPL, aes(year(NPL$ccDate))) + geom_histogram()+xlab("Construction Complete")+ylab("Count")
  ggsave(file=paste(path,'latex/countcc.png', sep=""),height = 6,width =8)
  
  ggplot(data=NPL, aes(year(NPL$deleteDate))) + geom_histogram()+xlab("Deletion from NPL")+ylab("Count")
  ggsave(file=paste(path,'latex/countdel.png', sep=""),height = 6,width =8)
 
  
##############################################################################################  
  #final date 
for(i in varchar){
  #i<-varchar[1]
  fit<-lm(NPL[[paste0(i)]]~as.factor(year(NPL$finalDate)))
  summary(fit)
  
  years<-levels(as.factor(year(NPL$finalDate)))[-1]
  allModelFrame <- data.frame(Variable = years,
                              Coefficient = coef(summary(fit))[,"Estimate"][-1],
                              SE = coef(summary(fit))[,"Std. Error"][-1],
                              modelName = 'Year')
  
  
  interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier
  #allModelFrame<-allModelFrame[allModelFrame$modelName!="2k",]
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
  
  zp1 <- zp1 + scale_x_discrete(breaks = seq(1985,2015, by=5))
  zp1 <- zp1  +  labs(color="Fixed Effects")
  print(zp1)  # The trick to these is position_dodge().
  
  ggsave(file=paste(path,'latex/characterchange',i, 'final.png', sep=""),height = 7,width =9)
  
  
}
  
 #cc date
  
for(i in varchar){
    #i<-varchar[1]
    #i<-3
    fit<-lm(NPL[[paste0(i)]]~as.factor(year(NPL$ccDate)))
    summary(fit)
    
    years<-levels(as.factor(year(NPL$ccDate)))[-1]
    allModelFrame <- data.frame(Variable = years,
                                Coefficient = coef(summary(fit))[,"Estimate"][-1],
                                SE = coef(summary(fit))[,"Std. Error"][-1],
                                modelName = 'Year')
    
    
    interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier
    #allModelFrame<-allModelFrame[allModelFrame$modelName!="2k",]
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
    
    zp1 <- zp1 + scale_x_discrete(breaks = seq(1985,2015, by=5))
    zp1 <- zp1  +  labs(color="Fixed Effects")
    print(zp1)  # The trick to these is position_dodge().
    
    ggsave(file=paste(path,'latex/characterchange',i, 'cc.png', sep=""),height = 7,width =9)
    
    
  }
  
  #delete date
  
for(i in varchar){
    #i<-varchar[1]
    fit<-lm(NPL[[paste0(i)]]~as.factor(year(NPL$deleteDate)))
    summary(fit)
    
    years<-levels(as.factor(year(NPL$deleteDate)))[-1]
    allModelFrame <- data.frame(Variable = years,
                                Coefficient = coef(summary(fit))[,"Estimate"][-1],
                                SE = coef(summary(fit))[,"Std. Error"][-1],
                                modelName = 'Year')
    
    
    interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier
    #allModelFrame<-allModelFrame[allModelFrame$modelName!="2k",]
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
    
    zp1 <- zp1 + scale_x_discrete(breaks = seq(1985,2015, by=5))
    zp1 <- zp1  +  labs(color="Fixed Effects")
    print(zp1)  # The trick to these is position_dodge().
    
    ggsave(file=paste(path,'latex/characterchange',i, 'del.png', sep=""),height = 7,width =9)
    
    
  }
  