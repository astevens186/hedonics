
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
packages <- c("readxl","xtable","splines","ck37r","data.table","matrixStats","tmle","xgboost", "MatchIt","gtools","statar","foreign","multiwayvcov","lmtest","readstata13","xlsx", "data.table","doSNOW","parallel","compare","doParallel","devtools","foreach","spdep","reshape2","sm","plyr","utils","tcltk","geosphere", "matrixcalc", "dplyr","ExPosition", "randomForest","lfe", "hdm", "rdrobust", "stargazer", "ggplot2", "outliers","rpart","e1071")
lapply(packages, pkgTest)

#library(statar)


## These lines set several options
options(scipen = 999) # Do not print scientific notation
options(stringsAsFactors = FALSE) ## Do not load strings as factors

memory.limit(10000000000000)

NPL<-readRDS(paste(path,'NPLfull.rds', sep=""), refhook = NULL)

pNPL<-subset(NPL,rat_name=="PROPOSAL TO NPL")
fNPL<-subset(NPL,rat_name=="FINAL LISTING ON NPL")
dNPL<-subset(NPL,rat_name=="DELETION FROM NPL"|rat_name=="PARTIAL NPL DELETION")
pdNPL<-subset(NPL,rat_name=="PARTIAL NPL DELETION")

opNPL<- pNPL#[order(pNPL$date),] 
ofNPL<- fNPL#[order(fNPL$date),] 
odNPL<- dNPL#[order(dNPL$date),] 
opdNPL<- pdNPL#[order(pdNPL$date),] 

##set up superlearner for TMLE
num_cores = RhpcBLASctl::get_num_cores()

# How many cores does this computer have?
num_cores

# Use all of those cores for parallel SuperLearner.
options(mc.cores = num_cores)

# Check how many parallel workers we are using: 
getOption("mc.cores")

# Set multicore compatible seed.
set.seed(1, "L'Ecuyer-CMRG")

# multicore superlearner

# different configurations.
xgboost.tune <- list(ntrees = c(50, 100),
                     max_depth = c(5,10,15),
                     shrinkage = c( 0.01,0.1),
                     minobspernode = c(10))

# Set detailed names = T so we can see the configuration for each function.
# Also shorten the name prefix.
xgboost <- create.Learner("SL.xgboost", tune = xgboost.tune, detailed_names = T, name_prefix = "xgb")

# configurations
length(xgboost$names)
xgboost$names

xgboost.tune2 <- list(ntrees = c(50, 100),
                     max_depth = c(5,15),
                     shrinkage = c( 0.01,0.1),
                     minobspernode = c(10))

# Set detailed names = T so we can see the configuration for each function.
# Also shorten the name prefix.
xgboost2 <- create.Learner("SL.xgboost", tune = xgboost.tune2, detailed_names = T, name_prefix = "xgb")

# configurations
length(xgboost2$names)
xgboost2$names

# different configurations.
glmnet.tune <- list(alpha = c(0,.1, .25,.5,.75,.9,1))

# Set detailed names = T so we can see the configuration for each function.
# Also shorten the name prefix.
glmnet <- create.Learner("SL.glmnet", tune = glmnet.tune, detailed_names = T, name_prefix = "glmnet")

# configurations
length(glmnet$names)
glmnet$names

glmnet.tune2 <- list(alpha = c(0,.5,1))

# Set detailed names = T so we can see the configuration for each function.
# Also shorten the name prefix.
glmnet2 <- create.Learner("SL.glmnet", tune = glmnet.tune2, detailed_names = T, name_prefix = "glmnet")

# configurations
length(glmnet2$names)
glmnet2$names

# different configurations.
randomForest.tune <- list(ntree = c(2000))

# Set detailed names = T so we can see the configuration for each function.
# Also shorten the name prefix.
randomForest <- create.Learner("SL.randomForest", tune = randomForest.tune, 
                               detailed_names = T, name_prefix = "randomForest")

# configurations
length(randomForest$names)
randomForest$names

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
SL.library<-expandingList()
PS.library<-expandingList()
SL.library$add("SL.mean")
PS.library$add("SL.mean")
SL.library2<-expandingList()
PS.library2<-expandingList()
SL.library2$add("SL.mean")
PS.library2$add("SL.mean")
SL.library$add("SL.randomForest")
SL.library$add("SL.xgboost")

for(i in 1:length(glmnet$names)){
  SL.library$add(glmnet$names[i])
}    
for(i in 1:length(glmnet$names)){
  PS.library$add(glmnet$names[i])
}  

for(i in 1:length(glmnet2$names)){
  SL.library2$add(glmnet2$names[i])
}    
for(i in 1:length(glmnet2$names)){
  PS.library2$add(glmnet2$names[i])
}

for(i in 1:length(xgboost$names)){
  SL.library$add(c(xgboost$names[i]))
}

for(i in 1:length(xgboost2$names)){
  SL.library2$add(c(xgboost2$names[i]))
}

for(i in 1:length(randomForest$names)){
  SL.library$add(c(randomForest$names[i]))
}

SL.library$as.list()
SL.library2$as.list()
PS.library$as.list()
PS.library2$as.list()


#############################################################################
#potential sites
psites<-expandingList()
sample<-NULL
for(k in c("full")){
  for(i in psite){
    #k<-"full"
    #i<-73
    if (file.exists(paste(path,k,'deletionbajgw',i,'.rds', sep=""))){
      sample.1<-readRDS(paste(path,k,'deletionbajgw',i,'.rds', sep=""), refhook = NULL)
      print(paste0(k,"and",i,"sample size = ", dim(sample.1[sample.1$treatmentgroup>0,])))
      print(paste0(k,"and",i,"treated = ", dim(sample.1[sample.1$treatdgw>0,])))
      print(paste0(k,"and",i,"control = ", dim(sample.1[sample.1$control>0,])))
      sample.1$treatgwWL<- sample.1$treatdgw * sample.1$WaterStndCode.fWL
      print(paste0(k,"and",i,"treat well = ", dim(sample.1[sample.1$treatgwWL>0,])[1]))
      sample.1$treatgwMU<- sample.1$treatdgw * sample.1$WaterStndCode.fMU
      print(paste0(k,"and",i,"treat public = ", dim(sample.1[sample.1$treatgwMU>0,])[1]))
      cut<-100
      if(dim(sample.1[sample.1$treatgwWL>0,])[1]-cut>0 & dim(sample.1[sample.1$treatgwMU>0,])[1]-cut>0 &
        dim(sample.1[sample.1$control>0,])[1]-cut>0 & 
        dim(sample.1[sample.1$treatmentgroup>0,])[1]-dim(sample.1[sample.1$treatdgw>0,])[1]-cut>0 &
        dim(sample.1)[2]-54275>0){
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
psite<-c("209","210","214","217","223","224","227","228","254","262")
#sample.1<-readRDS(paste(path,k,'deletiongw',i,'.rds', sep=""), refhook = NULL)
for(k in c("full")){
  sample.1<-readRDS(paste(path,k,'deletion',psite[1],'.rds', sep=""), refhook = NULL)
  for(i in psite[2:10]){
    if (file.exists(paste(path,k,'deletion',i,'.rds', sep=""))){
      sample.2<-readRDS(paste(path,k,'deletion',i,'.rds', sep=""), refhook = NULL)
      sample.1<-smartbind(sample.1,sample.2)
    }
  }
}
saveRDS(sample.1, file = paste(path,'repeat5wells.rds', sep=""), ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

####################################################################

#sample.1<-readRDS(paste(path,'repeat5wells.rds', sep=""), refhook = NULL)

sample.209<-readRDS(paste(path,'fulldeletiongw',209,'.rds', sep=""), refhook = NULL)
sample.210<-readRDS(paste(path,'fulldeletiongw',210,'.rds', sep=""), refhook = NULL)
sample.217<-readRDS(paste(path,'fulldeletiongw',217,'.rds', sep=""), refhook = NULL)
sample.228<-readRDS(paste(path,'fulldeletiongw',228,'.rds', sep=""), refhook = NULL)

sample.209.con<-sample.209[sample.209$control==1,]
sample.209.n<-rbind(sample.209[sample.209$treatdgw209==1,],sample.209[sample.209$treatmentgroup==1,],sample.209.con[1:dim(sample.209[sample.209$treatmentgroup==1,])[1],])
sample.209.n$timeFEgw<-sample.209.n$timeFEgw209

sample.210.con<-sample.210[sample.210$control==1,]
sample.210.n<-rbind(sample.210[sample.210$treatdgw210==1,],sample.210[sample.210$treatmentgroup==1,],sample.210.con[1:dim(sample.210[sample.210$treatmentgroup==1,])[1],])
sample.210.n$timeFEgw<-sample.210.n$timeFEgw210

sample.217.con<-sample.217[sample.217$control==1,]
sample.217.n<-rbind(sample.217[sample.217$treatdgw217==1,],sample.217[sample.217$treatmentgroup==1,],sample.217.con[1:dim(sample.217[sample.217$treatmentgroup==1,])[1],])
sample.217.n$timeFEgw<-sample.217.n$timeFEgw217

sample.228.con<-sample.228[sample.228$control==1,]
sample.228.n<-rbind(sample.228[sample.228$treatdgw228==1,],sample.228[sample.228$treatmentgroup==1,],sample.228.con[1:dim(sample.228[sample.228$treatmentgroup==1,])[1],])
sample.228.n$timeFEgw<-sample.228.n$timeFEgw228

sample<-rbindlist(list(sample.209.n, sample.210.n, sample.217.n, sample.228.n), fill = TRUE)
#sample<-smartbind(sample.209.n,sample.210.n,sample.217.n,sample.228.n)
#sample<-subset(sample,select=-c(63:321,325:585,610:1191))

saveRDS(sample, file = paste(path,'repeat5wellspres.rds', sep=""), ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

#######################################################################################

d.sample.data<-readRDS(paste(path,'repeat5wellspres.rds', sep=""), refhook = NULL)

#repeat sales Bajari
#d.sample.data<-sample1
#rm(sample1)
gc()
d.sample.data<-d.sample.data[d.sample.data$price>0,]
d.sample.data<-d.sample.data[!duplicated(d.sample.data[,c("date","HHID")]),]

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
A<-rep.row(as.numeric(d.sample.data$HHID),nrow(d.sample.data))
D<-t(A)-A
D[D>0]<-2
D[D<0]<-2
D[D==0]<-1
D[D==2]<-0
sameHouse<-D

A<-rep.row(as.numeric(d.sample.data$TransId),nrow(d.sample.data))
D<-t(A)-A
D[D>0]<-2
D[D<0]<-2
D[D==0]<-1
D[D==2]<-0
sameSale<-D

otherSales<-sameHouse-sameSale
rm(sameHouse,sameSale)
gc()
A<-rep.row(d.sample.data$date,nrow(d.sample.data))
D<-t(A)-A
D[D<0]<-0
diffDates<-D
rm(A,D)
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
d.sample.data$presstatus<-ifelse(d.sample.data$predate-as.numeric(odNPL$date[i])>0,1,0)


sample<-d.sample.data[d.sample.data$predate>0,]
#sample1<-sample1[sample1$presstatus<1 ,]
#sample1<-sample1[sample1$treatdgw<1 ,]

sample<-sample[!is.na(sqfeet),]
sample<-sample[!is.na(YearBuilt),]
sample<-sample[!is.na(LotSizeSquareFeet),]
sample<-sample[!is.na(PropertyAddressLatitude),]    
sample<-sample[!is.na(NoOfStories),] 
sample<-sample[!is.na(TotalRooms),] 
#sample<-subset(sample, TotalRooms>0)
sample<-sample[!is.na(TotalBedrooms),] 
library(dplyr)

sample$difflogprice<-sample$logprice-sample$prelogprice

#Total Average Treatment Effect
timefe<-select(sample, starts_with('timefe'))
treatgroupm<-select(sample, starts_with('treatmentgroup'))
year<-select(sample, starts_with('year'))
bin<-select(sample, starts_with('bin'))
sample$demlogprice<-demeanlist(sample$logprice,
                               list(as.factor(sample$PropertyAddressCensusTractAndBlock)))

#TATE<-as.formula(logprice ~ treatdgw+ treatmentgroup+data.matrix(treatgroupm)+ data.matrix(timefe)+ data.matrix(year[,3:25])+ as.factor(HHID))
feTATE<-model.matrix(~ treatmentgroup+data.matrix(treatgroupm)+ 
                       data.matrix(year[,4:25]),sample)
xTATE<-model.matrix(~ LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                      sqfeet+presstatus+prediffdate+predate+prelogprice,sample)
xcTATE<-model.matrix(~(data.matrix(bin))^3+poly(prediffdate,4)+poly(predate,4)+poly(prelogprice,4)+
                       poly(LotSizeSquareFeet,4) + poly(YearBuilt,4) +   
                       poly(sqfeet,4),sample)

lm.TATE<-felm(logprice ~treatdgw+ LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                sqfeet+ +presstatus+prediffdate+predate+prelogprice+data.matrix(treatgroupm)   + 
              + data.matrix(year[,3:25])|as.factor(PropertyAddressCensusTractAndBlock),sample)
tsample<-sample[control==0,]
year<-select(tsample, starts_with('year'))
t.lm.TATE<-felm(logprice ~treatdgw+ LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                sqfeet+ +presstatus+prediffdate+predate+prelogprice+ treatmentgroup   + 
                + data.matrix(year[,3:25])|as.factor(PropertyAddressCensusTractAndBlock),tsample)

#assign(paste0('lm.TATE.',k,i),lm.TATE)
summary(lm.TATE)
summary(t.lm.TATE)

#Well Average Treatment Effect
sample$treatdgwWL<- sample$treatdgw * sample$WaterStndCode.fWL
sample$treatgroupWL<-sample$treatmentgroup * sample$WaterStndCode.fWL
sample$controlWL<-sample$control*sample$WaterStndCode.fWL

sample$sample.WLATE<-sample$control+sample$treatgroupWL

sample.WLATE<-subset(sample, sample.WLATE==1)
timefe<-select(sample.WLATE, starts_with('timefe'))
#treatgroupm<-select(sample.WLATE, starts_with('treatmentgroup'))
year<-select(sample.WLATE, starts_with('year'))
bin<-select(sample.WLATE, starts_with('bin'))


sample.WLATE$demlogprice<-demeanlist(sample.WLATE$logprice,
                                     list(as.factor(sample.WLATE$PropertyAddressCensusTractAndBlock)))

#WLATE<-as.formula(logprice ~ treatgwWL+ treatgroupWL+data.matrix(treatgroupm)+ data.matrix(timefe)+ data.matrix(year[,3:25])+ as.factor(HHID))
feWLATE<-model.matrix(~ treatmentgroup+ 
                       data.matrix(year[,4:25]),sample.WLATE)
xWLATE<-model.matrix(~ LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                      sqfeet+presstatus+prediffdate+predate+prelogprice,sample.WLATE)
xcWLATE<-model.matrix(~(data.matrix(bin))^3+poly(prediffdate,4)+poly(predate,4)+poly(prelogprice,4)+
                       poly(LotSizeSquareFeet,4) + poly(YearBuilt,4) +   
                       poly(sqfeet,4),sample.WLATE)

lm.WLATE<-felm(logprice ~treatdgwWL+ LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                sqfeet+ treatmentgroup   + 
                + data.matrix(year[,3:25])+presstatus+prediffdate+predate+prelogprice
               |as.factor(PropertyAddressCensusTractAndBlock),sample.WLATE)
tsample.WLATE<-sample.WLATE[control==0,]
year<-select(tsample.WLATE, starts_with('year'))
t.lm.WLATE<-felm(logprice ~treatdgwWL+ LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                 sqfeet+ treatmentgroup   + 
                 + data.matrix(year[,3:25])+presstatus+prediffdate+predate+prelogprice
               |as.factor(PropertyAddressCensusTractAndBlock),tsample.WLATE)
#assign(paste0('lm.TATE.',k,i),lm.TATE)
summary(lm.WLATE)
summary(t.lm.WLATE)

#Municipal Utility Average Treatment Effect
sample$treatdgwMU<- sample$treatdgw * sample$WaterStndCode.fMU
sample$treatgroupMU<-sample$treatmentgroup * sample$WaterStndCode.fMU
sample$controlMU<-sample$control*sample$WaterStndCode.fMU

sample$sample.MUATE<-sample$control+sample$treatgroupMU

sample.MUATE<-subset(sample, sample.MUATE==1)
timefe<-select(sample.MUATE, starts_with('timefe'))
treatgroupm<-select(sample.MUATE, starts_with('treatmentgroup'))
year<-select(sample.MUATE, starts_with('year'))
bin<-select(sample.MUATE, starts_with('bin'))
sample.MUATE$demlogprice<-demeanlist(sample.MUATE$logprice,
                                     list(as.factor(sample.MUATE$PropertyAddressCensusTractAndBlock)))

#MUATE<-as.formula(logprice ~ treatgwMU+ treatgroupMU+data.matrix(treatgroupm)+ data.matrix(timefe)+ data.matrix(year[,3:25])+ as.factor(HHID))
feMUATE<-model.matrix(~ treatmentgroup+
                       data.matrix(year[,4:25]),sample.MUATE)
xMUATE<-model.matrix(~ LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                      sqfeet+presstatus+prediffdate+predate+prelogprice,sample.MUATE)
xcMUATE<-model.matrix(~(data.matrix(bin))^3+poly(prediffdate,4)+poly(predate,4)+poly(prelogprice,4)+
                       poly(LotSizeSquareFeet,4) + poly(YearBuilt,4) +   
                       poly(sqfeet,4),sample.MUATE)

lm.MUATE<-felm(logprice ~treatdgwMU+ LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                sqfeet+treatmentgroup   + 
                + data.matrix(year[,3:25])+presstatus+prediffdate+predate+prelogprice
               |as.factor(PropertyAddressCensusTractAndBlock),sample.MUATE)
tsample.MUATE<-sample.MUATE[control==0,]
year<-select(tsample.MUATE, starts_with('year'))
t.lm.MUATE<-felm(logprice ~treatdgwMU+ LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                 sqfeet+treatmentgroup   + 
                 + data.matrix(year[,3:25])+presstatus+prediffdate+predate+prelogprice
               |as.factor(PropertyAddressCensusTractAndBlock),tsample.MUATE)
#assign(paste0('lm.TATE.',k,i),lm.TATE)
summary(lm.MUATE)
summary(t.lm.MUATE)

olsdidt<-rbind(lm.TATE$beta[1],lm.TATE$pval[1])
olsdidwl<-rbind(lm.WLATE$beta[1],lm.WLATE$pval[1])
olsdidmu<-rbind(lm.MUATE$beta[1],lm.MUATE$pval[1])

olstab<-cbind(olsdidt,olsdidwl,olsdidmu)

print.xtable(xtable(olstab),
             type="latex", file=paste0(path,"olstab.tex"))

t.olsdidt<-rbind(t.lm.TATE$beta[1],t.lm.TATE$pval[1])
t.olsdidwl<-rbind(t.lm.WLATE$beta[1],t.lm.WLATE$pval[1])
t.olsdidmu<-rbind(t.lm.MUATE$beta[1],t.lm.MUATE$pval[1])

t.olstab<-cbind(t.olsdidt,t.olsdidwl,t.olsdidmu)

print.xtable(xtable(t.olstab),
             type="latex", file=paste0(path,"tolstab.tex"))

#######Estimation

xxTATE<-cbind(xcTATE,xTATE,feTATE)
xxWLATE<-cbind(xcWLATE,xWLATE,feWLATE)
xxMUATE<-cbind(xcMUATE,xMUATE,feMUATE)


FxcTATE<-logical(dim(xcTATE)[2])
TxTATE<-!logical(dim(xTATE)[2])
TfeTATE<-!logical(dim(feTATE)[2])
indTATE<-rbind(data.matrix(FxcTATE),data.matrix(TxTATE),data.matrix(TfeTATE))

fit.ds.os1 <- rlassoEffect(x=xxTATE, y=sample$logprice, d=sample$treatdgw,data = sample,
                           method = "double selection", I3=indTATE)
tsample<-sample[control==0,]
t.fit.ds.os1 <- rlassoEffect(x=xxTATE[sample$control==0,], y=tsample$logprice, d=tsample$treatdgw,data = tsample,
                             method = "double selection", I3=indTATE)
#assign(paste0('fit.ds.os1.',k,i),fit.ds.os1)
summary(fit.ds.os1)

FxcWLATE<-logical(dim(xcWLATE)[2])
TxWLATE<-!logical(dim(xWLATE)[2])
TfeWLATE<-!logical(dim(feWLATE)[2])
indWLATE<-rbind(data.matrix(FxcWLATE),data.matrix(TxWLATE),data.matrix(TfeWLATE))

fit.ds.os2 <- rlassoEffect(x=xxWLATE, y=sample.WLATE$logprice, d=sample.WLATE$treatdgwWL,data = sample.WLATE, 
                           method = "double selection", I3=indWLATE)
tsample.WLATE<-sample.WLATE[control==0,]
t.fit.ds.os2 <- rlassoEffect(x=xxWLATE[sample.WLATE$control==0,], y=tsample.WLATE$logprice, d=tsample.WLATE$treatdgwWL,data = tsample.WLATE, 
                             method = "double selection", I3=indWLATE)
#assign(paste0('fit.ds.os2.',k,i),fit.ds.os2)
summary(fit.ds.os2)

FxcMUATE<-logical(dim(xcMUATE)[2])
TxMUATE<-!logical(dim(xMUATE)[2])
TfeMUATE<-!logical(dim(feMUATE)[2])
indMUATE<-rbind(data.matrix(FxcMUATE),data.matrix(TxMUATE),data.matrix(TfeMUATE))


fit.ds.os3 <- rlassoEffect(x=xxMUATE, y=sample.MUATE$logprice, d=sample.MUATE$treatdgwMU, data = sample.MUATE,
                           method = "double selection", I3=indMUATE)
tsample.MUATE<-sample.MUATE[control==0,]
t.fit.ds.os3 <- rlassoEffect(x=xxMUATE[sample.MUATE$control==0,], y=tsample.MUATE$logprice, d=tsample.MUATE$treatdgwMU, data = tsample.MUATE,
                             method = "double selection", I3=indMUATE)
#assign(paste0('fit.ds.os3.',k,i),fit.ds.os3)
summary(fit.ds.os3)

dsdidt<-rbind(fit.ds.os1$alpha,fit.ds.os1$pval)
dsdidwl<-rbind(fit.ds.os2$alpha,fit.ds.os2$pval)
dsdidmu<-rbind(fit.ds.os3$alpha,fit.ds.os3$pval)

dstab<-cbind(dsdidt,dsdidwl,dsdidmu)

print.xtable(xtable(dstab),
             type="latex", file=paste0(path,"dstab.tex"))

t.dsdidt<-rbind(t.fit.ds.os1$alpha,t.fit.ds.os1$pval)
t.dsdidwl<-rbind(t.fit.ds.os2$alpha,t.fit.ds.os2$pval)
t.dsdidmu<-rbind(t.fit.ds.os3$alpha,t.fit.ds.os3$pval)

t.dstab<-cbind(t.dsdidt,t.dsdidwl,t.dsdidmu)

print.xtable(xtable(t.dstab),
             type="latex", file=paste0(path,"tdstab.tex"))

###################################################################
#TMLE
setup_parallel_tmle(parallel = "multicore", max_cores = 3,
                    allow_multinode = T, env = .GlobalEnv)
#simple
resultt<- tmle_parallel(Y=sample$demlogprice,A=sample$treatdgw,W=xxTATE,family="gaussian",
                Q.SL.library = SL.library$as.list(),
                g.SL.library = PS.library$as.list(),
               verbose = TRUE)
summary(resultt)

resultwl<- tmle_parallel(Y=sample.WLATE$demlogprice,A=sample.WLATE$treatdgwWL,W=xxWLATE,
                         family="gaussian",
                         Q.SL.library = SL.library$as.list(),
                         g.SL.library = PS.library$as.list(),
                verbose = TRUE)
summary(resultwl)

resultmu<- tmle_parallel(Y=sample.MUATE$demlogprice,A=sample.MUATE$treatdgwMU,W=xxMUATE,
                         family="gaussian",
                         Q.SL.library = SL.library$as.list(),
                         g.SL.library = PS.library$as.list(),
                verbose = TRUE)
summary(resultmu)


####
#Two-step
csample<-sample[control==1,]
resulttc<- tmle_parallel(Y=csample$demlogprice,A=csample$timeFEgw,W=xxTATE[sample$control==1,],family="gaussian",
                        Q.SL.library = SL.library$as.list(),
                        g.SL.library = PS.library$as.list(),
                        verbose = TRUE)
summary(resulttc)

tsample<-sample[control==0,]
resulttt<- tmle_parallel(Y=tsample$demlogprice,A=tsample$treatdgw,W=xxTATE[sample$control==0,],family="gaussian",
                        Q.SL.library = SL.library$as.list(),
                        g.SL.library = PS.library$as.list(),
                        verbose = TRUE)
summary(resulttt)

teffectt<-resulttt$estimates$ATE$psi-resulttc$estimates$ATE$psi
tp<-mean(tsample$treatdgw)
tn<-dim(tsample)[1]
tsx<-sd(tsample$treatdgw)
cp<-mean(csample$timeFEgw)
cn<-dim(csample)[1]
csx<-sd(csample$timeFEgw)
sres2<-(((tn-2)*resulttt$estimates$ATE$var.psi)+((cn-2)*resulttc$estimates$ATE$var.psi))/((tn-2)+(cn-2))
tt<-teffectt/sqrt(sres2*(1/(tsx^2*(tn-1)))*(1/(csx^2*(cn-1))))

csample<-sample.WLATE[control==1,]
resultwlc<- tmle_parallel(Y=csample$demlogprice,A=csample$timeFEgw,W=xxWLATE[sample.WLATE$control==1,],family="gaussian",
                         Q.SL.library = SL.library$as.list(),
                         g.SL.library = PS.library$as.list(),
                         verbose = TRUE)
summary(resultwlc)

tsample<-sample.WLATE[control==0,]
resultwlt<- tmle_parallel(Y=tsample$demlogprice,A=tsample$treatdgw,W=xxWLATE[sample.WLATE$control==0,],family="gaussian",
                         Q.SL.library = SL.library$as.list(),
                         g.SL.library = PS.library$as.list(),
                         verbose = TRUE)
summary(resultwlt)

teffectwl<-resultwlt$estimates$ATE$psi-resultwlc$estimates$ATE$psi
tp<-mean(tsample$treatdgwWL)
tn<-dim(tsample)[1]
tsx<-sd(tsample$treatdgwWL)
cp<-mean(csample$timeFEgw)
cn<-dim(csample)[1]
csx<-sd(csample$timeFEgw)
sres2<-(((tn-2)*resultwlt$estimates$ATE$var.psi)+((cn-2)*resultwlc$estimates$ATE$var.psi))/((tn-2)+(cn-2))
twl<-teffectwl/sqrt(sres2*(1/(tsx^2*(tn-1)))*(1/(csx^2*(cn-1))))

csample<-sample.MUATE[control==1,]
resultmuc<- tmle_parallel(Y=csample$demlogprice,A=csample$timeFEgw,W=xxTATE[sample$control==1,],family="gaussian",
                         Q.SL.library = SL.library$as.list(),
                         g.SL.library = PS.library$as.list(),
                         verbose = TRUE)
summary(resultmuc)

tsample<-sample.MUATE[control==0,]
resultmut<- tmle_parallel(Y=tsample$demlogprice,A=tsample$treatdgw,W=xxMUATE[sample.MUATE$control==0,],family="gaussian",
                         Q.SL.library = SL.library$as.list(),
                         g.SL.library = PS.library$as.list(),
                         verbose = TRUE)
summary(resultmut)

teffectmu<-resultmut$estimates$ATE$psi-resultmuc$estimates$ATE$psi
tp<-mean(tsample$treatdgwMU)
tn<-dim(tsample)[1]
tsx<-sd(tsample$treatdgwMU)
cp<-mean(csample$timeFEgw)
cn<-dim(csample)[1]
csx<-sd(csample$timeFEgw)
sres2<-(((tn-2)*resultmut$estimates$ATE$var.psi)+((cn-2)*resultmuc$estimates$ATE$var.psi))/((tn-2)+(cn-2))
tmu<-teffectmu/sqrt(sres2*(1/(tsx^2*(tn-1)))*(1/(csx^2*(cn-1))))


tmlet<-rbind(resultt$estimates$ATE$psi,resultt$estimates$ATE$psi/sqrt(resultt$estimates$ATE$var.psi/dim(sample)[1]))
tmlewl<-rbind(resultwl$estimates$ATE$psi,resultwl$estimates$ATE$psi/sqrt(resultwl$estimates$ATE$var.psi/dim(sample)[1]))
tmlemu<-rbind(resultmu$estimates$ATE$psi,resultmu$estimates$ATE$psi/sqrt(resultmu$estimates$ATE$var.psi/dim(sample)[1]))

tmles<-cbind(tmlet,tmlewl,tmlemu)

tmlett<-rbind(resulttt$estimates$ATE$psi,resulttt$estimates$ATE$psi/sqrt(resulttt$estimates$ATE$var.psi/dim(sample)[1]))
tmlewlt<-rbind(resultwlt$estimates$ATE$psi,resultwlt$estimates$ATE$psi/sqrt(resultwlt$estimates$ATE$var.psi/dim(sample)[1]))
tmlemut<-rbind(resultmut$estimates$ATE$psi,resultmut$estimates$ATE$psi/sqrt(resultmut$estimates$ATE$var.psi/dim(sample)[1]))

tmlest<-cbind(tmlett,tmlewlt,tmlemut)

tmledidt<-rbind(teffectt,tt)
tmledidwl<-rbind(teffectwl,twl)
tmledidmu<-rbind(teffectmu,tmu)

tmle<-cbind(tmledidt,tmledidwl,tmledidmu)

print.xtable(xtable(tmle),
             type="latex", file=paste0(path,"logtmle.tex"))

print.xtable(xtable(tmles),
             type="latex", file=paste0(path,"logtmles.tex"))


print.xtable(xtable(tmlest),
             type="latex", file=paste0(path,"logtmlest.tex"))


##################################################################
sdf<-50
lat<-sample$PropertyAddressLatitude
long<-sample$PropertyAddressLongitude
splat<-bs(lat, df = sdf)
splong<-bs(long, df = sdf)
spint<-model.matrix(~splat:splong)

xcTATE<-cbind(xcTATE,splat,splong,spint,lat,long)

lat<-sample.WLATE$PropertyAddressLatitude
long<-sample.WLATE$PropertyAddressLongitude
splat<-bs(lat, df = sdf)
splong<-bs(long, df = sdf)
spint<-model.matrix(~splat:splong)
xcWLATE<-cbind(xcWLATE,splat,splong,spint,lat,long)

lat<-sample.MUATE$PropertyAddressLatitude
long<-sample.MUATE$PropertyAddressLongitude
splat<-bs(lat, df = sdf)
splong<-bs(long, df = sdf)
spint<-model.matrix(~splat:splong)
xcMUATE<-cbind(xcMUATE,splat,splong,spint,lat,long)
#######Estimation

xxTATE<-cbind(xcTATE,xTATE,feTATE)
xxWLATE<-cbind(xcWLATE,xWLATE,feWLATE)
xxMUATE<-cbind(xcMUATE,xMUATE,feMUATE)

FxcTATE<-logical(dim(xcTATE)[2])
TxTATE<-!logical(dim(xTATE)[2])
TfeTATE<-!logical(dim(feTATE)[2])
indTATE<-rbind(data.matrix(FxcTATE),data.matrix(TxTATE),data.matrix(TfeTATE))

fit.ds.os1 <- rlassoEffect(x=xxTATE, y=sample$logprice, d=sample$treatdgw,data = sample,
                           method = "double selection", I3=indTATE)
tsample<-sample[control==0,]
t.fit.ds.os1 <- rlassoEffect(x=xxTATE[sample$control==0,], y=tsample$logprice, d=tsample$treatdgw,data = tsample,
                           method = "double selection", I3=indTATE)
#assign(paste0('fit.ds.os1.',k,i),fit.ds.os1)
summary(fit.ds.os1)

FxcWLATE<-logical(dim(xcWLATE)[2])
TxWLATE<-!logical(dim(xWLATE)[2])
TfeWLATE<-!logical(dim(feWLATE)[2])
indWLATE<-rbind(data.matrix(FxcWLATE),data.matrix(TxWLATE),data.matrix(TfeWLATE))

fit.ds.os2 <- rlassoEffect(x=xxWLATE, y=sample.WLATE$logprice, d=sample.WLATE$treatdgwWL,data = sample.WLATE, 
                           method = "double selection", I3=indWLATE)
tsample.WLATE<-sample.WLATE[control==0,]
t.fit.ds.os2 <- rlassoEffect(x=xxWLATE[sample.WLATE$control==0,], y=tsample.WLATE$logprice, d=tsample.WLATE$treatdgwWL,data = tsample.WLATE, 
                           method = "double selection", I3=indWLATE)
#assign(paste0('fit.ds.os2.',k,i),fit.ds.os2)
summary(fit.ds.os2)

FxcMUATE<-logical(dim(xcMUATE)[2])
TxMUATE<-!logical(dim(xMUATE)[2])
TfeMUATE<-!logical(dim(feMUATE)[2])
indMUATE<-rbind(data.matrix(FxcMUATE),data.matrix(TxMUATE),data.matrix(TfeMUATE))


fit.ds.os3 <- rlassoEffect(x=xxMUATE, y=sample.MUATE$logprice, d=sample.MUATE$treatdgwMU, data = sample.MUATE,
                           method = "double selection", I3=indMUATE)
tsample.MUATE<-sample.MUATE[control==0,]
t.fit.ds.os3 <- rlassoEffect(x=xxMUATE[sample.MUATE$control==0,], y=tsample.MUATE$logprice, d=tsample.MUATE$treatdgwMU, data = tsample.MUATE,
                           method = "double selection", I3=indMUATE)
#assign(paste0('fit.ds.os3.',k,i),fit.ds.os3)
summary(fit.ds.os3)

dsdidt<-rbind(fit.ds.os1$alpha,fit.ds.os1$pval)
dsdidwl<-rbind(fit.ds.os2$alpha,fit.ds.os2$pval)
dsdidmu<-rbind(fit.ds.os3$alpha,fit.ds.os3$pval)

dstab<-cbind(dsdidt,dsdidwl,dsdidmu)

print.xtable(xtable(dstab),
             type="latex", file=paste0(path,"dstabspatial.tex"))

t.dsdidt<-rbind(t.fit.ds.os1$alpha,t.fit.ds.os1$pval)
t.dsdidwl<-rbind(t.fit.ds.os2$alpha,t.fit.ds.os2$pval)
t.dsdidmu<-rbind(t.fit.ds.os3$alpha,t.fit.ds.os3$pval)

t.dstab<-cbind(t.dsdidt,t.dsdidwl,t.dsdidmu)

print.xtable(xtable(t.dstab),
             type="latex", file=paste0(path,"tdstabspatial.tex"))

######################

#TMLE
setup_parallel_tmle(parallel = "multicore", max_cores = 3,
                    allow_multinode = T, env = .GlobalEnv)

####
#Two-step
csample<-sample[control==1,]
resulttc<- tmle_parallel(Y=csample$logprice,A=csample$timeFEgw,W=xxTATE[sample$control==1,],family="gaussian",
                         Q.SL.library = SL.library2$as.list(),
                         g.SL.library = PS.library2$as.list(),
                         verbose = TRUE)
summary(resulttc)

tsample<-sample[control==0,]
resulttt<- tmle_parallel(Y=tsample$logprice,A=tsample$treatdgw,W=xxTATE[sample$control==0,],family="gaussian",
                         Q.SL.library = SL.library2$as.list(),
                         g.SL.library = PS.library2$as.list(),
                         verbose = TRUE)
summary(resulttt)

teffectt<-resulttt$estimates$ATE$psi-resulttc$estimates$ATE$psi
tp<-mean(tsample$treatdgw)
tn<-dim(tsample)[1]
tsx<-sd(tsample$treatdgw)
cp<-mean(csample$timeFEgw)
cn<-dim(csample)[1]
csx<-sd(csample$timeFEgw)
sres2<-(((tn-2)*resulttt$estimates$ATE$var.psi)+((cn-2)*resulttc$estimates$ATE$var.psi))/((tn-2)+(cn-2))
tt<-teffectt/sqrt(sres2*(1/(tsx^2*(tn-1)))*(1/(csx^2*(cn-1))))

csample<-sample.WLATE[control==1,]
resultwlc<- tmle_parallel(Y=csample$logprice,A=csample$timeFEgw,W=xxWLATE[sample.WLATE$control==1,],family="gaussian",
                          Q.SL.library = SL.library2$as.list(),
                          g.SL.library = PS.library2$as.list(),
                          verbose = TRUE)
summary(resultwlc)

tsample<-sample.WLATE[control==0,]
resultwlt<- tmle_parallel(Y=tsample$logprice,A=tsample$treatdgwWL,W=xxWLATE[sample.WLATE$control==0,],family="gaussian",
                          Q.SL.library = SL.library2$as.list(),
                          g.SL.library = PS.library2$as.list(),
                          verbose = TRUE)
summary(resultwlt)

teffectwl<-resultwlt$estimates$ATE$psi-resultwlc$estimates$ATE$psi
tp<-mean(tsample$treatdgwWL)
tn<-dim(tsample)[1]
tsx<-sd(tsample$treatdgwWL)
cp<-mean(csample$timeFEgw)
cn<-dim(csample)[1]
csx<-sd(csample$timeFEgw)
sres2<-(((tn-2)*resultwlt$estimates$ATE$var.psi)+((cn-2)*resultwlc$estimates$ATE$var.psi))/((tn-2)+(cn-2))
twl<-teffectwl/sqrt(sres2*(1/(tsx^2*(tn-1)))*(1/(csx^2*(cn-1))))

csample<-sample.MUATE[control==1,]
resultmuc<- tmle_parallel(Y=csample$logprice,A=csample$timeFEgw,W=xxTATE[sample.MUATE$control==1,],family="gaussian",
                          Q.SL.library = SL.library2$as.list(),
                          g.SL.library = PS.library2$as.list(),
                          verbose = TRUE)
summary(resultmuc)

tsample<-sample.MUATE[control==0,]
resultmut<- tmle_parallel(Y=tsample$logprice,A=tsample$treatdgwMU,W=xxMUATE[sample.MUATE$control==0,],family="gaussian",
                          Q.SL.library = SL.library2$as.list(),
                          g.SL.library = PS.library2$as.list(),
                          verbose = TRUE)
summary(resultmut)

teffectmu<-resultmut$estimates$ATE$psi-resultmuc$estimates$ATE$psi
tp<-mean(tsample$treatdgwMU)
tn<-dim(tsample)[1]
tsx<-sd(tsample$treatdgwMU)
cp<-mean(csample$timeFEgw)
cn<-dim(csample)[1]
csx<-sd(csample$timeFEgw)
sres2<-(((tn-2)*resultmut$estimates$ATE$var.psi)+((cn-2)*resultmuc$estimates$ATE$var.psi))/((tn-2)+(cn-2))
tmu<-teffectmu/sqrt(sres2*(1/(tsx^2*(tn-1)))*(1/(csx^2*(cn-1))))


tmlet<-rbind(resulttt$estimates$ATE$psi,resulttt$estimates$ATE$psi/sqrt(resulttt$estimates$ATE$var.psi/dim(sample)[1]))
tmlewl<-rbind(resultwlt$estimates$ATE$psi,resultwlt$estimates$ATE$psi/sqrt(resultwlt$estimates$ATE$var.psi/dim(sample)[1]))
tmlemu<-rbind(resultmut$estimates$ATE$psi,resultmut$estimates$ATE$psi/sqrt(resultmut$estimates$ATE$var.psi/dim(sample)[1]))

tmles<-cbind(tmlet,tmlewl,tmlemu)
tmledidt<-rbind(teffectt,tt)
tmledidwl<-rbind(teffectwl,twl)
tmledidmu<-rbind(teffectmu,tmu)

tmle<-cbind(tmledidt,tmledidwl,tmledidmu)

print.xtable(xtable(tmle),
             type="latex", file=paste0(path,"tmlespatial.tex"))

print.xtable(xtable(tmles),
             type="latex", file=paste0(path,"tmlesspatial.tex"))


##################################################################


