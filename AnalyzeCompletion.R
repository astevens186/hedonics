
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
packages <- c("readxl","matrixStats","tmle","xgboost", "MatchIt","gtools","statar","foreign","multiwayvcov","lmtest","readstata13","xlsx", "data.table","doSNOW","parallel","compare","doParallel","devtools","foreach","spdep","reshape2","sm","plyr","utils","tcltk","geosphere", "matrixcalc", "dplyr","ExPosition", "randomForest","lfe", "hdm", "rdrobust", "stargazer", "ggplot2", "outliers","rpart","e1071")
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

opNPL<- pNPL[order(pNPL$date),] 
ofNPL<- fNPL[order(fNPL$date),] 
odNPL<- dNPL[order(dNPL$date),] 
opdNPL<- dNPL[order(pdNPL$date),] 

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

# different configurations.
glmnet.tune <- list(alpha = c(0,.1, .25,.5,.75,.9,1))

# Set detailed names = T so we can see the configuration for each function.
# Also shorten the name prefix.
glmnet <- create.Learner("SL.glmnet", tune = glmnet.tune, detailed_names = T, name_prefix = "glmnet")

# configurations
length(glmnet$names)
glmnet$names

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

SL.library$add("SL.mean")
SL.library$add("SL.randomForest")
SL.library$add("SL.xgboost")

for(i in 1:length(glmnet$names)){
  SL.library$add(glmnet$names[i])
}    

for(i in 1:length(xgboost$names)){
  SL.library$add(c(xgboost$names[i]))
}

for(i in 1:length(randomForest$names)){
  SL.library$add(c(randomForest$names[i]))
}

SL.library$as.list()

###############################################################################################
for(k in c("full")){
  for(i in 1:dim(odNPL)[1]){
    #k<-"full"
    #i<-73
    if (file.exists(paste(path,k,'ControlsComplete',i,'.rds', sep=""))){
      sample.1<-readRDS(paste(path,k,'ControlsComplete',i,'.rds', sep=""), refhook = NULL)
      
      sample<-sample.1
      #sample<-subset(sample,select=-c(63:321,325:585,610:1191))
      
      library(dplyr)
      
      #Total Average Treatment Effect
      timefe<-select(sample, starts_with('timefe'))
      treatgroupm<-select(sample, starts_with('treatmentgroup'))
      year<-select(sample, starts_with('year'))
      bin<-select(sample, starts_with('bin'))
      
      #TATE<-as.formula(logprice ~ treatdgw+ treatmentgroup+data.matrix(treatgroupm)+ data.matrix(timefe)+ data.matrix(year[,3:25])+ as.factor(HHID))
      feTATE<-model.matrix(~ treatmentgroup+data.matrix(treatgroupm)+ 
                             data.matrix(year[,4:25]),sample)
      xTATE<-model.matrix(~as.numeric(preprice)+as.numeric(prelogprice)+
                            as.numeric(predate)+as.numeric(prediffdate)+
                            LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                            sqfeet,sample)
      xcTATE<-model.matrix(~(poly(as.numeric(preprice),2)+poly(as.numeric(prelogprice),2)+
                               poly(as.numeric(predate),2)+poly(as.numeric(prediffdate),2)+ 
                               data.matrix(bin))^3+
                             poly(LotSizeSquareFeet,4) + poly(YearBuilt,4) +   
                             poly(sqfeet,4),sample)
      
      lm.TATE<-felm(logprice ~treatControlsComplete+ as.numeric(preprice)+as.numeric(prelogprice)+
                      as.numeric(predate)+as.numeric(prediffdate)+
                      LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                      sqfeet+ data.matrix(treatgroupm)
                    + data.matrix(year[,3:25])|as.factor(PropertyAddressCensusTractAndBlock),sample)
      assign(paste0('lm.TATE.',k,i),lm.TATE)
      summary(lm.TATE)
      
      #Well Average Treatment Effect
      sample$treatgwWL<- sample$treatControlsComplete * sample$WaterStndCode.fWL
      sample$treatgroupWL<-sample$treatmentgroup * sample$WaterStndCode.fWL
      sample$controlWL<-sample$control*sample$WaterStndCode.fWL
      
      sample$sample.WLATE<-sample$control+sample$treatgroupWL
      
      sample.WLATE<-subset(sample, sample.WLATE==1)
      timefe<-select(sample.WLATE, starts_with('timefe'))
      treatgroupm<-select(sample.WLATE, starts_with('treatmentgroup'))
      year<-select(sample.WLATE, starts_with('year'))
      bin<-select(sample.WLATE, starts_with('bin'))
      
      #WLATE<-as.formula(logprice ~ treatgwWL+ treatgroupWL+data.matrix(treatgroupm)+ data.matrix(timefe)+ data.matrix(year[,3:25])+ as.factor(HHID))
      feWLATE<-model.matrix(~ treatgroupWL+data.matrix(treatgroupm)+ 
                              data.matrix(year[,4:25]),sample.WLATE)
      xWLATE<-model.matrix(~as.numeric(preprice)+as.numeric(prelogprice)+
                             as.numeric(predate)+as.numeric(prediffdate)+
                             LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                             sqfeet,sample.WLATE)
      xcWLATE<-model.matrix(~(poly(as.numeric(preprice),2)+poly(as.numeric(prelogprice),2)+
                                poly(as.numeric(predate),2)+poly(as.numeric(prediffdate),2)+ 
                                data.matrix(bin))^3+
                              poly(LotSizeSquareFeet,4) + poly(YearBuilt,4) +   
                              poly(sqfeet,4),sample.WLATE)
      
      lm.WLATE<-felm(logprice ~ treatgwWL+as.numeric(preprice)+as.numeric(prelogprice)+
                       as.numeric(predate)+as.numeric(prediffdate)+
                       LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                       sqfeet+treatgroupWL+ data.matrix(treatgroupm)
                     +data.matrix(year[,3:25])|as.factor(PropertyAddressCensusTractAndBlock),sample.WLATE)
      assign(paste0('lm.WLATE.',k,i),lm.WLATE)
      summary(lm.WLATE)
      
      #Municipal Utility Average Treatment Effect
      sample$treatgwMU<- sample$treatControlsComplete * sample$WaterStndCode.fMU
      sample$treatgroupMU<-sample$treatmentgroup * sample$WaterStndCode.fMU
      sample$controlMU<-sample$control*sample$WaterStndCode.fMU
      
      sample$sample.MUATE<-sample$control+sample$treatgroupMU
      
      sample.MUATE<-subset(sample, sample.MUATE==1)
      timefe<-select(sample.MUATE, starts_with('timefe'))
      treatgroupm<-select(sample.MUATE, starts_with('treatmentgroup'))
      year<-select(sample.MUATE, starts_with('year'))
      bin<-select(sample.MUATE, starts_with('bin'))
      
      #MUATE<-as.formula(logprice ~ treatgwMU+ treatgroupMU+data.matrix(treatgroupm)+ data.matrix(timefe)+ data.matrix(year[,3:25])+ as.factor(HHID))
      feMUATE<-model.matrix(~ treatgroupMU+data.matrix(treatgroupm)+ 
                              data.matrix(year[,4:25]),sample.MUATE)
      xMUATE<-model.matrix(~as.numeric(preprice)+as.numeric(prelogprice)+
                             as.numeric(predate)+as.numeric(prediffdate)+
                             LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                             sqfeet,sample.MUATE)
      xcMUATE<-model.matrix(~(poly(as.numeric(preprice),2)+poly(as.numeric(prelogprice),2)+
                                poly(as.numeric(predate),2)+poly(as.numeric(prediffdate),2)+ 
                                data.matrix(bin))^3+
                              poly(LotSizeSquareFeet,4) + poly(YearBuilt,4) +   
                              poly(sqfeet,4),sample.MUATE)
      
      lm.MUATE<-felm(logprice ~treatgwMU+as.numeric(preprice)+as.numeric(prelogprice)+
                       as.numeric(predate)+as.numeric(prediffdate)+
                       LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                       sqfeet+  treatgroupMU+ data.matrix(treatgroupm)
                     + data.matrix(year[,3:25])|as.factor(PropertyAddressCensusTractAndBlock),sample.MUATE)
      assign(paste0('lm.MUATE.',k,i),lm.MUATE)
      summary(lm.MUATE)
      #######Estimation
      
      xxTATE<-cbind(xcTATE,xTATE,feTATE)
      xxWLATE<-cbind(xcWLATE,xWLATE,feWLATE)
      xxMUATE<-cbind(xcMUATE,xMUATE,feMUATE)
      
      sample$demlogprice<-demeanlist(sample$logprice,
                                     list(as.factor(sample$HHID)))
      sample.WLATE$demlogprice<-demeanlist(sample.WLATE$logprice,
                                           list(as.factor(sample.WLATE$HHID)))
      sample.MUATE$demlogprice<-demeanlist(sample.MUATE$logprice,
                                           list(as.factor(sample.MUATE$HHID)))
      
      FxcTATE<-logical(dim(xcTATE)[2])
      TxTATE<-!logical(dim(xTATE)[2])
      TfeTATE<-!logical(dim(feTATE)[2])
      indTATE<-rbind(data.matrix(FxcTATE),data.matrix(TxTATE),data.matrix(TfeTATE))
      
      fit.ds.os1 <- rlassoEffect(x=xxTATE, y=sample$demlogprice, d=sample$treatControlsComplete,data = sample,
                                 method = "double selection", I3=indTATE)
      assign(paste0('fit.ds.os1.',k,i),fit.ds.os1)
      summary(fit.ds.os1)
      
      FxcWLATE<-logical(dim(xcWLATE)[2])
      TxWLATE<-!logical(dim(xWLATE)[2])
      TfeWLATE<-!logical(dim(feWLATE)[2])
      indWLATE<-rbind(data.matrix(FxcWLATE),data.matrix(TxWLATE),data.matrix(TfeWLATE))
      
      fit.ds.os2 <- rlassoEffect(x=xxWLATE, y=sample.WLATE$demlogprice, d=sample.WLATE$treatgwWL,data = sample.WLATE, 
                                 method = "double selection", I3=indWLATE)
      assign(paste0('fit.ds.os2.',k,i),fit.ds.os2)
      summary(fit.ds.os2)
      
      FxcMUATE<-logical(dim(xcMUATE)[2])
      TxMUATE<-!logical(dim(xMUATE)[2])
      TfeMUATE<-!logical(dim(feMUATE)[2])
      indMUATE<-rbind(data.matrix(FxcMUATE),data.matrix(TxMUATE),data.matrix(TfeMUATE))
      
      fit.ds.os3 <- rlassoEffect(x=xxMUATE, y=sample.MUATE$demlogprice, d=sample.MUATE$treatgwMU, data = sample.MUATE,
                                 method = "double selection", I3=indMUATE)
      assign(paste0('fit.ds.os3.',k,i),fit.ds.os3)
      summary(fit.ds.os3)
      
      ## TMLE
      
      xfTATE<-cbind(xTATE,feTATE)
      xfWLATE<-cbind(xWLATE,feWLATE)
      xfMUATE<-cbind(xMUATE,feMUATE)
      
      result1 <- tmle(Y=sample$demlogprice,A=sample$treatControlsComplete,W=xfTATE,
                      Q.SL.library = SL.library$as.list(),
                      g.SL.library = SL.library$as.list())
      assign(paste0('fit.tmle.os1.',k,i),result1)
      summary(result1)
      
      result2 <- tmle(Y=sample.WLATE$demlogprice,A=sample.WLATE$treatgwWL,W=xfWLATE,
                      Q.SL.library = SL.library$as.list(),
                      g.SL.library = SL.library$as.list())
      assign(paste0('fit.tmle.os2.',k,i),result2)
      summary(result2)
      
      result3 <- tmle(Y=sample.MUATE$demlogprice,A=sample.MUATE$treatgwMU,W=xfMUATE,
                      Q.SL.library = SL.library$as.list(),
                      g.SL.library = SL.library$as.list())
      assign(paste0('fit.tmle.os3.',k,i),result3)
      summary(result3)
    }
    print(paste0(k,"and",i))
  }
}


coef(summary(lm.TATE.full73))['treatControlsComplete',]
summary(fit.ds.os1.full73)
coef(summary(lm.WLATE.full73))['treatgwWL',]
summary(fit.ds.os2.full73)
coef(summary(lm.MUATE.full73))['treatgwMU',]
summary(fit.ds.os3.full73)

coef(summary(lm.TATE.full84))['treatControlsComplete',]
summary(fit.ds.os1.full84)
coef(summary(lm.WLATE.full84))['treatgwWL',]
summary(fit.ds.os2.full84)
coef(summary(lm.MUATE.full84))['treatgwMU',]
summary(fit.ds.os3.full84)

coef(summary(lm.TATE.full88))['treatControlsComplete',]
summary(fit.ds.os1.full88)
coef(summary(lm.WLATE.full88))['treatgwWL',]
summary(fit.ds.os2.full88)
coef(summary(lm.MUATE.full88))['treatgwMU',]
summary(fit.ds.os3.full88)

coef(summary(lm.TATE.full91))['treatControlsComplete',]
summary(fit.ds.os1.full91)
coef(summary(lm.WLATE.full91))['treatgwWL',]
summary(fit.ds.os2.full91)
coef(summary(lm.MUATE.full91))['treatgwMU',]
summary(fit.ds.os3.full91)

#sample.1<-readRDS(paste(path,k,'deletiongw',i,'.rds', sep=""), refhook = NULL)
for(k in c("full","mdm","psm")){
  sample.1<-readRDS(paste(path,k,'ControlsComplete',73,'.rds', sep=""), refhook = NULL)
  for(i in 89:dim(odNPL)[1]){
    if (file.exists(paste(path,k,'ControlsComplete',i,'.rds', sep=""))){
      sample.2<-readRDS(paste(path,k,'ControlsComplete',i,'.rds', sep=""), refhook = NULL)
      sample.1<-smartbind(sample.1,sample.2)
    }
  }
}


sample<-sample.1
#sample<-subset(sample,select=-c(63:321,325:585,610:1191))

library(dplyr)

#Total Average Treatment Effect
timefe<-select(sample, starts_with('timefe'))
treatgroupm<-select(sample, starts_with('treatmentgroup'))
year<-select(sample, starts_with('year'))
bin<-select(sample, starts_with('bin'))

#TATE<-as.formula(logprice ~ treatdgw+ treatmentgroup+data.matrix(treatgroupm)+ data.matrix(timefe)+ data.matrix(year[,3:25])+ as.factor(HHID))
feTATE<-model.matrix(~ treatmentgroup+data.matrix(treatgroupm)+ 
                       data.matrix(year[,4:25]),sample)
xTATE<-model.matrix(~as.numeric(preprice)+as.numeric(prelogprice)+
                      as.numeric(predate)+as.numeric(prediffdate)+
                      LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                      sqfeet,sample)
xcTATE<-model.matrix(~(poly(as.numeric(preprice),2)+poly(as.numeric(prelogprice),2)+
                         poly(as.numeric(predate),2)+poly(as.numeric(prediffdate),2)+ 
                         data.matrix(bin))^3+
                       poly(LotSizeSquareFeet,4) + poly(YearBuilt,4) +   
                       poly(sqfeet,4),sample)

lm.TATE<-felm(logprice ~treatControlsComplete+ as.numeric(preprice)+as.numeric(prelogprice)+
                as.numeric(predate)+as.numeric(prediffdate)+
                LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                sqfeet+ data.matrix(treatgroupm)   + 
              + data.matrix(year[,3:25])|as.factor(PropertyAddressCensusTractAndBlock),sample)
assign(paste0('lm.TATE.',k,i),lm.TATE)
summary(lm.TATE)

#Well Average Treatment Effect
sample$treatgwWL<- sample$treatControlsComplete * sample$WaterStndCode.fWL
sample$treatgroupWL<-sample$treatmentgroup * sample$WaterStndCode.fWL
sample$controlWL<-sample$control*sample$WaterStndCode.fWL

sample$sample.WLATE<-sample$control+sample$treatgroupWL

sample.WLATE<-subset(sample, sample.WLATE==1)
timefe<-select(sample.WLATE, starts_with('timefe'))
treatgroupm<-select(sample.WLATE, starts_with('treatmentgroup'))
year<-select(sample.WLATE, starts_with('year'))
bin<-select(sample.WLATE, starts_with('bin'))

#WLATE<-as.formula(logprice ~ treatgwWL+ treatgroupWL+data.matrix(treatgroupm)+ data.matrix(timefe)+ data.matrix(year[,3:25])+ as.factor(HHID))
feWLATE<-model.matrix(~ treatgroupWL+data.matrix(treatgroupm)+ 
                        data.matrix(year[,4:25]),sample.WLATE)
xWLATE<-model.matrix(~as.numeric(preprice)+as.numeric(prelogprice)+
                       as.numeric(predate)+as.numeric(prediffdate)+
                       LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                       sqfeet,sample.WLATE)
xcWLATE<-model.matrix(~(poly(as.numeric(preprice),2)+poly(as.numeric(prelogprice),2)+
                          poly(as.numeric(predate),2)+poly(as.numeric(prediffdate),2)+ 
                          data.matrix(bin))^3+
                        poly(LotSizeSquareFeet,4) + poly(YearBuilt,4) +   
                        poly(sqfeet,4),sample.WLATE)


lm.WLATE<-felm(logprice ~ treatgwWL+as.numeric(preprice)+as.numeric(prelogprice)+
                 as.numeric(predate)+as.numeric(prediffdate)+
                 LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                 sqfeet+treatgroupWL+ data.matrix(treatgroupm)   + 
               + data.matrix(year[,3:25])|as.factor(PropertyAddressCensusTractAndBlock),sample.WLATE)
assign(paste0('lm.WLATE.',k,i),lm.WLATE)
summary(lm.WLATE)

#Municipal Utility Average Treatment Effect
sample$treatgwMU<- sample$treatControlsComplete * sample$WaterStndCode.fMU
sample$treatgroupMU<-sample$treatmentgroup * sample$WaterStndCode.fMU
sample$controlMU<-sample$control*sample$WaterStndCode.fMU

sample$sample.MUATE<-sample$control+sample$treatgroupMU

sample.MUATE<-subset(sample, sample.MUATE==1)
timefe<-select(sample.MUATE, starts_with('timefe'))
treatgroupm<-select(sample.MUATE, starts_with('treatmentgroup'))
year<-select(sample.MUATE, starts_with('year'))
bin<-select(sample.MUATE, starts_with('bin'))

#MUATE<-as.formula(logprice ~ treatgwMU+ treatgroupMU+data.matrix(treatgroupm)+ data.matrix(timefe)+ data.matrix(year[,3:25])+ as.factor(HHID))
feMUATE<-model.matrix(~ treatgroupMU+data.matrix(treatgroupm)+ 
                        data.matrix(year[,4:25]),sample.MUATE)
xMUATE<-model.matrix(~as.numeric(preprice)+as.numeric(prelogprice)+
                       as.numeric(predate)+as.numeric(prediffdate)+
                       LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                       sqfeet,sample.MUATE)
xcMUATE<-model.matrix(~(poly(as.numeric(preprice),2)+poly(as.numeric(prelogprice),2)+
                          poly(as.numeric(predate),2)+poly(as.numeric(prediffdate),2)+ 
                        data.matrix(bin))^3+
                        poly(LotSizeSquareFeet,4) + poly(YearBuilt,4) +   
                        poly(sqfeet,4),sample.MUATE)

lm.MUATE<-felm(logprice ~treatgwMU+as.numeric(preprice)+as.numeric(prelogprice)+
                 as.numeric(predate)+as.numeric(prediffdate)+
                 LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                 sqfeet+  treatgroupMU+ data.matrix(treatgroupm)   + 
               + data.matrix(year[,3:25])|as.factor(PropertyAddressCensusTractAndBlock),sample.MUATE)
assign(paste0('lm.MUATE.',k,i),lm.MUATE)
summary(lm.MUATE)
#######Estimation

xxTATE<-cbind(xcTATE,xTATE,feTATE)
xxWLATE<-cbind(xcWLATE,xWLATE,feWLATE)
xxMUATE<-cbind(xcMUATE,xMUATE,feMUATE)

sample$demlogprice<-demeanlist(sample$logprice,
                               list(as.factor(sample$PropertyAddressCensusTractAndBlock)))
sample.WLATE$demlogprice<-demeanlist(sample.WLATE$logprice,
                                     list(as.factor(sample.WLATE$PropertyAddressCensusTractAndBlock)))
sample.MUATE$demlogprice<-demeanlist(sample.MUATE$logprice,
                                     list(as.factor(sample.MUATE$PropertyAddressCensusTractAndBlock)))
FxcTATE<-logical(dim(xcTATE)[2])
TxTATE<-!logical(dim(xTATE)[2])
TfeTATE<-!logical(dim(feTATE)[2])
indTATE<-rbind(data.matrix(FxcTATE),data.matrix(TxTATE),data.matrix(TfeTATE))

fit.ds.os1 <- rlassoEffect(x=xxTATE, y=sample$demlogprice, d=sample$treatControlsComplete,data = sample,
                           method = "double selection", I3=indTATE)
assign(paste0('fit.ds.os1.',k,i),fit.ds.os1)
summary(fit.ds.os1)

FxcWLATE<-logical(dim(xcWLATE)[2])
TxWLATE<-!logical(dim(xWLATE)[2])
TfeWLATE<-!logical(dim(feWLATE)[2])
indWLATE<-rbind(data.matrix(FxcWLATE),data.matrix(TxWLATE),data.matrix(TfeWLATE))

fit.ds.os2 <- rlassoEffect(x=xxWLATE, y=sample.WLATE$demlogprice, d=sample.WLATE$treatgwWL,data = sample.WLATE, 
                           method = "double selection", I3=indWLATE)
assign(paste0('fit.ds.os2.',k,i),fit.ds.os2)
summary(fit.ds.os2)

FxcMUATE<-logical(dim(xcMUATE)[2])
TxMUATE<-!logical(dim(xMUATE)[2])
TfeMUATE<-!logical(dim(feMUATE)[2])
indMUATE<-rbind(data.matrix(FxcMUATE),data.matrix(TxMUATE),data.matrix(TfeMUATE))


fit.ds.os3 <- rlassoEffect(x=xxMUATE, y=sample.MUATE$demlogprice, d=sample.MUATE$treatgwMU, data = sample.MUATE,
                           method = "double selection", I3=indMUATE)
assign(paste0('fit.ds.os3.',k,i),fit.ds.os3)
summary(fit.ds.os3)



