
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
packages <- c("readxl","rdd","psych","xtable","splines","ck37r","data.table","matrixStats","tmle","xgboost", "MatchIt","gtools","statar","foreign","multiwayvcov","lmtest","readstata13","xlsx", "data.table","doSNOW","parallel","compare","doParallel","devtools","foreach","spdep","reshape2","sm","plyr","utils","tcltk","geosphere", "matrixcalc", "dplyr","ExPosition", "randomForest","lfe", "hdm", "rdrobust", "stargazer", "ggplot2", "outliers","rpart","e1071")
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
#if(FALSE){
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
#}

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
psitel<-c(2,4,11,12,15,16,19,20,21)
psitel<-c(11,12,15,16,19,20,21)

#psite<-16
for(psite in psitel){
samplefull<-readRDS(paste(path,'fulldeletionbajgw',psite,'.rds', sep=""), refhook = NULL)

sample<-samplefull
sample$treatd0gw<-sample$treatdgw
#dist<-c('10k','8k','6k','5k','4k','3k','2k')#,'1k','500m')

dist<-c('10k','8k','6k','4k','2k')#,'1k','500m')
laglead<-c("")
#di<-4
#ll<-1

for(buf in 1:2){
#buf<-2
  sample$buffer<-ifelse(sample$date-odNPL$date[psite]-(buf*365)<0&sample$date-odNPL$date[psite]>0,1,0)

  sample<-sample[buffer<1,]
 


#matrices
for(i in c("lm","tmle")){
    for(j in c("t","wl","mu")){
        for(k in c("did","es")){
            assign(paste0('betas.',i,'.',j,'.',k),matrix(nrow = length(dist),ncol=length(laglead)))
            assign(paste0('ses.',i,'.',j,'.',k),matrix(nrow = length(dist),ncol=length(laglead)))
            assign(paste0('ps.',i,'.',j,'.',k),matrix(nrow = length(dist),ncol=length(laglead)))
        
            assign(paste0('cc.betas.',i,'.',j,'.',k),matrix(nrow = length(dist),ncol=length(laglead)))
            assign(paste0('cc.ses.',i,'.',j,'.',k),matrix(nrow = length(dist),ncol=length(laglead)))
            assign(paste0('cc.ps.',i,'.',j,'.',k),matrix(nrow = length(dist),ncol=length(laglead)))
            
              }
      }
  }

for(ll in 1:length(laglead)){
for(di in 1:length(dist)){
    dic<-dist[[di]]
    llc<-laglead[[ll]]
    sample<-samplefull[samplefull[[paste0('dist',dic)]]>0,]
    
    #Total Average Treatment Effect
    timefe<-select(sample, starts_with('timefe'))
    treatgroupm<-select(sample, starts_with('treatmentgroup'))
    year<-select(sample, starts_with('year'))
    bin<-select(sample, starts_with('bin'))
    sample$demlogprice<-demeanlist(sample$logprice,
                                   list(as.factor(sample$PropertyAddressCensusTractAndBlock)))
    sample$aftpropnpl<-sample[[paste0('aftpropnpl',psite)]]
    
    #sample$aftfinalnpl<-sample[[paste0('aftfinalnpl',fnplsite)]]
    #sample$timefinalnplfe<-sample[[paste0('timefinalnplfe',fnplsite)]]
    sample$treatControlsComplete<-sample[[paste0('treatControlsComplete',psite)]]
    sample$timefed<-sample[[paste0('timefed',psite)]]
    sample$timefedControlsComplete<-sample[[paste0('timefedControlsComplete',psite)]]
    sample$treatst<-sample[[paste0('treatd',llc,'gw',psite)]] #*sample[[paste0('dist',dist[[5]])]]
    sample$treatexCC<-ifelse(sample$treatControlsComplete==1 &sample$treatst==0,1,0)
    
    
    
    
    #TATE<-as.formula(logprice ~ treatdgw+ treatmentgroup+data.matrix(treatgroupm)+ data.matrix(timefe)+ data.matrix(year[,3:25])+ as.factor(HHID))
    
    sdf<-5
    lat<-sample$PropertyAddressLatitude
    long<-sample$PropertyAddressLongitude
    splat<-bs(lat, df = sdf)
    splong<-bs(long, df = sdf)
    spint<-model.matrix(~splat:splong)
    
    spTATE<-cbind(splat,splong,spint,lat,long)
    
    
    xTATE<-model.matrix(~ LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                          sqfeet+day+
                          data.matrix(year[,4:25])-1,sample)
    xcTATE<-model.matrix(~(data.matrix(bin)+
                             poly(LotSizeSquareFeet,4) + poly(YearBuilt,4) +   
                             poly(sqfeet,4))^3-1,sample)
    
    X<-model.matrix(~ treatst+treatexCC+treatmentgroup+
                             LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                             sqfeet+ 
                             data.matrix(year[,3:25]),sample)
    
    qr.X <- qr(X, tol=1e-3, LAPACK = FALSE)
    (rnkX <- qr.X$rank)  ## 4 (number of non-collinear columns)
    (keep <- qr.X$pivot[seq_len(rnkX)])
    ## 1 2 4 5 
    X <- X[,keep]
    
    results.lm.t.did<-felm(logprice ~X|as.factor(PropertyAddressCensusTractAndBlock),sample)
    summary(results.lm.t.did)
    nocc<-0
    notg<-0
    if(is.na(results.lm.t.did$coefficients["Xtreatst",])){
      X<-model.matrix(~ treatst+treatmentgroup+
                        LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                        sqfeet+ 
                        data.matrix(year[,3:25]),sample)
      
      qr.X <- qr(X, tol=1e-3, LAPACK = FALSE)
      (rnkX <- qr.X$rank)  ## 4 (number of non-collinear columns)
      (keep <- qr.X$pivot[seq_len(rnkX)])
      ## 1 2 4 5 
      X <- X[,keep]
      
      results.lm.t.did<-felm(logprice ~X|as.factor(PropertyAddressCensusTractAndBlock),sample)
      summary(results.lm.t.did)
      nocc<-1
      notg<-0
    }
    if(is.na(results.lm.t.did$coefficients["Xtreatst",])){
      X<-model.matrix(~ treatst+
                        LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                        sqfeet+ 
                        data.matrix(year[,3:25]),sample)
      
      qr.X <- qr(X, tol=1e-3, LAPACK = FALSE)
      (rnkX <- qr.X$rank)  ## 4 (number of non-collinear columns)
      (keep <- qr.X$pivot[seq_len(rnkX)])
      ## 1 2 4 5 
      X <- X[,keep]
      
      results.lm.t.did<-felm(logprice ~X|as.factor(PropertyAddressCensusTractAndBlock),sample)
      summary(results.lm.t.did)
      nocc<-1
      notg<-1
    }
    
    tsample<-sample[control==0,]
    year<-select(tsample, starts_with('year'))
    
    X<-model.matrix(~ treatst+treatexCC+
                               LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                               sqfeet+ 
                               data.matrix(year[,3:25]) ,tsample)
    qr.X <- qr(X, tol=1e-3, LAPACK = FALSE)
    (rnkX <- qr.X$rank)  ## 4 (number of non-collinear columns)
    (keep <- qr.X$pivot[seq_len(rnkX)])
    ## 1 2 4 5 
    X <- X[,keep]
    
    results.lm.t.es<-felm(logprice ~X|as.factor(PropertyAddressCensusTractAndBlock),tsample)
    summary(results.lm.t.es)
    
    if(is.na(results.lm.t.es$coefficients["Xtreatst",])){
      X<-model.matrix(~ treatst+
                        LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                        sqfeet+ 
                        data.matrix(year[,3:25]) ,tsample)
      qr.X <- qr(X, tol=1e-3, LAPACK = FALSE)
      (rnkX <- qr.X$rank)  ## 4 (number of non-collinear columns)
      (keep <- qr.X$pivot[seq_len(rnkX)])
      ## 1 2 4 5 
      X <- X[,keep]
      
      results.lm.t.es<-felm(logprice ~X|as.factor(PropertyAddressCensusTractAndBlock),tsample)
      summary(results.lm.t.es)
    }
    
    
    #assign(paste0('lm.TATE.',k,i),lm.TATE)
    summary(results.lm.t.did)
    summary(results.lm.t.es)
    
    
    sdf<-5
    lat<-sample$PropertyAddressLatitude
    long<-sample$PropertyAddressLongitude
    splat<-bs(lat, df = sdf)
    splong<-bs(long, df = sdf)
    spint<-model.matrix(~splat:splong)
    
    xcTATE<-cbind(splat,splong,spint,lat,long,poly(sample$day,5),bs(sample$day, df = 10))
    year<-select(sample, starts_with('year'))
   if(nocc==0&notg==0){
     feTATE<-model.matrix(~ treatmentgroup+#timefed +
                           treatexCC,sample)#+#timefedControlsComplete+
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
   }
    if(nocc==1&notg==0){
      feTATE<-model.matrix(~ treatmentgroup-1,sample)#+#timefed +
                             #timefedControlsComplete+
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
    }
    if(nocc==1&notg==1){
      feTATE<-model.matrix(~ #treatmentgroup+#timefed +
                             #treatexCC+#timefedControlsComplete+
                             #aftfinalnpl+timefinalnplfe+
                             #data.matrix(exdum)+
                             data.matrix(year[,25]),sample)
      feTATE<-as.matrix(feTATE[,SD(feTATE)>0])
      feTATE<-as.matrix(feTATE[,!duplicated(cor(feTATE))])
      qr.X <- qr(feTATE, tol=1e-3, LAPACK = FALSE)
      (rnkX <- qr.X$rank)  ## 4 (number of non-collinear columns)
      (keep <- qr.X$pivot[seq_len(rnkX)])
      ## 1 2 4 5 
      feTATE <- feTATE[,keep]
    }
    
    W<-cbind(xcTATE,xTATE,feTATE)
    A<-sample$treatst
    V<-feTATE
    Time<-sample$day
    
    results.tmle.t.did <- tmleMSM(Y = sample$logprice, A = A, W = W, V = V, T= Time,
                                  MSM = "A + V",family="gaussian", 
                                  Q.SL.library = SL.library2$as.list(),
                                  g.SL.library = PS.library2$as.list(),
                                  #Qform = Y ~ A+V+W,
                                  gform = A~1,
                                  #hAVform = A~ as.matrix(xTATE),
                                  ub = 20,
                                  V_SL =5,
                                  alpha = 0.90,
                                  inference = TRUE,
                                  verbose=TRUE)
    print(results.tmle.t.did)
    
    
    
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
    if(nocc==0){
    W<-cbind(xcTATE[sample$control==0,],xTATE[sample$control==0,],feTATEes)
    
    A<-sample[control==0,treatst]
    V<-feTATEes
    Time<-sample[control==0,day]
    
    results.tmle.t.es <- tmleMSM(Y = tsample$logprice, A = A, W = W, V = V, #T= Time,
                                 MSM = "A + V",family="gaussian", 
                                 Q.SL.library = SL.library2$as.list(),
                                 g.SL.library = PS.library2$as.list(),
                                 #Qform = Y ~ A+V+W,
                                 gform = A~1,
                                 #hAVform = A~ as.matrix(xTATE),
                                 ub = 20,
                                 V_SL =5,
                                 alpha = 0.90,
                                 inference = TRUE,
                                 verbose=TRUE)
    print(results.tmle.t.es)
    
    betas.tmle.t.es[di,ll]<-results.tmle.t.es$psi["A"]
    ses.tmle.t.es[di,ll]<-results.tmle.t.es$se["A"]
    ps.tmle.t.es[di,ll]<-results.tmle.t.es$pvalue["A"]
    
    cc.betas.tmle.t.es[di,ll]<-results.tmle.t.es$psi["V"]
    cc.ses.tmle.t.es[di,ll]<-results.tmle.t.es$se["V"]
    cc.ps.tmle.t.es[di,ll]<-results.tmle.t.es$pvalue["V"]
    
    
    }
    
    
    
    betas.tmle.t.did[di,ll]<-results.tmle.t.did$psi["A"]
    ses.tmle.t.did[di,ll]<-results.tmle.t.did$se["A"]
    ps.tmle.t.did[di,ll]<-results.tmle.t.did$pvalue["A"]
    
    cc.betas.tmle.t.did[di,ll]<-results.tmle.t.did$psi["VtreatexCC"]
    cc.ses.tmle.t.did[di,ll]<-results.tmle.t.did$se["VtreatexCC"]
    cc.ps.tmle.t.did[di,ll]<-results.tmle.t.did$pvalue["VtreatexCC"]
    
    
    betas.lm.t.did[di,ll]<-as.numeric(coef(summary(results.lm.t.did))[,"Estimate"][2])
    betas.lm.t.es[di,ll]<-as.numeric(coef(summary(results.lm.t.es))[,"Estimate"][2])
    ses.lm.t.did[di,ll]<-as.numeric(coef(summary(results.lm.t.did))[,"Std. Error"][2])
    ses.lm.t.es[di,ll]<-as.numeric(coef(summary(results.lm.t.es))[,"Std. Error"][2])
    ps.lm.t.did[di,ll]<-as.numeric(coef(summary(results.lm.t.did))[,"Pr(>|t|)"][2])
    ps.lm.t.es[di,ll]<-as.numeric(coef(summary(results.lm.t.es))[,"Pr(>|t|)"][2])
    
    cc.betas.lm.t.did[di,ll]<-as.numeric(coef(summary(results.lm.t.did))[,"Estimate"][3])
    cc.betas.lm.t.es[di,ll]<-as.numeric(coef(summary(results.lm.t.es))[,"Estimate"][3])
    cc.ses.lm.t.did[di,ll]<-as.numeric(coef(summary(results.lm.t.did))[,"Std. Error"][3])
    cc.ses.lm.t.es[di,ll]<-as.numeric(coef(summary(results.lm.t.es))[,"Std. Error"][3])
    cc.ps.lm.t.did[di,ll]<-as.numeric(coef(summary(results.lm.t.did))[,"Pr(>|t|)"][3])
    cc.ps.lm.t.es[di,ll]<-as.numeric(coef(summary(results.lm.t.es))[,"Pr(>|t|)"][3])
    
    print(paste0('distance = ',dic))
    print(paste0('laglead = ',llc))
    
  }
}

for(meth in c('lm','tmle')){
  for(inf in c('did','es')){
  #meth<-'tmle'
p<-get(paste0('ps.',meth,'.t.',inf))
mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ", ifelse(p < .1, "^\\bullet  ", " "))))

pb<-exp(get(paste0('betas.',meth,'.t.',inf)))-1
rpb<-round(pb,3)
se<-round(exp(get(paste0('ps.',meth,'.t.',inf)))-1,3)

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

#results.mat<-cbind(c('10k','','8k','','6k','','4k','','2k',''),results.mat)
#results.mat<-rbind(c('',laglead),results.mat)
xtable(results.mat)
print.xtable(xtable(results.mat),include.rownames=FALSE, 
             include.colnames=FALSE, sanitize.text.function = identity,
             type="latex", file=paste0(path,meth,inf,psite,"buffer",buf,".tex"))
}
}
}
}

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
