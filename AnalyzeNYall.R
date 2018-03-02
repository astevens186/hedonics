
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
packages <- c("readxl","Hmisc","sphet","mgcv","McSpatial","pastecs","rdd","Matrix","psych","xtable","splines","ck37r","data.table","matrixStats","tmle","xgboost", "MatchIt","gtools","statar","foreign","multiwayvcov","lmtest","readstata13","xlsx", "data.table","doSNOW","parallel","compare","doParallel","devtools","foreach","spdep","reshape2","sm","plyr","utils","tcltk","geosphere", "matrixcalc", "dplyr","ExPosition", "randomForest","lfe", "hdm", "rdrobust", "stargazer", "ggplot2", "outliers","rpart","e1071")
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
odNPL<-odNPL[odNPL$rstate_code == "NY",]
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

xgboost.tune2 <- list(ntrees = c(50),
                      max_depth = c(5,15),
                      shrinkage = c( 0.1),
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

glmnet.tune2 <- list(alpha = c(0,1))

# Set detailed names = T so we can see the configuration for each function.
# Also shorten the name prefix.
glmnet2 <- create.Learner("SL.glmnet", tune = glmnet.tune2, detailed_names = T, name_prefix = "glmnet")

# configurations
length(glmnet2$names)
glmnet2$names

# different configurations.
randomForest.tune <- list(ntree = c(500))

# Set detailed names = T so we can see the configuration for each function.
# Also shorten the name prefix.
randomForest <- create.Learner("SL.randomForest", tune = randomForest.tune, 
                               detailed_names = T, name_prefix = "randomForest")

# configurations
length(randomForest$names)
randomForest$names


# different configurations.
gam.tune <- list(spatialsp= c("ts","gp","tp"))#, cts.num = 10)

# Set detailed names = T so we can see the configuration for each function.
# Also shorten the name prefix.
gam <- create.Learner("SL.gam", tune = gam.tune, detailed_names = T, name_prefix = "gam")

# configurations
length(gam$names)
gam$names

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

if(FALSE){
  
  
  results.gam<-mgcv::gam(logprice~treatind+X+s(day,lat,long,bs="tp",m=3,k=50),data=sample)
  mgcv::summary.gam(results.gam)   
  Xs<-cbind(X,lat,long)
}
  SL.gam <- function(Y, X, newX, family, obsWeights,ms=3,ks=50,ksi=10,kt=10,spatialsp,temporalsp="gp",
                     deg.gam =2 , cts.num=4 ,slat=lat,slong=long, ...) {
    .SL.require('mgcv')
    s=mgcv:::s
    cts.x <- apply(X, 2, function(x) (length(unique(x)) > cts.num))
    cts.x["lat"] <- FALSE
    cts.x["long"]<- FALSE
    cts.x["day"]<- FALSE
    if (sum(!cts.x) > 0) { 
      gam.model <- as.formula(paste("Y~", paste(colnames(X[, cts.x, drop = FALSE]), 
                                                collapse = "+"), "+s(lat,long,bs='",spatialsp,"',m=",ms,",k=",ks,")",sep=""))
                                               # "+s(day,bs='",temporalsp,"',m=",ms,",k=",kt,")+",
                                    #paste(paste("s(", colnames(X[, cts.x, drop = FALSE]),",lat,long,bs='",spatialsp,"',m=",ms,",k=",ksi,")", sep=""), 
                                               # collapse = "+"),
                                                #"+s(day,lat,long,bs='",spatialsp,"',m=",ms,",k=",ks,")", sep=""))
    } else {
      gam.model <- as.formula(paste("Y~", paste(paste("s(", colnames(X[, cts.x, drop = FALSE]), ",", deg.gam, ")", sep=""), collapse = "+")))
    }
    # fix for when all variables are binomial
    if (sum(!cts.x) == length(cts.x)) {
      gam.model <- as.formula(paste("Y~", paste(colnames(X), collapse = "+"), sep = ""))
    }
    fit.gam <- mgcv::gam(gam.model, data = X, family = family)
    pred <-mgcv::predict.gam(fit.gam, newdata = newX, type = "response")
    fit <- list(object = fit.gam)
    out <- list(pred = pred, fit = fit)
    class(out$fit) <- c("SL.gam")
    return(out)
  }
  
  predict.SL.gam <- function(object, newdata, ...){
    .SL.require('gam')
    pred <- gam::predict.gam(object = object$object, newdata = newdata, type = "response")
    return(pred)
  }
  
  .SL.require <- function(package, message = paste('loading required package (', package, ') failed', sep = '')) {
    if(!require(package, character.only = TRUE)) {
      stop(message, call. = FALSE)
    }
    invisible(TRUE)
  }


SL.library<-expandingList()
PS.library<-expandingList()
SL.library$add("SL.mean")
PS.library$add("SL.mean")
SL.library2<-expandingList()
PS.library2<-expandingList()
SL.library2$add("SL.mean")
PS.library2$add("SL.mean")
#SL.library2$add(randomForest$names)
#SL.library2$add("SL.nnls")
#SL.library2$add("SL.gam")
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
#for(i in 1:length(gam$names)){
#  SL.library$add(gam$names[i])
#}
for(i in 1:length(gam$names)){
 SL.library2$add(gam$names[i])
}

for(i in 1:length(xgboost$names)){
  SL.library$add(c(xgboost$names[i]))
}

#for(i in 1:length(xgboost2$names)){
#  SL.library2$add(c(xgboost2$names[i]))
#}

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
dNPL$row<-seq(1:dim(dNPL)[1])
psitel<-dNPL[dNPL$rsitinc_desc=="LANDFILL","row"]

psitel<-c(2,4,11,12,15,16,19,20,21)
psitel<-c(2,12,15,16)
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
  sample<-sample[,c("TransId","preprice","prelogprice","predate","prediffdate","presstatusd","presstatuscc","difflogprice")]
  data.dt<-data.table(data)
  sample.dt<-data.table(sample)
  sample.new<-merge(data.dt,sample.dt,all.x=TRUE,by="TransId")
  
  saveRDS(sample.new, file = paste(path,'fullbaj',psite,'.rds', sep=""), ascii = FALSE, version = NULL,
          compress = TRUE, refhook = NULL)
}
for(treat in  1:length(treatl)){
  
  treatc<-treatl[[treat]]
  #for(ll in 1:length(laglead)){
  for(di in 1:length(dist)){
    dic<-dist[[di]]
    ll<-1
    llc<-laglead[[ll]]
    sample<-samplefull[samplefull[[paste0('dist',dic)]]>0,]
    #sample$treatst<-sample$treatst-sample$presstatusd
    #sample<-sample[treatst>0,]
    sample<-sample[presstatusd==0,]
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
        t<-d.sample.data$date
        rep.row<-function(x,n){
          matrix(rep(x,each=n),nrow=n)
        }
        t<-d.sample.data$date
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
      
      
      saveRDS(Wst, file = path, ascii = FALSE, version = NULL,
              compress = TRUE, refhook = NULL)
      #stop cluster
      
    }
    start.time <- Sys.time()
    
    W.trend.lag.variables(upper.spatial.range,lower.spatial.range,spatial.power.range,
                          temporal.cut.range,temporal.power.range,paste0(path,"Wmat",treatc,dic,".rds"))
    
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    time.taken
  }
}

psitel<-c(2,15,16)
psitel<-c(2,11,12,15,16,19,20,21)


dNPL$row<-seq(1:dim(dNPL)[1])
psitel<-dNPL[dNPL$rsitinc_desc=="LANDFILL","row"]

psite<-psitel[1]
samplefull<-readRDS(paste(path,'fullbaj',psite,'.rds', sep=""), refhook = NULL)
samplefull<-samplefull[predate>0,]
samplefull$treatst<-samplefull[[paste0('treatdgw',psite)]] #*sample[[paste0('dist',dist[[5]])]]
buf<-1
#samplefull$buffer<-ifelse(samplefull$date-odNPL$date[psite]+(buf*365)>0&samplefull$date-odNPL$date[psite]<0,1,0)
samplefull$buffer<-ifelse(abs(samplefull$date-odNPL$date[psite])-(buf*365)<0,1,0)
samplefull<-samplefull[buffer<1,]
samplefull$timetotreat<-samplefull$date-odNPL$date[psite]

for(psite in psitel[2:length(psitel)]){
  sample2<-readRDS(paste(path,'fullbaj',psite,'.rds', sep=""), refhook = NULL)
  sample2<-sample2[predate>0,]
  sample2$treatst<-sample2[[paste0('treatdgw',psite)]] #*sample[[paste0('dist',dist[[5]])]]
  buf<-1
  #sample2$buffer<-ifelse(sample2$date-odNPL$date[psite]+(buf*365)>0&sample2$date-odNPL$date[psite]<0,1,0)
  sample2$buffer<-ifelse(abs(sample2$date-odNPL$date[psite])-(buf*365)<0,1,0)
  sample2<-sample2[buffer<1,]
  sample2$timetotreat<-sample2$date-odNPL$date[psite]
  
  samplefull<-rbind(samplefull,sample2)
}

#psite<-2
samplefull$timetotreat<-samplefull$timetotreat/365
samplefull<-samplefull[abs(timetotreat)-10<0,]
samplefull$treatd0gw<-samplefull$treatdgw

quant<-10
qcut<-cut2(samplefull$preprice, g=quant, onlycuts = TRUE)

price.ag <- aggregate(price ~ treatmentgroup, data=samplefull, FUN=mean, na.rm=TRUE)
sqfeet.ag <- aggregate(sqfeet ~ treatmentgroup, data=samplefull, FUN=mean, na.rm=TRUE)
TotalRooms.ag <- aggregate(TotalRooms ~ treatmentgroup, data=samplefull, FUN=mean, na.rm=TRUE)
YearBuilt.ag <- aggregate(YearBuilt ~ treatmentgroup, data=samplefull, FUN=mean, na.rm=TRUE)
FullBath.ag <- aggregate(FullBath ~ treatmentgroup, data=samplefull, FUN=mean, na.rm=TRUE)

sumvar<-c("price","sqfeet","TotalRooms","YearBuilt","FullBath")
ts<-samplefull[treatmentgroup==1,c("price","sqfeet","YearBuilt","RecordingDate")]
ts$year <- factor(format(as.Date(ts$RecordingDate),'%Y'))
cs<-samplefull[treatmentgroup==0,c("price","sqfeet","YearBuilt","RecordingDate")]
cs$year <- factor(format(as.Date(cs$RecordingDate),'%Y'))
B <- matrix( c(min(ts$price), min(cs$price), median(ts$price), median(cs$price),
              mean(ts$price), mean(cs$price),max(ts$price),max(cs$price),
              min(ts$sqfeet), min(cs$sqfeet), median(ts$sqfeet), median(cs$sqfeet),
              mean(ts$sqfeet), mean(cs$sqfeet),max(ts$sqfeet),max(cs$sqfeet),
              min(as.numeric(as.character(ts$year))), min(as.numeric(as.character(cs$year))), 
              median(as.numeric(as.character(ts$year))), median(as.numeric(as.character(cs$year))),
              mean(as.numeric(as.character(ts$year))), mean(as.numeric(as.character(cs$year))),
              max(as.numeric(as.character(ts$year))),max(as.numeric(as.character(cs$year))),
              length(ts$year),length(cs$year)), nrow=2, ncol=13)
sumtab<-t(floor(B))
sumtab[9:12,]<-as.Date(sumtab[9:12,])
colnames(sumtab)<-c("Treatment Group","Control Group")
rownames(sumtab)<-c("Min. Price","Med. Price", "Mean Price", "Max. Price",
                    "Min. Sq. Feet","Med. Sq. Feet", "Mean Sq. Feet", "Max. Sq. Feet",
                    "Min. Year","Med. Year", "Mean Year", "Max. Year","Obs")

xtable(sumtab)
print.xtable(xtable(floor(sumtab),digits=c(0,0,0)),include.rownames=TRUE, 
             include.colnames=TRUE, sanitize.text.function = identity,
             type="latex", file=paste0(path,'latex/sumtab.tex'))

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



#matrices

for(i in c("lm","gam","sp")){
  for(j in c("t","wl","mu")){
    for(k in c("did")){
      for(treat in treatl){
      assign(paste0('betas.',i,'.',j,'.',k,'.',treat),matrix(ncol = length(dist),nrow=quant))
      assign(paste0('ses.',i,'.',j,'.',k,'.',treat),matrix(ncol = length(dist),nrow=quant))
      assign(paste0('ps.',i,'.',j,'.',k,'.',treat),matrix(ncol = length(dist),nrow=quant))
      
      assign(paste0('betas.',i,'.',k),matrix(ncol = length(dist),nrow=length(treatl)))
      assign(paste0('ses.',i,'.',k),matrix(ncol = length(dist),nrow=length(treatl)))
      assign(paste0('ps.',i,'.',k),matrix(ncol = length(dist),nrow=length(treatl)))
       
      
}
    }
  }
}
di<-3
treat<-1

for(treat in  1:length(treatl)){
  
  treatc<-treatl[[treat]]
  #for(ll in 1:length(laglead)){
  for(di in 1:length(dist)){
    dic<-dist[[di]]
    ll<-1
    llc<-laglead[[ll]]
    sample<-samplefull[samplefull[[paste0('dist',dic)]]>0,]
    #sample$treatst<-sample$treatst-sample$presstatusd
    #sample<-sample[treatst>0,]
    sample<-sample[presstatusd==0,]
    # sample$logprice<-sample$logprice.x
    if(treatc=='TATE'){
      #Total Average Treatment Effect
      timefe<-dplyr::select(sample, starts_with('timefe'))
      treatgroupm<-dplyr::select(sample, starts_with('treatmentgroup'))
      year<-dplyr::select(sample, starts_with('year'))
      bin<-dplyr::select(sample, starts_with('bin'))
      sample$demlogprice<-demeanlist(sample$logprice,
                                     list(as.factor(sample$PropertyAddressCensusTractAndBlock)))
      sample$aftpropnpl<-sample[[paste0('aftpropnpl',psite)]]
      
      #sample$aftfinalnpl<-sample[[paste0('aftfinalnpl',fnplsite)]]
      #sample$timefinalnplfe<-sample[[paste0('timefinalnplfe',fnplsite)]]
      sample$treatControlsComplete<-sample[[paste0('treatControlsComplete',psite)]]
      sample$timefed<-sample[[paste0('timefed',psite)]]
      sample$timefedControlsComplete<-sample[[paste0('timefedControlsComplete',psite)]]
      sample$treatexCC<-ifelse(sample$treatControlsComplete==1 &sample$treatst==0,1,0)
    }
    if(treatc=='MUATE'){
      #Municipal ATE
      sample$treatdgwMU<- sample[[paste0('treatd',llc,'gw',psite)]]* sample$WaterStndCode.fMU
      sample$treatgroupMU<-sample$treatmentgroup * sample$WaterStndCode.fMU
      sample$controlMU<-sample$control*sample$WaterStndCode.fMU
      
      sample$sample.MUATE<-sample$control+sample$treatgroupMU
      
      sample<-subset(sample, sample.MUATE==1)
      timefe<-dplyr::select(sample, starts_with('timefe'))
      treatgroupm<-dplyr::select(sample, starts_with('treatmentgroup'))
      year<-dplyr::select(sample, starts_with('year'))
      bin<-dplyr::select(sample, starts_with('bin'))
      sample$demlogprice<-demeanlist(sample$logprice,
                                     list(as.factor(sample$PropertyAddressCensusTractAndBlock)))
      sample$aftpropnpl<-sample[[paste0('aftpropnpl',psite)]]
      
      #sample$aftfinalnpl<-sample[[paste0('aftfinalnpl',fnplsite)]]
      #sample$timefinalnplfe<-sample[[paste0('timefinalnplfe',fnplsite)]]
      sample$treatControlsComplete<-sample[[paste0('treatControlsComplete',psite)]]
      sample$timefed<-sample[[paste0('timefed',psite)]]
      sample$timefedControlsComplete<-sample[[paste0('timefedControlsComplete',psite)]]
      #sample$treatst<-sample[[paste0('treatd',llc,'gw',psite)]] #*sample[[paste0('dist',dist[[5]])]]
      sample$treatexCC<-ifelse(sample$treatControlsComplete==1 &sample$treatst==0,1,0)
    }
    if(treatc=='WLATE'){
      #Municipal ATE
      sample$treatdgwWL<- sample[[paste0('treatd',llc,'gw',psite)]]* sample$WaterStndCode.fWL
      sample$treatgroupWL<-sample$treatmentgroup * sample$WaterStndCode.fWL
      sample$controlWL<-sample$control*sample$WaterStndCode.fWL
      
      sample$sample.WLATE<-sample$control+sample$treatgroupWL
      
      sample<-subset(sample, sample.WLATE==1)
      timefe<-dplyr::select(sample, starts_with('timefe'))
      treatgroupm<-dplyr::select(sample, starts_with('treatmentgroup'))
      year<-dplyr::select(sample, starts_with('year'))
      bin<-dplyr::select(sample, starts_with('bin'))
      sample$demlogprice<-demeanlist(sample$logprice,
                                     list(as.factor(sample$PropertyAddressCensusTractAndBlock)))
      sample$aftpropnpl<-sample[[paste0('aftpropnpl',psite)]]
      
      #sample$aftfinalnpl<-sample[[paste0('aftfinalnpl',fnplsite)]]
      #sample$timefinalnplfe<-sample[[paste0('timefinalnplfe',fnplsite)]]
      sample$treatControlsComplete<-sample[[paste0('treatControlsComplete',psite)]]
      sample$timefed<-sample[[paste0('timefed',psite)]]
      sample$timefedControlsComplete<-sample[[paste0('timefedControlsComplete',psite)]]
      #sample$treatst<-sample[[paste0('treatd',llc,'gw',psite)]] #*sample[[paste0('dist',dist[[5]])]]
      sample$treatexCC<-ifelse(sample$treatControlsComplete==1 &sample$treatst==0,1,0)
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
      
      
      xTATE<-model.matrix(~ LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                            sqfeet+day+prediffdate+predate+prelogprice+
                            data.matrix(year[,4:25])
                          -1,sample)
      xcTATE<-model.matrix(~(data.matrix(bin)+
                               poly(LotSizeSquareFeet,4) + poly(YearBuilt,4) +   
                               poly(sqfeet,4))^3-1,sample)
      
      sample$tprelprice<-sample$prelogprice*ifelse(sample$treatst==1 & sample$presstatusd==0,1,0)
      sample$treatfin<-ifelse(sample$treatst==1 & sample$presstatusd==0,1,0)
      indx<- factor(as.numeric(cut2(sample$preprice,cuts=qcut,minmax=TRUE)))
                    
      treatind<-model.matrix(~treatst:indx-1,sample)
      treatind<-treatind[,1:(quant)]
      
      X<-model.matrix(~ treatmentgroup+#:indx+#bs(timetotreat,5)+bs(day,5)+#as.factor(round(timetotreat,1))+
                        #poly(timetotreat,3)+poly(day,3)+
                        #treatmentgroup*timetotreat+#treatmentgroup*day+
                        LotSizeSquareFeet + YearBuilt + FullBath + HalfBath + 
                        sqfeet+ prediffdate+predate+prelogprice+ 
                        indx+#presstatusd+
                        #as.factor(floor(timetotreat))+
                        data.matrix(year[,4:25])-1,sample)
      
      qr.X <- qr(X, tol=1e-3, LAPACK = FALSE)
      (rnkX <- qr.X$rank)  ## 4 (number of non-collinear columns)
      (keep <- qr.X$pivot[seq_len(rnkX)])
      ## 1 2 4 5 
      X <- X[,keep]
      
      results.lm.t.did<-felm(logprice ~treatind+X|as.factor(PropertyAddressCensusTractAndBlock),sample) #as.factor(HHID)+
      summary(results.lm.t.did)
        
      resid.lm.t.did<-felm(logprice ~X[,-1]|as.factor(PropertyAddressCensusTractAndBlock),sample) #as.factor(HHID)+
      summary(resid.lm.t.did)
      
      coefdid<-as.matrix(resid.lm.t.did$coefficients)
      coefdid[is.na(coefdid),]<-0
      #Xp<-X
      #Xp[,1]<-0
      #sample$crdid<-sample$demlogprice-cbind(treatind,X)%*%coefdid
      sample$crdid<-sample$demlogprice-X[,-1]%*%coefdid
      sample$crdid<-sample$crdid-mean(sample$crdid)
      sample$post<-ifelse(sample$timetotreat>0,1,0)
      treatment <- aggregate(crdid ~ treatmentgroup+as.factor(round(timetotreat,2)), data=sample, FUN=mean, na.rm=TRUE)
      names(treatment)[2]<-"timetotreat"
      treatment$timetotreat<-as.numeric(as.character(treatment$timetotreat))
      #treatment<-treatment[as.numeric(treatment$timetotreat)!=0.5,]
     
      ggplot() +
        geom_point(data=subset(treatment,treatmentgroup==1), aes(x=timetotreat, y=crdid, color= "Treatment Group")) +
        geom_point(data=subset(treatment,treatmentgroup==0), aes(x=timetotreat, y=crdid, color= "Control Group")) +
        stat_smooth(method = 'loess', formula = y ~ x  ,data=subset(sample,treatmentgroup==1&post==1),se= FALSE,  aes(x=timetotreat, y=crdid),color= "blue") +
        stat_smooth(method = 'loess', formula = y ~ x  ,data=subset(sample,treatmentgroup==1&post==0),se= FALSE,  aes(x=timetotreat, y=crdid),color= "blue") +
        
        stat_smooth(method = 'loess', formula = y ~ x  ,data=subset(sample,treatmentgroup==0&post==1),se= FALSE,  aes(x=timetotreat, y=crdid),color= "red") +
        stat_smooth(method = 'loess', formula = y ~ x  ,data=subset(sample,treatmentgroup==0&post==0),se= FALSE,  aes(x=timetotreat, y=crdid),color= "red") +
        ggtitle("Common Trends Assumption")+
        xlab('Years from Deletion') + labs(color="Legend") +  geom_vline(xintercept=0)+
        ylab('Monthly Average Residuals')
      
      ggsave(file=paste(path,'latex/','lmdidavg',treatc,dic, '.png', sep=""),height = 6,width =10)
      
      results.lm.t.did.agg<-felm(logprice ~treatst+X|as.factor(PropertyAddressCensusTractAndBlock),sample) #as.factor(HHID)+
      summary(results.lm.t.did.agg)
      
      betas.lm.did[treat,di]<-as.numeric(coef(summary(results.lm.t.did.agg))[,"Estimate"][1])
      ses.lm.did[treat,di]<-as.numeric(coef(summary(results.lm.t.did.agg))[,"Std. Error"][1])
      ps.lm.did[treat,di]<-as.numeric(coef(summary(results.lm.t.did.agg))[,"Pr(>|t|)"][1])
      
      results.gam<-mgcv::gam(logprice~treatst+X+#s(prelogprice,bs="cr")+
                               s(day,bs="gp")+s(lat,long,bs="tp",m=3,k=300),data=sample)
      gam.model<-mgcv::summary.gam(results.gam)   
      
      betas.gam.did[treat,di]<-as.numeric(gam.model$p.table[,"Estimate"])[2]
      ses.gam.did[treat,di]<-as.numeric(gam.model$p.table[,"Std. Error"])[2]
      ps.gam.did[treat,di]<-as.numeric(gam.model$p.table[,"Pr(>|t|)"])[2]
      
      
      
      
      #results.spline<-lm(logprice~treatind+X+spTATE,sample)
      #summary(results.spline)
      
      #results.semip<-semip(logprice~treatind+X[,1:8]+preprice,nonpar=~lat+long,window1 = .5, window2 = .5,
            #               kern="tcub",distance="Mahal",targetfull=NULL, print.summary=TRUE, data=sample)
      #summary(results.semip)
      #library(gam)
      library(mgcv)
      s=mgcv:::s
      
      if(treatc=='TATE' & dic=='6k'){
      results.gam<-mgcv::gam(logprice~treatind+X+#s(prelogprice,bs="cr")+
                               s(day,bs="gp")+s(lat,long,bs="tp",m=3,k=300),data=sample)
      gam.model<-mgcv::summary.gam(results.gam)   
      
      
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
      
      
      if(treatc=="TATE"){
        betas.lm.t.did.TATE[,di]<-as.numeric(coef(summary(results.lm.t.did))[,"Estimate"][1:quant])
        betas.gam.t.did.TATE[,di]<-as.numeric(gam.model$p.table[,"Estimate"])[2:(1+quant)]
        #betas.lm.t.es.TATE[,di]<-as.numeric(coef(summary(results.lm.t.es))[,"Estimate"][1:quant])
        ses.lm.t.did.TATE[,di]<-as.numeric(coef(summary(results.lm.t.did))[,"Std. Error"][1:quant])
        ses.gam.t.did.TATE[,di]<-as.numeric(gam.model$p.table[,"Std. Error"])[2:(1+quant)]
        #ses.lm.t.es.TATE[,di]<-as.numeric(coef(summary(results.lm.t.es))[,"Std. Error"][1:quant])
        ps.lm.t.did.TATE[,di]<-as.numeric(coef(summary(results.lm.t.did))[,"Pr(>|t|)"][1:quant])
        ps.gam.t.did.TATE[,di]<-as.numeric(gam.model$p.table[,"Pr(>|t|)"])[2:(1+quant)]
        #ps.lm.t.es.TATE[,di]<-as.numeric(coef(summary(results.lm.t.es))[,"Pr(>|t|)"][1:quant])
        }
        
        if(treatc=="MUATE"){
          betas.lm.t.did.MUATE[,di]<-as.numeric(coef(summary(results.lm.t.did))[,"Estimate"][1:quant])
          betas.gam.t.did.MUATE[,di]<-as.numeric(gam.model$p.table[,"Estimate"])[2:(1+quant)]
          #betas.lm.t.es.TATE[,di]<-as.numeric(coef(summary(results.lm.t.es))[,"Estimate"][1:quant])
          ses.lm.t.did.MUATE[,di]<-as.numeric(coef(summary(results.lm.t.did))[,"Std. Error"][1:quant])
          ses.gam.t.did.MUATE[,di]<-as.numeric(gam.model$p.table[,"Std. Error"])[2:(1+quant)]
          #ses.lm.t.es.TATE[,di]<-as.numeric(coef(summary(results.lm.t.es))[,"Std. Error"][1:quant])
          ps.lm.t.did.MUATE[,di]<-as.numeric(coef(summary(results.lm.t.did))[,"Pr(>|t|)"][1:quant])
          ps.gam.t.did.MUATE[,di]<-as.numeric(gam.model$p.table[,"Pr(>|t|)"])[2:(1+quant)]
          #ps.lm.t.es.TATE[,di]<-as.numeric(coef(summary(results.lm.t.es))[,"Pr(>|t|)"][1:quant])
        }
        
        if(treatc=="WLATE"){
          betas.lm.t.did.WLATE[,di]<-as.numeric(coef(summary(results.lm.t.did))[,"Estimate"][1:quant])
          betas.gam.t.did.WLATE[,di]<-as.numeric(gam.model$p.table[,"Estimate"])[2:(1+quant)]
          #betas.lm.t.es.TATE[,di]<-as.numeric(coef(summary(results.lm.t.es))[,"Estimate"][1:quant])
          ses.lm.t.did.WLATE[,di]<-as.numeric(coef(summary(results.lm.t.did))[,"Std. Error"][1:quant])
          ses.gam.t.did.WLATE[,di]<-as.numeric(gam.model$p.table[,"Std. Error"])[2:(1+quant)]
          #ses.lm.t.es.TATE[,di]<-as.numeric(coef(summary(results.lm.t.es))[,"Std. Error"][1:quant])
          ps.lm.t.did.WLATE[,di]<-as.numeric(coef(summary(results.lm.t.did))[,"Pr(>|t|)"][1:quant])
          ps.gam.t.did.WLATE[,di]<-as.numeric(gam.model$p.table[,"Pr(>|t|)"])[2:(1+quant)]
          #ps.lm.t.es.TATE[,di]<-as.numeric(coef(summary(results.lm.t.es))[,"Pr(>|t|)"][1:quant])
        }
      #assign(paste0('lm.TATE.',k,i),lm.TATE)
      #summary(results.lm.t.did)
      #summary(results.lm.t.es)
      
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
      print(paste0('site = ',psite))
    }
  }
  }
}
for(statchange in c('')){
  for(meth in c('lm')){
    for(inf in c('did')){
      for(treat in  treatl){
      #treat<-"TATE"
        
      p<-get(paste0(statchange,'ps.',meth,'.t.',inf,'.',treat))
      mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ", ifelse(p < .1, "^\\bullet  ", " "))))
      
      #pb<-exp(get(paste0(statchange,'betas.',meth,'.t.',inf,'.',treat)))-1
      #rpb<-round(pb,3)
      #se<-round(exp(get(paste0(statchange,'ses.',meth,'.t.',inf,'.',treat)))-1,3)
      
      pb<-get(paste0(statchange,'betas.',meth,'.t.',inf,'.',treat))
      rpb<-round(pb,3)
      se<-round(get(paste0(statchange,'ses.',meth,'.t.',inf,'.',treat)),3)
      
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
      
      colnames(results.mat)<-c('10k','8k','6k','4k','2k')
     # rn<-c(paste0("(",qcut[1],","),paste0(qcut[2],"]"),paste0("(",qcut[2],","),paste0(qcut[3],"]"),
      #      paste0("(",qcut[3],","),paste0(qcut[4],"]"),
       #     paste0("(",qcut[4],","),paste0(qcut[5],"]"),paste0("(",qcut[5],","),
        #    paste0(qcut[6],"]"),paste0("(",qcut[6],","),paste0(qcut[7],"]"),paste0("(",qcut[7],","),
         #   paste0(qcut[8],"]"),paste0("(",qcut[8],","),paste0(qcut[9],"]"),paste0("(",qcut[9],","),
          #  paste0(qcut[10],"]"),
           # paste0("(",qcut[10],","),paste0(qcut[11],"]"))
      rn<-c(paste0("(",qcut[1],",",qcut[2],"]"),"1" ,paste0("(",qcut[2],",",qcut[3],"]"),"2 " ,
            paste0("(",qcut[3],",",qcut[4],"]"),"3"  ,
            paste0("(",qcut[4],",",qcut[5],"]"), "4",paste0("(",qcut[5],",",qcut[6],"]"),"5" ,
            paste0("(",qcut[6],",",qcut[7],"]"),"6" ,paste0("(",qcut[7],",",qcut[8],"]"),"7" , 
            paste0("(",qcut[8],",",qcut[9],"]"), "8 ",paste0("(",qcut[9],",",qcut[10],"]"), "9",
            paste0("(",qcut[10],",",qcut[11],"]"), "10")
      rn2<-c(paste0("(",qcut[1],",",qcut[2],"]"),paste0("(",qcut[2],",",qcut[3],"]"),
            paste0("(",qcut[3],",",qcut[4],"]"),
            paste0("(",qcut[4],",",qcut[5],"]"),paste0("(",qcut[5],",",qcut[6],"]"),
            paste0("(",qcut[6],",",qcut[7],"]"),paste0("(",qcut[7],",",qcut[8],"]"),
            paste0("(",qcut[8],",",qcut[9],"]"), paste0("(",qcut[9],",",qcut[10],"]"),
            paste0("(",qcut[10],",",qcut[11],"]"))
      rownames(results.mat)<-rn
      xtable(results.mat)
      print.xtable(xtable(results.mat),include.rownames=TRUE, 
                   include.colnames=TRUE, sanitize.text.function = identity,
                   type="latex", file=paste0(path,'latex/',meth,inf,statchange,treat,"buffer",buf,".tex"))
     
      allModelFrame <- data.frame(Variable = 1:10,
                               Coefficient = pb[,1],
                               SE = se[, 1],
                               modelName = "10k")
      for(i in 2:length(dist)){
      di<-dist[i]
      modelFrame <- data.frame(Variable =  1:10,
                                Coefficient = pb[,i],
                                SE = se[, i],
                                modelName = di)
      allModelFrame <- data.frame(rbind(allModelFrame,modelFrame))
      
    }
      qu<-c("1","2","3","4","5","6","7","8","9","10")
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
        zp1 <- zp1 + scale_x_continuous(breaks=seq(0,10,1))
      print(zp1)  # The trick to these is position_dodge().
      
      ggsave(file=paste(path,'latex/','coeff',treat, '.png', sep=""),height = 7,width =9)
      
      }
  }
}
}

        
        p<-ps.lm.did
        mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ", ifelse(p < .1, "^\\bullet  ", " "))))
        
        #pb<-exp(get(paste0(statchange,'betas.',meth,'.t.',inf,'.',treat)))-1
        #rpb<-round(pb,3)
        #se<-round(exp(get(paste0(statchange,'ses.',meth,'.t.',inf,'.',treat)))-1,3)
        
        pb<-betas.lm.did
        rpb<-round(pb,3)
        se<-round(ses.lm.did,3)
        
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
        
        colnames(results.mat)<-c('10k','8k','6k','4k','2k')
        # rn<-c(paste0("(",qcut[1],","),paste0(qcut[2],"]"),paste0("(",qcut[2],","),paste0(qcut[3],"]"),
        #      paste0("(",qcut[3],","),paste0(qcut[4],"]"),
        #     paste0("(",qcut[4],","),paste0(qcut[5],"]"),paste0("(",qcut[5],","),
        #    paste0(qcut[6],"]"),paste0("(",qcut[6],","),paste0(qcut[7],"]"),paste0("(",qcut[7],","),
        #   paste0(qcut[8],"]"),paste0("(",qcut[8],","),paste0(qcut[9],"]"),paste0("(",qcut[9],","),
        #  paste0(qcut[10],"]"),
        # paste0("(",qcut[10],","),paste0(qcut[11],"]"))
        rn<-c("Total","  ","Municipal Water","   ","Well Water","    ")
       
        rownames(results.mat)<-rn
        xtable(results.mat)
        print.xtable(xtable(results.mat),include.rownames=TRUE, 
                     include.colnames=TRUE, sanitize.text.function = identity,
                     type="latex", file=paste0(path,'latex/ATEcomp.tex'))
        
        allModelFrame <- data.frame(Variable = 1:10,
                                    Coefficient = pb[,1],
                                    SE = se[, 1],
                                    modelName = "10k")
        for(i in 2:length(dist)){
          di<-dist[i]
          modelFrame <- data.frame(Variable =  1:10,
                                   Coefficient = pb[,i],
                                   SE = se[, i],
                                   modelName = di)
          allModelFrame <- data.frame(rbind(allModelFrame,modelFrame))
          
        }
        qu<-c("1","2","3","4","5","6","7","8","9","10")
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
        zp1 <- zp1 + scale_x_continuous(breaks=seq(0,10,1))
        print(zp1)  # The trick to these is position_dodge().
        
        ggsave(file=paste(path,'latex/','coeff',treat, '.png', sep=""),height = 7,width =9)


