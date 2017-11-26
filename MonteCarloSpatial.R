
## Preliminaries
rm(list=ls())

path<- "P:/Peter/Hedonics/"

## This function will check if a package is installed, and if not, install it
pkgTest <- function(x) {
  if (!require(x, character.only = TRUE))
  {
    install.packages(x, dep = TRUE)
    if(!require(x, character.only = TRUE)) stop("Package not found")
  }
}

## These lines load the required packages
packages <- c("readxl","MASS","rpgm","splines","Hmisc","grf","hdi","sphet","spdep","tmle","xgboost", "data.table","magrittr","Hmisc","DT","h2o","nnet","caret","quantmod","neuralnet","doSNOW","parallel","compare","doParallel","devtools","foreach","spdep","reshape2","sm","plyr","utils","tcltk","geosphere", "matrixcalc", "dplyr","ExPosition", "randomForest","lfe", "hdm", "rdrobust", "stargazer", "ggplot2", "outliers","rpart","e1071")
lapply(packages, pkgTest)


## These lines set several options
options(scipen = 999) # Do not print scientific notation
options(stringsAsFactors = FALSE) ## Do not load strings as factors

#install.packages("mvtnorm")
library(mvtnorm)
# Setup parallel computation - use all cores on our computer.
# (Install "parallel" and "RhpcBLASctl" if you don't already have those packages.)
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

# different configurations.
bart.tune <- list(num_burn_in = c(100),
                  alpha = c(0.95,0.8,0.5), beta = c(1,2,3),
                  num_iterations_after_burn_in = c(500,1000))

# Set detailed names = T so we can see the configuration for each function.
# Also shorten the name prefix.
bart <- create.Learner("SL.bartMachine", tune = bart.tune, detailed_names = T, name_prefix = "bart")

# configurations
length(bart$names)
bart$names

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
#SL.library$add("SL.ksvm")
#SL.library$add("SL.polymars")
#SL.library$add("SL.kernelKnn")
#SL.library$add("SL.randomForest")
SL.library$add("SL.xgboost")
#SL.library$add(c("SL.nnet","screen.glmnet"))
#SL.library$add("SL.caret.rpart")
#SL.library$add("SL.caret" )
#SL.library$add("SL.knn")
#SL.library$add("SL.bartMachine")
#SL.library$add("SL.dbarts")
#SL.library$add("SL.loess" )

for(i in 1:length(glmnet$names)){
  SL.library$add(glmnet$names[i])
}    

#for(i in 1:length(xgboost$names)){
#  SL.library$add(c(xgboost$names[i],"screen.glmnet"))
#}

#for(i in 1:length(randomForest$names)){
#  SL.library$add(c(randomForest$names[i],"screen.glmnet"))
#}

#for(i in 1:length(bart$names)){
#  SL.library$add(c(bart$names[i],"screen.glmnet"))
#}

for(i in 1:length(xgboost$names)){
  SL.library$add(c(xgboost$names[i]))
}

#for(i in 1:length(randomForest$names)){
#  SL.library$add(c(randomForest$names[i]))
#}

#for(i in 1:length(bart$names)){
#  SL.library$add(c(bart$names[i]))
#}
SL.library$as.list()

set.seed(12345)
#set n=200 and number of replications=200
N <- 2000
P<-40
n.sims <- 20



#creates matrix for the outputs for the tests from 2a
betahat.tmle <- matrix(NA,nrow=n.sims,ncol=3)
betahat.ds <- matrix(NA,nrow=n.sims,ncol=1)
betahat.ols.os <- matrix(NA,nrow=n.sims,ncol=3)
betahat.ds.os <- matrix(NA,nrow=n.sims,ncol=1)
type1.tmle <- matrix(NA,nrow=n.sims,ncol=1)
type1.ds <- matrix(NA,nrow=n.sims,ncol=1)
type1.ols.os <- matrix(NA,nrow=n.sims,ncol=1)
type1.ds.os <- matrix(NA,nrow=n.sims,ncol=1)
type2.tmle <- matrix(NA,nrow=n.sims,ncol=1)
type2.ds <- matrix(NA,nrow=n.sims,ncol=1)
type2.ols.os <- matrix(NA,nrow=n.sims,ncol=1)
type2.ds.os <- matrix(NA,nrow=n.sims,ncol=1)

#creates loops for each value of rho, each test of beta (rows in betam)

for(j in 1:n.sims){
  # Create data  
  x <- matrix(runif(N*P), nrow = N, ncol = P)
  nvar<-P
  quants<-5
  dumx<-matrix(nrow=N,ncol=quants*nvar)
  for(i in 1:nvar){
    indx<- factor(as.numeric(cut2(x[,i], g=quants)))
    dumx[,((quants*i)-(quants-1)):(quants*i)]<-model.matrix(~indx-1)
  }
  
  ps<-x[,11] * 2+ x[,12] * 2.5 + 60*x[,13]-2*x[,13]^2-2*x[,13]^3+
    10*dumx[,(quants*15)]+9*dumx[,(quants*15)+1]+3*dumx[,(quants*15)+2]+
    10*dumx[,(quants*18)+1]+5*dumx[,(quants*18)+2]+7*dumx[,(quants*18)+3]+
    30*rnorm(N)
  #Treatment<- ifelse(ps-median(ps)>0,1,0) 
  polyx<-cbind(x^2,x^3)
  
  I<-diag(N)
  alpha<-.75
  lat<-.08*runif(N)-87
  long<-1*runif(N)+41
  time<-1990+(10*runif(N))
  psn<-rnorm(N)
  Treatment<- ifelse(psn-median(psn)>0,1,0) 
  #Treatment<- ifelse(time-median(time)>0,1,0) 
  
  truegridquants<-15
  indlat<- factor(as.numeric(cut2(lat, g=modelgridquants)))
  indlong<- factor(as.numeric(cut2(long, g=modelgridquants)))
  blockfe<-model.matrix(~indmlat:indmlong-1)
  
  modelgridquants<-10
  indmlat<- factor(as.numeric(cut2(lat, g=modelgridquants)))
  indmlong<- factor(as.numeric(cut2(long, g=modelgridquants)))
  blockmfe<-model.matrix(~indmlat:indmlong-1)
  
  treatlat<-median(lat)
  treatlong<-median(long)
  disttreat<-distm(cbind(long, lat),c(treatlong,treatlat), fun = distHaversine)
  
  nslat<-max(lat)
  nslong<-max(long)
  distns<-distm(cbind(long, lat),c(nslong,nslat), fun = distHaversine)
  
  TEi<-pmin(100000*(1/(disttreat+2000)),35*rep(1,N))*Treatment
  TE<-ifelse(TEi<20,0,TEi)
  summary(TE)
  #Create a function to generate a continuous color palette
  rbPal <- colorRampPalette(c('red','blue'))
  
  #This adds a column of color values
  # based on the y values
  Col <- rbPal(10)[as.numeric(cut(TE,breaks = 10))]
  
  plot(long,lat,pch = 20,col = Col)
  legend("topleft",title="Decile",legend=levels(cut(TE,breaks = 10)),col =rbPal(10),pch=20)
  
  TEi<-pmin(100000*(1/(disttreat+2000)),35*rep(1,N))
  TE2<-ifelse(TEi<20,0,TEi)
  summary(TE2)
  #Create a function to generate a continuous color palette
  rbPal <- colorRampPalette(c('red','blue'))
  
  #This adds a column of color values
  # based on the y values
  Col <- rbPal(10)[as.numeric(cut(TE2,breaks = 10))]
  
  plot(long,lat,pch = 20,col = Col)
  legend("topleft",title="Decile",legend=levels(cut(TE2,breaks = 10)),col =rbPal(10),pch=20)
  
  
  dW<-distm(cbind(long, lat), fun = distHaversine)
  dW<-ifelse(dW<2000, dW, -10)
  Ws<-1/(1+dW)
  Ws<-ifelse(Ws>0, Ws, 0)
  
  rep.row<-function(x,n){
    matrix(rep(x,each=n),nrow=n)
  }
  t<-time
  T1<-rep.row(t,N)
  T2<-T1
  T1<-t(T1)
  Tdp<-T1-T2
  Tdp[Tdp >= 5*365] = -2
  Tdp[Tdp <=0] = -2
  Wtp<-1/(1+Tdp)
  Wtp[Wtp<0]<-0
  Wf<-Wtp*Ws
  W<-ifelse(Wf==0,0,Wf/rowSums(Wf))
  
  # Create y  
  prey <- x[,1] * 2+ x[,2] * 2.5 + 60*x[,3]-2*x[,3]^2-2*x[,3]^3+
    #5*(x[,4]*x[,5])+2*(x[,4]*x[,5])^2+
    6*x[,7]^2+3*x[,8]^2- 2*x[,8]^3+
    10*dumx[,(quants*5)]+9*dumx[,(quants*5)+1]+3*dumx[,(quants*5)+2]+
    10*dumx[,(quants*8)+1]+5*dumx[,(quants*8)+2]+7*dumx[,(quants*8)+3]+
    1*dumx[,(quants*9)]+6*dumx[,(quants*9)+3]-6*dumx[,(quants*9)+5]+
    .8*time+.000002*time^2+.0000002*time^3
  
  y <- TE+5*Treatment*rnorm(N) +.5*Treatment*(time-1989)^2+
    #(10+50*x[,1])*Treatment+
    x[,2] * 2.5 + 60*x[,3]-2*x[,3]^2-2*x[,3]^3+
    #5*(x[,4]*x[,5])+2*(x[,4]*x[,5])^2+
    6*x[,7]^2+3*x[,8]^2- 2*x[,8]^3+
    10*dumx[,(quants*5)]+9*dumx[,(quants*5)+1]+3*dumx[,(quants*5)+2]+
    10*dumx[,(quants*8)+1]+5*dumx[,(quants*8)+2]+7*dumx[,(quants*8)+3]+
    1*dumx[,(quants*9)]+6*dumx[,(quants*9)+3]-6*dumx[,(quants*9)+5]+
    -1.6*time+.000004*time^2+.0000004*time^3+(90000/(distns+2000))+
      (90/(distns+2000))*time+
    rnorm(N)
  #y<-ginv(I-alpha*W)%*%prey
  
  #### Estimating true relationships
  
  #set points  
  at <- as.matrix(seq(min(lat), max(lat), length.out = 10))
  bt <- as.matrix(seq(min(long), max(long), length.out = 10))
  setp<-matrix(ncol=2,nrow =110 )
  for(i in 1:dim(at)[1]){
    for(j in 1:dim(bt)[1]){
      setp[(i*dim(bt)[1])+j-1,1:2]<-cbind(at[i],bt[j])
    }
  }
  setp<-setp[10:109,]
  #dsetp<-distm( cbind(long, lat),setp,fun = distHaversine)
  ## Overspecify
  sdf<-20
  splat<-bs(lat, df = sdf)
  splong<-bs(long, df = sdf)
  spint<-model.matrix(~splat:splong)
  
  xds<-cbind(lat,long,lat,long,lat,long,
             lat,long,lat,long,lat,long,
             #lat+long,lat*long,lat^2*long^2,lat/long,
             #1/(lat*long),lat^2+long^2,
             blockmfe,time,bs(time,df =sdf),x)
  feds<-cbind(model.matrix(~as.factor(floor(time))))
  xtl<-cbind(x,W%*%y,lat,long,time)
  xgrf<-cbind(xds,feds,bs(disttreat,df =sdf),spint,splat,splong,disttreat)
  
  #xds2<-cbind(lat,long, blockmfe,time,x,bs(time, df = sdf),polyx,dumx )
  #feds<-cbind(model.matrix(~as.factor(floor(time))))
 # xgrf<-cbind(xds2,feds,spint,splat,splong,disttreat*Treatment,bs(disttreat,df =sdf)*Treatment,disttreat,bs(disttreat,df =sdf))
  
  ###########################################################################################
  #LASSO
  #glmfit<-glmnet()
  
  ############################################################################################
  
  #GRF
  c.forest<-causal_forest(X=xgrf, Y=y, W=Treatment,
                          num.trees=10000,#min.node.size=10,num.threads =100,
                          mtry = ceiling(ncol(xgrf)*.3) )
  
  xtest<- matrix(0, 10100, dim(xgrf)[2])
  a <- as.matrix(seq(min(long), max(long), length.out = 100))
  b <- as.matrix(seq(min(lat), max(lat), length.out = 100))
  
  for(i in 1:dim(a)[1]){
    for(j in 1:dim(b)[1]){
      xtest[(i*dim(b)[1])+j,1:2]<-cbind(a[i],b[j])
    }
  }
  xtest<-xtest[101:dim(xtest)[1],]
  
  xtest[,3]<-xtest[,1]
  xtest[,4]<-xtest[,2]
  xtest[,5]<-xtest[,1]
  xtest[,6]<-xtest[,2]
  xtest[,5]<-xtest[,1]
  xtest[,6]<-xtest[,2]
  xtest[,7]<-xtest[,1]
  xtest[,8]<-xtest[,2]
  xtest[,9]<-xtest[,1]
  xtest[,10]<-xtest[,2]
  xtest[,11]<-xtest[,1]
  xtest[,12]<-xtest[,2]
  
  
  indmlat<- factor(as.numeric(cut2(xtest[,1], g=modelgridquants)))
  indmlong<- factor(as.numeric(cut2(xtest[,2], g=modelgridquants)))
  xtest[,13:(dim(blockmfe)[2]+12)]<-model.matrix(~indmlat:indmlong-1)
  
  xtest1<- matrix(0, 100, dim(xgrf)[2])
  xtest1[,(dim(blockmfe)[2]+13)] <- as.matrix(seq(min(time), max(time), length.out = 100))
  xtest1[,(dim(xds)[2]+1):(dim(xds)[2]+dim(feds)[2])]<-model.matrix(~as.factor(floor(xtest1[,(dim(blockmfe)[2]+13)] )))
  xtest1[,(dim(blockmfe)[2]+14):(dim(blockmfe)[2]+13+dim(bs(time,df =sdf))[2])]<-bs(xtest1[,(dim(blockmfe)[2]+13)],df =sdf)
  xtest<-rbind(xtest,xtest1)
  
  xtest[,dim(xgrf)[2]]<-distm(cbind(xtest[,1], xtest[,2]),c(treatlong,treatlat), fun = distHaversine)
  
  #dsetptest<-distm( cbind(b, a),setp,fun = distHaversine)
  #xtest[,(dim(blockmfe)[2]+13):((dim(blockmfe)[2]+13)+dim(dsetptest)[2])]<-dsetptest
  
  xtest[,(dim(xgrf)[2]-dim(splat)[2]-dim(splong)[2]):(dim(xgrf)[2]-1-dim(splong)[2])]<-bs(xtest[,2], df = sdf)
  xtest[,(dim(xgrf)[2]-dim(splong)[2]):(dim(xgrf)[2]-1)]<-bs(xtest[,1], df = sdf)
  int<-(dim(xgrf)[2]-dim(splat)[2]-dim(splong)[2]-dim(spint)[2]):(dim(xgrf)[2]-1-dim(splong)[2]-dim(splat)[2])
  
  xtest[,int] <-model.matrix(~bs(xtest[,2], df = sdf):bs(xtest[,1], df = sdf))
  
  disttn<-(dim(xgrf)[2]-dim(splat)[2]-dim(splong)[2]-dim(spint)[2]-dim(bs(disttreat,df =sdf))[2]):(dim(xgrf)[2]-dim(splat)[2]-dim(splong)[2]-dim(spint)[2]-1)
  
  xtest[,disttn]<-bs(xtest[,dim(xtest)[2]],df =sdf)
  
  if(FALSE){
  for(i in 1:dim(xtest)[2]){
    if(sum(xtest[,i])==0){
      xtest[,i]<-median(xgrf[,i])
    }
  }
  }
  
  tau.hat <- predict(c.forest, xtest, estimate.variance = TRUE)
  t<-tau.hat$predictions/sqrt(tau.hat$variance.estimates/N)
  p<-ifelse(abs(t)>1.96,1,0)
  
  #Create a function to generate a continuous color palette
  rbPal <- colorRampPalette(c('red','blue'))
  
  #This adds a column of color values
  # based on the y values
  Col <- rbPal(10)[as.numeric(cut(tau.hat$predictions[1:9999],breaks = 10))]
  
  plot(xtest[1:9999,1],xtest[1:9999,2],pch = 20,col = Col)
  legend("topleft",title="Decile",legend=levels(cut(tau.hat$predictions[1:9999],breaks = 10)),col =rbPal(10),pch=20)
  
  #Create a function to generate a continuous color palette
  rbPal2 <- colorRampPalette(c('yellow','green'))
  
  #This adds a column of color values
  # based on the y values
  Col2 <- rbPal2(2)[as.numeric(cut(p[1:9999],breaks = 2))]
  
  
  plot(xtest[1:9999,1],xtest[1:9999,2],pch = 20,col = Col2)
  legend("topleft",title="Significance",legend=c("p>0.05","p<0.05"),col =rbPal2(2),pch=20)
  
  rbPal3 <- colorRampPalette(c('blue','red'))
  
  #This adds a column of color values
  # based on the y values
  pos<-ifelse(t>1.96,1,0)
  neg<-ifelse(t<(-1.96),-1,0)
  qual<-pos+neg
  Col3 <- rbPal3(3)[as.numeric(cut(qual[1:10000],breaks = 3))]
  
  plot(xtest[1:9999,1],xtest[1:9999,2],pch = 20,col = Col3)
  legend("top",title="Significance",legend=c("t< -1.96","-1.96<t<1.96","t>1.96"),col =rbPal3(3),pch=20)
  
  
  rbPal3 <- colorRampPalette(c('blue','red'))
  
  #This adds a column of color values
  # based on the y values
  
  plot(xtest[10001:10100,(dim(blockmfe)[2]+13)],tau.hat$predictions[10001:10100,],pch = 20,col = Col3)
  legend("top",title="Significance",legend=c("t< -1.96","-1.96<t<1.96","t>1.96"),col =rbPal3(3),pch=20)
  
  
  xtest[,dim(xgrf)[2]]<-distm(cbind(xtest[,1], xtest[,2]),c(treatlong,treatlat), fun = distHaversine)
  tTEi<-pmin(100000*(1/(xtest[1:9999,dim(xgrf)[2]]+2000)),35*rep(1,N))
  TTE<-ifelse(tTEi<20,0,tTEi)
  
  #TTE<-100000*(1/(xtest[1:9999,dim(xgrf)[2]]+2000))
  
  #Create a function to generate a continuous color palette
  rbPal <- colorRampPalette(c('red','blue'))
  
  #This adds a column of color values
  # based on the y values
  Col <- rbPal(10)[as.numeric(cut(tau.hat$predictions[1:9999]-TTE,breaks = 10))]
  
  plot(xtest[1:9999,1],xtest[1:9999,2],pch = 20,col = Col)
  legend("topleft",title="Decile",legend=levels(cut(tau.hat$predictions[1:9999]-TTE,breaks = 10)),col =rbPal(10),pch=20)
  
  grf.y.hat <- predict(c.forest, xgrf, estimate.variance = FALSE)
  
  #TE<-100000*Treatment*(1/(disttreat+2000))
  
  #Create a function to generate a continuous color palette
  rbPal <- colorRampPalette(c('red','blue'))
  
  #This adds a column of color values
  # based on the y values
  er<-grf.y.hat$predictions[Treatment==1,]-TTE[Treatment==1,]
  Col <- rbPal(10)[as.numeric(cut(er,breaks = 10))]
  
  plot(long,lat,pch = 20,col = Col)
  legend("topleft",title="Decile",legend=levels(cut(er,breaks = 10)),col =rbPal(10),pch=20)
  
  summary(er)
  
  #######################################################################
  #FE
  feds<-model.matrix(~as.factor(floor(time)-1))
  xfe<-model.matrix(~Treatment*blockmfe+x+feds-1)
  xfe<-cbind(Treatment,blockmfe,model.matrix(~Treatment:blockmfe-1),x,feds)#,disttreat,Treatment*disttreat)
  fit<-lm(y~xfe)
  summary(fit)
  
  xtest3<- matrix(0, 10000, dim(xfe)[2])
  a <- as.matrix(seq(min(long), max(long), length.out = 100))
  b <- as.matrix(seq(min(lat), max(lat), length.out = 100))
  xtest2<-matrix(0, 10100, dim(xfe)[2])
  for(i in 1:dim(a)[1]){
    for(j in 1:dim(b)[1]){
      xtest2[(i*dim(b)[1])+j,1:2]<-cbind(a[i],b[j])
    }
  }
  xtest2<-xtest2[101:dim(xtest2)[1],]
  
  
  indmlat<- factor(as.numeric(cut2(xtest2[,2], g=modelgridquants)))
  indmlong<- factor(as.numeric(cut2(xtest2[,1], g=modelgridquants)))
  xtest3[,2:(dim(blockmfe)[2]+1)]<-model.matrix(~indmlat:indmlong-1)
  xtest3[,(dim(blockmfe)[2]+2):(2*dim(blockmfe)[2]+1)]<-model.matrix(~indmlat:indmlong-1)
  xtest3[,(2*dim(blockmfe)[2]+1)]<-rep(0,10000)
    
  xtest3[,1]<-rep(1,10000)

  #xtest3[,dim(xfe)[2]-1]<-distm(cbind(xtest3[,1], xtest3[,2]),c(treatlong,treatlat), fun = distHaversine)
 # xtest3[,dim(xfe)[2]]<-xtest3[,dim(xfe)[2]-1]
  if(FALSE){
  for(i in 1:dim(xtest3)[2]){
    if(sum(xtest3[,i])==0){
      xtest3[,i]<-median(xfe[,i])
    }
  }
  }
  
  xtest3<-cbind(rep(1,10000),xtest3)
  fit$coefficients[is.na(as.numeric(coef(fit)))]<-0
  lm.pred1<-xtest3 %*% as.matrix(fit$coefficients )
  #lm.pred1<-predict(object=fit,newdata=data.frame(xtest3))
  xtest3[,1]<-rep(0,10000)
  xtest3[,dim(xfe)[2]]<-rep(0,10000)
  xtest3[,(dim(blockmfe)[2]+2):(2*dim(blockmfe)[2]+1)]<-rep(0,10000)
  #lm.pred0<-predict(object=fit,newdata=data.frame(xtest3))
  lm.pred0<-xtest3 %*% as.matrix(fit$coefficients )
  lm.te<-lm.pred1-lm.pred0
  
  
  xtest2[,dim(xfe)[2]]<-distm(cbind(xtest2[,1], xtest2[,2]),c(treatlong,treatlat), fun = distHaversine)
  
  tTEi<-pmin(100000*(1/(xtest2[,dim(xfe)[2]]+2000)),35*rep(1,N))
  TTE2<-ifelse(tTEi<20,0,tTEi)
  
  #TTE<-100000*(1/(xtest2[,dim(xfe)[2]]+2000))
  
  #Create a function to generate a continuous color palette
  rbPal <- colorRampPalette(c('red','blue'))
  
  #This adds a column of color values
  # based on the y values
  Col <- rbPal(10)[as.numeric(cut(lm.te,breaks = 10))]
  
  plot(xtest2[,1],xtest2[,2],pch = 20,col = Col)
  legend("topleft",title="Decile",legend=levels(cut(lm.te,breaks = 10)),col =rbPal(10),pch=20)
  
  
  #Create a function to generate a continuous color palette
  rbPal <- colorRampPalette(c('red','blue'))
  
  #This adds a column of color values
  # based on the y values
  Col <- rbPal(10)[as.numeric(cut(lm.te-TTE2,breaks = 10))]
  
  plot(xtest2[,1],xtest2[,2],pch = 20,col = Col)
  legend("topleft",title="Decile",legend=levels(cut(lm.te-TTE2,breaks = 10)),col =rbPal(10),pch=20)
  
  summary(lm.te-TTE2)
  summary(tau.hat$predictions[1:9999]-TTE)
  
  
  print(paste0(j, ' of ',n.sims))
}

