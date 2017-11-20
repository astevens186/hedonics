
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
packages <- c("readxl","MASS","grf","hdi","sphet","spdep","tmle","xgboost", "data.table","magrittr","Hmisc","DT","h2o","nnet","caret","quantmod","neuralnet","doSNOW","parallel","compare","doParallel","devtools","foreach","spdep","reshape2","sm","plyr","utils","tcltk","geosphere", "matrixcalc", "dplyr","ExPosition", "randomForest","lfe", "hdm", "rdrobust", "stargazer", "ggplot2", "outliers","rpart","e1071")
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
N <- 3000
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
  Treatment<- ifelse(ps-median(ps)>0,1,0) 
  polyx<-cbind(x^2,x^3)
  
  I<-diag(N)
  alpha<-.75
  lat<-.1*runif(N)-87
  long<-.1*runif(N)+41
  time<-1990+(10*runif(N))
  
  gridquants<-5
  indlat<- factor(as.numeric(cut2(lat, g=gridquants)))
  indlong<- factor(as.numeric(cut2(long, g=gridquants)))
  blockfe<-model.matrix(~indlat:indlong-1)
  gblockfe<-blockfe[,sample.int(n=(gridquants^2), size=((gridquants^2)/4))]%*%as.matrix(sample.int(n=50,size=((gridquants^2)/4)))
  
  modelgridquants<-25
  indmlat<- factor(as.numeric(cut2(lat, g=modelgridquants)))
  indmlong<- factor(as.numeric(cut2(long, g=modelgridquants)))
  blockmfe<-model.matrix(~indmlat:indmlong-1)
  
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
    gblockfe-.8*time+.000002*time^2+.0000002*time^3
  y <- -20*Treatment*blockfe[,7]-20*Treatment*blockfe[,19]+#500*Treatment*((lat-41.04)+(long+86.95))+
    20*Treatment*blockfe[,9]+20*Treatment*blockfe[,17]+20*Treatment*x[,4]^2+
    #5*Treatment+
    x[,1] * 2+ x[,2] * 2.5 + 60*x[,3]-2*x[,3]^2-2*x[,3]^3+
    #5*(x[,4]*x[,5])+2*(x[,4]*x[,5])^2+
    6*x[,7]^2+3*x[,8]^2- 2*x[,8]^3+
    10*dumx[,(quants*5)]+9*dumx[,(quants*5)+1]+3*dumx[,(quants*5)+2]+
    10*dumx[,(quants*8)+1]+5*dumx[,(quants*8)+2]+7*dumx[,(quants*8)+3]+
    1*dumx[,(quants*9)]+6*dumx[,(quants*9)+3]-6*dumx[,(quants*9)+5]+
    alpha*W%*%prey-.8*time+.000002*time^2+.0000002*time^3+gblockfe+
    rnorm(N)
  #y<-ginv(I-alpha*W)%*%prey
  
  #### Estimating true relationships
    
  #set points  
    at <- as.matrix(seq(-87, -86.9, length.out = 10))
    bt <- as.matrix(seq(41, 41.1, length.out = 10))
    setp<-matrix(ncol=2,nrow =110 )
    for(i in 1:dim(at)[1]){
      for(j in 1:dim(bt)[1]){
        setp[(i*dim(bt)[1])+j-1,1:2]<-cbind(at[i],bt[j])
      }
    }
    setp<-setp[10:109,]
    #dsetp<-distm( cbind(long, lat),setp,fun = distHaversine)
  ## Overspecify
  xds<-cbind(lat,long,lat,long,lat,long,
             lat,long,lat,long,lat,long,
             #lat+long,lat*long,lat^2*long^2,lat/long,
             #1/(lat*long),lat^2+long^2,
             blockmfe,x,W%*%y,time,time^2,time^3,polyx,dumx)
  feds<-cbind(model.matrix(~as.factor(floor(time))))
  xtl<-cbind(x,W%*%y,lat,long,time)
  xgrf<-cbind(xds,feds)
  
  #GRF
  c.forest<-causal_forest(X=xgrf, Y=y, W=Treatment,
                          num.trees=1000000,min.node.size=10,num.threads =100,
                          mtry = ceiling(ncol(xgrf)*.5) )
  
  xtest<- matrix(0, 10100, dim(xgrf)[2])
  a <- as.matrix(seq(-87, -86.9, length.out = 100))
  b <- as.matrix(seq(41, 41.1, length.out = 100))
 
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
  xtest1[,(dim(blockmfe)[2]+13)] <- as.matrix(seq(min(x[,4]), max(x[,4]), length.out = 100))
  
  xtest<-rbind(xtest,xtest1)
  #dsetptest<-distm( cbind(b, a),setp,fun = distHaversine)
  #xtest[,(dim(blockmfe)[2]+13):((dim(blockmfe)[2]+13)+dim(dsetptest)[2])]<-dsetptest
  
  tau.hat <- predict(c.forest, xtest, estimate.variance = TRUE)
  t<-tau.hat$predictions/sqrt(tau.hat$variance.estimates/N)
  p<-ifelse(abs(t)>1.96,1,0)
  
  #Create a function to generate a continuous color palette
  rbPal <- colorRampPalette(c('red','blue'))
  
  #This adds a column of color values
  # based on the y values
  Col <- rbPal(10)[as.numeric(cut(tau.hat$predictions[1:9999],breaks = 10))]
  
  plot(xtest[1:9999,1],xtest[1:9999,2],pch = 20,col = Col)
  legend("topleft",title="Decile",legend=c(1:10),col =rbPal(10),pch=20)
  
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
  
  
  ###################################################################################
  #DML
  
  xdml<-cbind(xds,feds)
  
  # = Cross-fitting DML = #
  # = Split sample = #
  I<-sort(sample(1:N,N/2))
  IC<-setdiff(1:N,I)
  # = compute ghat on both sample = #
  
  model11<-mcSuperLearner(Y=y[IC& Treatment==1],X=xdml[IC & Treatment==1,],
                         SL.library=SL.library$as.list(), family=gaussian(),
                         #obsWeights = train$as.matrix.weight.d.weightdlog,
                         method="method.NNLS", verbose=TRUE)
  model10<-mcSuperLearner(Y=y[IC& Treatment==0],X=xdml[IC & Treatment==0,],
                          SL.library=SL.library$as.list(), family=gaussian(),
                          #obsWeights = train$as.matrix.weight.d.weightdlog,
                          method="method.NNLS", verbose=TRUE)
  model21<-mcSuperLearner(Y=y[I& Treatment==1],X=xdml[I & Treatment==1,],
                         SL.library=SL.library$as.list(), family=gaussian(),
                         #obsWeights = train$as.matrix.weight.d.weightdlog,
                         method="method.NNLS", verbose=TRUE)
  model20<-mcSuperLearner(Y=y[I& Treatment==0],X=xdml[I & Treatment==0,],
                          SL.library=SL.library$as.list(), family=gaussian(),
                          #obsWeights = train$as.matrix.weight.d.weightdlog,
                          method="method.NNLS", verbose=TRUE)
  
  #model1<-randomForest(xdml[IC,],y[IC],maxnodes = 10)
  #model2<-randomForest(xdml[I,],y[I], maxnodes = 10)
  G11<-predict(model11,xdml[I& Treatment==1,], onlySL = T)
  G21<-predict(model21,xdml[IC& Treatment==1,], onlySL = T)
  G10<-predict(model10,xdml[I& Treatment==0,], onlySL = T)
  G20<-predict(model20,xdml[IC& Treatment==0,], onlySL = T)
  # = Compute mhat and vhat on both samples = #
  
  modeld1<-mcSuperLearner(Y=Treatment[IC],X=xdml[IC,],
                          SL.library=SL.library$as.list(), family=gaussian(),
                          #obsWeights = train$as.matrix.weight.d.weightdlog,
                          method="method.NNLS", verbose=TRUE)
  modeld2<-mcSuperLearner(Y=Treatment[I],X=xdml[I,],
                          SL.library=SL.library$as.list(), family=gaussian(),
                          #obsWeights = train$as.matrix.weight.d.weightdlog,
                          method="method.NNLS", verbose=TRUE)
  #modeld1<-randomForest(xdml[IC,],Treatment[IC],maxnodes = 10)
  #modeld2<-randomForest(xdml[I,],Treatment[I],maxnodes = 10)
  M1<-predict(modeld1,xdml[I,], onlySL = T)
  M2<-predict(modeld2,xdml[IC,], onlySL = T)
  #V1<-Treatment[I]-M1$pred
  #V2<-Treatment[IC]-M2$pred
  
  # = Compute Cross-Fitting DML theta
  theta1<-G11$pred-G10$pred + (Treatment[I]*(y[I]-G11$pred)/M1$pred)-((1-Treatment[I])*(y[I]-G10$pred)/(1-M1$pred))
  theta2<-G21$pred-G20$pred + (Treatment[IC]*(y[IC]-G21$pred)/M2$pred)-((1-Treatment[IC])*(y[IC]-G20$pred)/(1-M2$pred))
  theta_cf<-mean(c(theta1,theta2))
  #theta1<-mean(V1*(y[I]-G1$pred))/mean(V1*Treatment[I])
  #theta2<-mean(V2*(y[IC]-G2$pred))/mean(V2*Treatment[IC])
  #theta_cf<-mean(c(theta1,theta2))
  zeta1<-y[I]-Treatment[I]*theta_cf-G1$pred
  zeta2<-y[IC]-Treatment[IC]*theta_cf-G2$pred
  sigma1<-(mean(V1*V1)^(-2))*mean(V1^2*zeta1^2)
  sigma2<-(mean(V2*V2)^(-2))*mean(V2^2*zeta2^2)
  sigma_cf<-mean(c(sigma1,sigma2))
  t1<-(theta_cf-5)/sqrt(sigma_cf/N)
  t2<-(theta_cf-0)/sqrt(sigma_cf/N)
  betahat.dml[j,]<-theta_cf
  
  
  
  
  print(paste0(j, ' of ',n.sims))
}

# Print type 1 and type 2 errors
print(paste0('type 1 error naive OLS rate = ',colMeans(type1.tmle)))
print(paste0('type 1 error naive DS rate = ',colMeans(type1.ds)))
print(paste0('type 1 error overspecified OLS rate = ',colMeans(type1.ols.os)))
print(paste0('type 1 error overspecified DS rate = ',colMeans(type1.ds.os)))

print(paste0('type 2 error naive OLS rate = ',colMeans(type2.tmle)))
print(paste0('type 2 error naive DS rate = ',colMeans(type2.ds)))
print(paste0('type 2 error overspecified OLS rate = ',colMeans(type2.ols.os)))
print(paste0('type 2 error overspecified DS rate = ',colMeans(type2.ds.os)))


plot (density(betahat.ds),col="red",main = "",xlab="Treatment Coefficient Estimate (True = 1)")
title(main = "Density Plot of naive Double Selection")

plot (density(betahat.tmle),col="red",main = "",xlab="Treatment Coefficient Estimate (True = 1)")
title(main = "Density Plot of TMLE")

plot (density(betahat.ds.os),col="red",main = "",xlab="Treatment Coefficient Estimate (True = 1)")
title(main = "Density Plot of overspecified Double Selection")

plot (density(betahat.ols.os),col="red",main = "",xlab="Treatment Coefficient Estimate (True = 1)")
title(main = "Density Plot of overspecified OLS")




















