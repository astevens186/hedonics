
## Preliminaries
rm(list=ls())
gc()

path<- "P:/HousingProject/"

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

set.seed(12345)
#set n=200 and number of replications=200
N <- 2000
P<-15
n.sims <- 20


# Create data  
x <- matrix(runif(N*P), nrow = N, ncol = P)
nvar<-P
quants<-5
dumx<-matrix(nrow=N,ncol=quants*nvar)
for(i in 1:nvar){
  indx<- factor(as.numeric(cut2(x[,i], g=quants)))
  dumx[,((quants*i)-(quants-1)):(quants*i)]<-model.matrix(~indx-1)
}



I<-diag(N)
alpha<-.75
lat<-.08*runif(N)-87
long<-1*runif(N)+41
time<-1990+(10*runif(N))
psn<-rnorm(N)
Treatment<- ifelse(psn-median(psn)>0,1,0) 
#Treatment<- ifelse(time-median(time)>0,1,0) 

truegridquants<-15
indlat<- factor(as.numeric(cut2(lat, g=truegridquants)))
indlong<- factor(as.numeric(cut2(long, g=truegridquants)))
blockfe<-model.matrix(~indlat:indlong-1)

modelgridquants<-20
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
TE<-ifelse(TEi<25,0,TEi)
summary(TE)
#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))

#This adds a column of color values
# based on the y values
Col <- rbPal(10)[as.numeric(cut(TE,breaks = 10))]

plot(lat,long,pch = 20,col = Col)
legend("topleft",title="Decile",legend=levels(cut(TE,breaks = 10)),col =rbPal(10),pch=20)

TEi<-pmin(100000*(1/(disttreat+2000)),35*rep(1,N))
TE2<-ifelse(TEi<25,0,TEi)
summary(TE2)
#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))

#This adds a column of color values
# based on the y values
Col <- rbPal(10)[as.numeric(cut(TE2,breaks = 10))]
jpeg(paste0(path,'trueeffect1.jpg'))

plot(long,lat,pch = 20,col = Col)
legend("topleft",title="Decile",legend=levels(cut(TE2,breaks = 10)),col =rbPal(10),pch=20)
dev.off()

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
if(FALSE){
  maty<-matrix(nrow = N, ncol = dim(blockfe)[2])
  for(i in 1:dim(blockfe)[2]){
    
    py<-x[,2] * 2.5*runif(1) + 60*x[,3]*runif(1)-2*x[,3]^2*runif(1)-2*x[,3]^3*runif(1)+
      #5*(x[,4]*x[,5])+2*(x[,4]*x[,5])^2+
      6*x[,7]^2*runif(1)+30*x[,8]^2*runif(1)- 20*x[,8]^3*runif(1)+
      10*dumx[,(quants*5)]*runif(1)+9*dumx[,(quants*5)+1]*runif(1)+3*dumx[,(quants*5)+2]*runif(1)+
      10*dumx[,(quants*8)+1]*runif(1)+5*dumx[,(quants*8)+2]*runif(1)+7*dumx[,(quants*8)+3]*runif(1)+
      1*dumx[,(quants*9)]*runif(1)+6*dumx[,(quants*9)+3]*runif(1)-6*dumx[,(quants*9)+5]*runif(1)+
      -.01*time*runif(1)+.0000001*runif(1)*time^2+.00000001*runif(1)*time^3+(90000/(distns+2000))*runif(1)+
      (90/(distns+2000))*time*runif(1)+
      rnorm(N)
    maty[,i]<-py*blockfe[,i]
    
  }
  y<-rowSums(maty)
  
}

py<- x[,2] * 2.5 + 60*x[,3]-2*x[,3]^2-2*x[,3]^3+
  #5*(x[,4]*x[,5])+2*(x[,4]*x[,5])^2+
  6*x[,7]^2+30*x[,8]^2- 20*x[,8]^3+
  10*dumx[,(quants*5)]+9*dumx[,(quants*5)+1]+3*dumx[,(quants*5)+2]+
  10*dumx[,(quants*8)+1]+5*dumx[,(quants*8)+2]+7*dumx[,(quants*8)+3]+
  1*dumx[,(quants*9)]+6*dumx[,(quants*9)+3]-6*dumx[,(quants*9)+5]+
  -.01*time+.0000001*time^2+.00000001*time^3+(90000/(distns+2000))+
  (90/(distns+2000))*time+
  rnorm(N)



y <- TE+ 5*Treatment*rnorm(N)+#2*Treatment*(time-1989)+ 
  py #+ 10*rnorm(N)
#y<-ginv(I-alpha*W)%*%prey
rm(blockfe,dumx,dW,I,T1,T2,Tdp,W,Wf,Ws,Wtp)
gc()
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
sdf<-40
splat<-bs(lat, df = sdf)
splong<-bs(long, df = sdf)
spint<-model.matrix(~splat:splong)

xds<-cbind(lat,long,lat,long,lat,long,
           lat,long,lat,long,lat,long,
           #lat+long,lat*long,lat^2*long^2,lat/long,
           #1/(lat*long),lat^2+long^2,
           blockmfe,time,bs(time,df =sdf),x)
feds<-cbind(model.matrix(~as.factor(floor(time))))
#xtl<-cbind(x,W%*%y,lat,long,time)
xgrf<-cbind(xds,feds,bs(disttreat,df =sdf),spint,splat,splong,disttreat)

#xds2<-cbind(lat,long, blockmfe,time,x,bs(time, df = sdf),polyx,dumx )
#feds<-cbind(model.matrix(~as.factor(floor(time))))
# xgrf<-cbind(xds2,feds,spint,splat,splong,disttreat*Treatment,bs(disttreat,df =sdf)*Treatment,disttreat,bs(disttreat,df =sdf))

###########################################################################################
#semi-parametric causal forest
zspgrf<-cbind(lat,long,bs(disttreat,df =sdf),spint,splat,splong,disttreat)

pcolf<-.5
trees<-1000
c.forest.ey <-causal_forest(X=zspgrf, Y=y, W=Treatment,
                            num.trees=trees,#min.node.size=10,num.threads =100,
                            mtry = ceiling(ncol(xgrf)*pcolf) )
beta.mat<-matrix(nrow = dim(x)[2],ncol = 1)
for(i in 1:dim(x)[2]){
  yx<-x[,i]
  c.forest.ex <-causal_forest(X=zspgrf, Y=yx, W=Treatment,
                              num.trees=trees,#min.node.size=10,num.threads =100,
                              mtry = ceiling(ncol(xgrf)*pcolf/2) )
  ey<-y-c.forest.ey$Y.hat
  ex<-y-c.forest.ex$Y.hat
  
  lm.beta<-lm(ey~ex)
  beta.mat[i,]<-lm.beta$coefficients[2]
}

e.para<-y-(x%*% beta.mat)

c.forest.sp <-causal_forest(X=zspgrf, Y=e.para, W=Treatment,
                            num.trees=trees,#min.node.size=10,num.threads =100,
                            mtry = ceiling(ncol(xgrf)*pcolf) )


ln<-100
nt<-ln*ln

xtest<- matrix(0, nt+ln, dim(zspgrf)[2])
a <- as.matrix(seq(min(long), max(long), length.out = ln))
b <- as.matrix(seq(min(lat), max(lat), length.out = ln))

for(i in 1:dim(a)[1]){
  for(j in 1:dim(b)[1]){
    xtest[(i*dim(b)[1])+j,1:2]<-cbind(a[i],b[j])
  }
}
xtest<-xtest[(ln+1):dim(xtest)[1],]

xtest[,dim(zspgrf)[2]]<-distm(cbind(xtest[,1], xtest[,2]),c(treatlong,treatlat), fun = distHaversine)

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
tau.hat <- predict(c.forest.sp, xtest, estimate.variance = TRUE)
upper<-tau.hat$predictions+1.96*sqrt(tau.hat$variance.estimates)
lower<-tau.hat$predictions-1.96*sqrt(tau.hat$variance.estimates)
con<-ifelse(sign(upper)>0 &sign(lower)>0,1,0)
con<-ifelse(sign(upper)<0 &sign(lower)<0,-1,con)
con<-ifelse(sign(upper)>0 &sign(lower)<0,0,con)

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))

#This adds a column of color values
# based on the y values
Col <- rbPal(10)[as.numeric(cut(tau.hat$predictions[1:nt],breaks = 10))]

jpeg(paste0(path,'cfpredeffectsemipar1.jpg'))

plot(xtest[1:nt,2],xtest[1:nt,1],pch = 20,col = Col)
legend("topleft",title="Decile",legend=levels(cut(tau.hat$predictions[1:nt],breaks = 10)),col =rbPal(10),pch=20)
dev.off()

rbPal3 <- colorRampPalette(c('blue','red'))

#This adds a column of color values
# based on the y values
Col3<-ifelse(con>0, "#FF0000","#7F007F")
COl3<-ifelse(con<0,"#0000FF",Col3)

jpeg(paste0(path,'cfsigneffectsemipar1.jpg'))

plot(xtest[1:nt,2],xtest[1:nt,1],pch = 20,col = Col3)
legend("topleft",title="Significance",legend=c("negative","insignificant","positive"),col =rbPal3(3),pch=20)

dev.off()
rbPal3 <- colorRampPalette(c('blue','red'))

#This adds a column of color values
# based on the y values

rm(Col3,con)
gc()

xtest[,dim(zspgrf)[2]]<-distm(cbind(xtest[,1], xtest[,2]),c(treatlong,treatlat), fun = distHaversine)
tTEi<-pmin(100000*(1/(xtest[1:nt,dim(zspgrf)[2]])),35*rep(1,N))
TTE<-ifelse(tTEi<25,0,tTEi)

#TTE<-100000*(1/(xtest[1:9999,dim(xgrf)[2]]+2000))

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))

#This adds a column of color values
# based on the y values
Col <- rbPal(10)[as.numeric(cut(tau.hat$predictions[1:nt]-TTE,breaks = 10))]

jpeg(paste0(path,'cfbiassemipar1.jpg'))

plot(xtest[1:nt,2],xtest[1:nt,1],pch = 20,col = Col)
legend("topleft",title="Decile",legend=levels(cut(tau.hat$predictions[1:nt]-TTE,breaks = 10)),col =rbPal(10),pch=20)
dev.off()

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))

#This adds a column of color values
# based on the y values
Col <- rbPal(10)[as.numeric(cut(TTE,breaks = 10))]

jpeg(paste0(path,'trueeffecttsemipar1.jpg'))

plot(xtest[1:nt,2],xtest[1:nt,1],pch = 20,col = Col)
legend("topleft",title="Decile",legend=levels(cut(TTE,breaks = 10)),col =rbPal(10),pch=20)

dev.off()



rm(xtest)
gc()

############################################################################################

#GRF
c.forest<-causal_forest(X=xgrf, Y=y, W=Treatment,
                        num.trees=1000,#min.node.size=10,num.threads =100,
                        mtry = ceiling(ncol(xgrf)*.3) )
ln<-100
nt<-ln*ln

xtest<- matrix(0, nt+ln, dim(xgrf)[2])
a <- as.matrix(seq(min(long), max(long), length.out = ln))
b <- as.matrix(seq(min(lat), max(lat), length.out = ln))

for(i in 1:dim(a)[1]){
  for(j in 1:dim(b)[1]){
    xtest[(i*dim(b)[1])+j,1:2]<-cbind(a[i],b[j])
  }
}
xtest<-xtest[(ln+1):dim(xtest)[1],]

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

xtest1<- matrix(0, ln, dim(xgrf)[2])
xtest1[,(dim(blockmfe)[2]+13)] <- as.matrix(seq(min(time), max(time), length.out = ln))
xtest1[,(dim(xds)[2]+1):(dim(xds)[2]+dim(feds)[2])]<-model.matrix(~as.factor(floor(xtest1[,(dim(blockmfe)[2]+13)] )))
xtest1[,(dim(blockmfe)[2]+14):(dim(blockmfe)[2]+13+dim(bs(time,df =sdf))[2])]<-bs(xtest1[,(dim(blockmfe)[2]+13)],df =sdf)
xtest<-rbind(xtest,xtest1)
rm(xtest1)
gc()

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
upper<-tau.hat$predictions+1.96*sqrt(tau.hat$variance.estimates)
lower<-tau.hat$predictions-1.96*sqrt(tau.hat$variance.estimates)
con<-ifelse(sign(upper)>0 &sign(lower)>0,1,0)
con<-ifelse(sign(upper)<0 &sign(lower)<0,-1,con)
con<-ifelse(sign(upper)>0 &sign(lower)<0,0,con)

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))

#This adds a column of color values
# based on the y values
Col <- rbPal(10)[as.numeric(cut(tau.hat$predictions[1:nt],breaks = 10))]

jpeg(paste0(path,'cfpredeffect1.jpg'))

plot(xtest[1:nt,2],xtest[1:nt,1],pch = 20,col = Col)
legend("topleft",title="Decile",legend=levels(cut(tau.hat$predictions[1:nt],breaks = 10)),col =rbPal(10),pch=20)
dev.off()

rbPal3 <- colorRampPalette(c('blue','red'))

#This adds a column of color values
# based on the y values
Col3<-ifelse(con>0, "#FF0000","#7F007F")
COl3<-ifelse(con<0,"#0000FF",Col3)

jpeg(paste0(path,'cfsigneffect1.jpg'))

plot(xtest[1:nt,2],xtest[1:nt,1],pch = 20,col = Col3)
legend("topleft",title="Significance",legend=c("negative","insignificant","positive"),col =rbPal3(3),pch=20)

dev.off()
rbPal3 <- colorRampPalette(c('blue','red'))

#This adds a column of color values
# based on the y values
jpeg(paste0(path,'cftimeeffect1.jpg'))

plot(xtest[(nt+1):(nt+ln),(dim(blockmfe)[2]+13)],tau.hat$predictions[(nt+1):(nt+ln),],pch = 20,col = Col3)
legend("top",title="Significance",legend=c("t< -1.96","-1.96<t<1.96","t>1.96"),col =rbPal3(3),pch=20)
rm(Col3,con)
gc()
dev.off()
xtest[,dim(xgrf)[2]]<-distm(cbind(xtest[,1], xtest[,2]),c(treatlong,treatlat), fun = distHaversine)
tTEi<-pmin(100000*(1/(xtest[1:nt,dim(xgrf)[2]]+2000)),35*rep(1,N))
TTE<-ifelse(tTEi<25,0,tTEi)

#TTE<-100000*(1/(xtest[1:9999,dim(xgrf)[2]]+2000))

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))

#This adds a column of color values
# based on the y values
Col <- rbPal(10)[as.numeric(cut(tau.hat$predictions[1:nt]-TTE,breaks = 10))]

jpeg(paste0(path,'cfbias1.jpg'))

plot(xtest[1:nt,2],xtest[1:nt,1],pch = 20,col = Col)
legend("topleft",title="Decile",legend=levels(cut(tau.hat$predictions[1:nt]-TTE,breaks = 10)),col =rbPal(10),pch=20)
dev.off()

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))

#This adds a column of color values
# based on the y values
Col <- rbPal(10)[as.numeric(cut(TTE,breaks = 10))]

jpeg(paste0(path,'trueeffectt1.jpg'))

plot(xtest[1:nt,2],xtest[1:nt,1],pch = 20,col = Col)
legend("topleft",title="Decile",legend=levels(cut(TTE,breaks = 10)),col =rbPal(10),pch=20)

dev.off()



rm(xtest)
gc()
#######################################################################
#FE
feds<-model.matrix(~as.factor(floor(time)-1))
xfe<-model.matrix(~Treatment*blockmfe+x+feds-1)
xfe<-cbind(Treatment,blockmfe,model.matrix(~Treatment:blockmfe-1),x,feds)#,disttreat,Treatment*disttreat)
fit<-lm(y~xfe)
summary(fit)
ln<-1000
nt<-ln*ln

xtest3<- matrix(0, nt, dim(xfe)[2])
a <- as.matrix(seq(min(long), max(long), length.out = ln))
b <- as.matrix(seq(min(lat), max(lat), length.out = ln))
xtest2<-matrix(0, nt+ln, dim(xfe)[2])
for(i in 1:dim(a)[1]){
  for(j in 1:dim(b)[1]){
    xtest2[(i*dim(b)[1])+j,1:2]<-cbind(a[i],b[j])
  }
}
xtest2<-xtest2[(ln+1):dim(xtest2)[1],]


indmlat<- factor(as.numeric(cut2(xtest2[,2], g=modelgridquants)))
indmlong<- factor(as.numeric(cut2(xtest2[,1], g=modelgridquants)))
xtest3[,2:(dim(blockmfe)[2]+1)]<-model.matrix(~indmlat:indmlong-1)
xtest3[,(dim(blockmfe)[2]+2):(2*dim(blockmfe)[2]+1)]<-model.matrix(~indmlat:indmlong-1)
xtest3[,(2*dim(blockmfe)[2]+1)]<-rep(0,nt)

xtest3[,1]<-rep(1,nt)

#xtest3[,dim(xfe)[2]-1]<-distm(cbind(xtest3[,1], xtest3[,2]),c(treatlong,treatlat), fun = distHaversine)
# xtest3[,dim(xfe)[2]]<-xtest3[,dim(xfe)[2]-1]
if(FALSE){
  for(i in 1:dim(xtest3)[2]){
    if(sum(xtest3[,i])==0){
      xtest3[,i]<-median(xfe[,i])
    }
  }
}

xtest3<-cbind(rep(1,nt),xtest3)
fit$coefficients[is.na(as.numeric(coef(fit)))]<-0
lm.pred1<-xtest3 %*% as.matrix(fit$coefficients )
#lm.pred1<-predict(object=fit,newdata=data.frame(xtest3))
#xtest3[,1]<-rep(0,nt)
#xtest3[,dim(xfe)[2]]<-rep(0,nt)
#xtest3[,(dim(blockmfe)[2]+2):(2*dim(blockmfe)[2]+1)]<-rep(0,nt)
#lm.pred0<-predict(object=fit,newdata=data.frame(xtest3))
#lm.pred0<-xtest3 %*% as.matrix(fit$coefficients )
#lm.te<-lm.pred1-lm.pred0
rm(xtest3)
gc()
xtest2[,dim(xfe)[2]]<-distm(cbind(xtest2[,1], xtest2[,2]),c(treatlong,treatlat), fun = distHaversine)

tTEi<-pmin(100000*(1/(xtest2[,dim(xfe)[2]]+2000)),35*rep(1,N))
TTE2<-ifelse(tTEi<25,0,tTEi)

#TTE<-100000*(1/(xtest2[,dim(xfe)[2]]+2000))

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))

#This adds a column of color values
# based on the y values
Col <- rbPal(10)[as.numeric(cut(lm.pred1,breaks = 10))]


jpeg(paste0(path,'lmpredeffect1.jpg'))

plot(xtest2[,2],xtest2[,1],pch = 20,col = Col)
legend("topleft",title="Decile",legend=levels(cut(lm.pred1,breaks = 10)),col =rbPal(10),pch=20)

dev.off()
#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('red','blue'))

#This adds a column of color values
# based on the y values
Col <- rbPal(10)[as.numeric(cut(lm.pred1-TTE2,breaks = 10))]

jpeg(paste0(path,'lmbias1.jpg'))

plot(xtest2[,2],xtest2[,1],pch = 20,col = Col)
legend("topleft",title="Decile",legend=levels(cut(lm.pred1-TTE2,breaks = 10)),col =rbPal(10),pch=20)

dev.off()
summary(lm.pred1-TTE2)
summary(tau.hat$predictions[1:nt]-TTE)

