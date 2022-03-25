#### Gustavo Brusse - EDSD 2018/2019 ###
#### Assignment Event History 1 ########

## Creating the database and a survival object
ti<-c(rep(7,30), rep(7,211), rep(14,355), rep(21,369))
delta.i<-c(rep(1,30), rep(1,211), rep(1,355), rep(0,369))

## 1.1
## Weibull likelihood function
like.weibull<- function(par){
  a <- par[1]
  b <- par[2]
  logL <- 30*log(pweibull(7,shape=a,scale=b))+211*log(pweibull(14,shape=a,scale=b)-pweibull(7,shape=a,scale=b))+355*log(pweibull(21,shape=a,scale=b)-pweibull(14,shape=a,scale=b))+369*log(1-pweibull(21,shape=a,scale=b))
    return(logL)
}

## Optimization 
p.start <- c(0.1, 0.1)
weibull.MLE <- optim(p.start, like.weibull, control=list(fnscale=-1),hessian=T)
weibull.MLE

# check convergence
weibull.MLE$conv  

a <- weibull.MLE$par[1]
b <- weibull.MLE$par[2]
a
b
## a=3.04; b=21.25

## 1.2
## Hessian[1,1]=-69.58 Hessian[2,2]=-11.86
H.a<--69.58
H.b<--11.86

se.a<-sqrt(-1/H.a)
se.b<-sqrt(-1/H.b)

## Confidence interval for a and b; 
upper.CIa=3.04+(1.96*se.a)
lower.CIa=3.04-(1.96*se.a)
upper.CIb=21.25+(1.96*se.b)
lower.CIb=21.25-(1.96*se.b)

## 1.3
## Reading the flies complete dataset
setwd("C:\\Users\\gbrusse\\Desktop\\Event History")
flies <- read.table("flies.txt",header=T)

## proportion full dataset
age2 <- flies$days[flies$days > 2] - 2   ## age2 is age in days after day 2
alive.30<- sum(age2>=30)/length(age2)

## proportion censoring data
p.30<-1-(pweibull(30,a,b))

## Delta-method for estimating S.E
Var<-solve(-weibull.MLE$hessian)
g.a<-(-(1-pweibull(30,a,b))*((30/b)^a)*log(30/b))
g.b<-((1-pweibull(30,a,b))*a*((30/b)^a))/b
se30<- sqrt(t(c(g.a, g.b)) %*% Var %*% c(g.a, g.b))

# Confidence Interval for S(30)
upper.CI30=p.30+(1.96*se30)
lower.CI30=p.30-(1.96*se30)

## Question 2
## reading the data ##
# Reading input dataset
setwd("C:\\Users\\gbrusse\\Desktop\\Event History")
data.meno <- read.table("C:\\Users\\gbrusse\\Desktop\\Event History\\meno.txt", header=TRUE, sep=",")

# 2.1
# Hazard and Survival function for a = -20 and b = 0.5
a=-20
b=0.5
t=data.meno$age.interview

ht<-(b*exp(a+b*t))/(1+exp(a+b*t))
Ht<-log((1+exp(a+b*t))/(1+exp(a)))
St<-exp(-Ht)
Ft<-1-(((1+exp(a+b*t))/(1+exp(a)))^(-1))

# 2.3
# Log-Likelihood
log.like<-function(par,data){
  a=par[1]
  b=par[2]
  t=data[,1]
  delta=data[,2]
  Ht<-log((1+exp(a+b*t))/(1+exp(a)))  
  St<-exp(-Ht) 
  Ft<-1-(((1+exp(a+b*t))/(1+exp(a)))^(-1))
  log.like<-sum(delta*log(Ft)+(1-delta)*(-Ht)) 
  return(log.like)
}

# 2.4
p.start <- c(0.1,0.1)
MLE<-optim(p.start,log.like,data=data.meno,control=list(fnscale=-1),hessian=T)

# Parameters
a=MLE$par[1]
b=MLE$par[2]

# CI 
upper.a<--21.47+(1.96*1.15)
lower.a<--21.47-(1.96*1.15)
upper.b<-0.43+(1.96*0.02)
lower.b<-0.43-(1.96*0.02)

# 2.5
# Plots
graphs<-function(t,par){
  a=par[1];b=par[2]
  h<-(b*exp(a+b*t))/(1+exp(a+b*t)) 
  H<-log((1+exp(a+b*t))/(1+exp(a))) 
  S<-exp(-H) 
  f<-h*S  
  
  plot(t,h,main="hazard",xlab="age",ylab="")
  plot(t,S,main="survival",xlab="age",ylab="")
  plot(t,f,main="density",xlab="age",ylab="")
  
}

p<-graphs(data.meno$age.interview,MLE$par)

# 2.6
# Percentage of women who has not experienced menopause by age 55
age <- 55
sum(data$menop[data.meno$age.interview>55])/length(data.meno$age.interview)*100 
# Maximum Likelihood Estimate
MLE<-exp(-log((1+exp(a+b*age))/(1+exp(a))))
MLE2 <- Est*100

# CI for estimated parameters
# Delta method
g.a<-exp(a)*(exp(b*age)-1)/((exp(a+b*age)+1)^2)
g.b<-((exp(a)+1)*age*exp(a+b*age))/((exp(a+b*age)+1)^2)
# SE
se.Est <- sqrt(t(c(g.a, g.b)) %*% V %*% c(g.a, g.b))
# CI
lb.est<- Est-qnorm(.975)*se.Est
ub.est<- Est+qnorm(.975)*se.Est
lb.est
ub.est
# Estimates, SE and CI 
Finalest<-cbind(c(Est*100,se.Est,lb.est*100,ub.est*100))
row.names(Finalest)<-c("estimate","SE","Lower","Upper")
colnames(Finalest)<-"Percentage"
t(Finalest)
