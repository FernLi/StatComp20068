## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----warning=FALSE------------------------------------------------------------
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
lm.D9 <- lm(weight ~ group)
plot(lm.D9)

## -----------------------------------------------------------------------------
knitr::kable (head(iris))

## -----------------------------------------------------------------------------
set.seed(20)
n=100
u=runif(n)
x=2/sqrt(1-u)  #F(x)=1-(2/x)^2 ,x>=b>0,a>0
hist(x,prob=TRUE,main = expression(f(x)==8/x^3))
y=seq(2,max(x),0.01)
lines(y,8*(1/(y^3)))

## -----------------------------------------------------------------------------
n=10000
u1=2*runif(n)-1;u2=2*runif(n)-1;u3=2*runif(n)-1;
x=vector();
for (i in 1:n){
  if (abs(u3[i])>abs(u2[i]) && abs(u3[i])>abs(u1[i])){
    x[i]=u2[i]
    }else
    x[i]=u3[i]
}
hist(x,prob=TRUE,main = expression(f(x)==3/4*(1-x^2)))
y=seq(-1,1,0.01)
lines(y,3/4*(1-y^2))

## -----------------------------------------------------------------------------
F=-0.25*y^3+0.75*y+0.5
ks.test(x,F)

## -----------------------------------------------------------------------------
set.seed(70)
n=1000
u=runif(n)
x=2/((1-u)^(1/4))-2  #F(x)=1-(2/x)^2 ,x>=b>0,a>0
hist(x,prob=TRUE,main = expression(f(x)==64*(2+x)^{-5}))
y=seq(0,max(x),0.01)
lines(y,64*(2+y)^{-5})

## -----------------------------------------------------------------------------
set.seed(66)
number=1000000
exact_value=1/2
u=runif(number,min = 0,max = pi/3)
MC_value=pi/3*mean(sin(u))
print(c(MC_value,exact_value,abs(MC_value-exact_value)))

## -----------------------------------------------------------------------------
set.seed(66)
number=1000
exact_value=exp(1)-1
n=1000
MC_value1=vector()
MC_value2=vector()
for (i in 1:n){
  u1=runif(number)
  u2=runif(number)
  MC_value1[i]=1/2*mean(exp(u1)+exp(1-u1))
  MC_value2[i]=mean((exp(u1)+exp(u2))/2)
}
sd1=sd(MC_value1)
sd2=sd(MC_value2)
print(c(exact_value,mean(MC_value1),mean(MC_value2),sd1,sd2))

## -----------------------------------------------------------------------------
print(1-var(MC_value1)/var(MC_value2))

## -----------------------------------------------------------------------------
set.seed(0)
m=10000
theta.hat=var.hat=numeric(2)
g=function(x){
  (x^2*exp(-x^2/2)/sqrt(2*pi))*(x>1)
}

u=runif(m)                    ##inverse transform method 
x=sqrt(1-2*log(1-u))    
fg1=g(x)/(exp(1/2-x^2/2)*x)
theta.hat[1]=mean(fg1)
var.hat[1]=var(fg1)

x=rexp(m,1)
fg2=g(x)/(exp(-x))
theta.hat[2]=mean(fg2)
var.hat[2]=var(fg2)

print(rbind(theta.hat,var.hat))

## -----------------------------------------------------------------------------
set.seed(0)
m=100;k=5;r=m/k;n=50;t2=numeric(k)
estimates=matrix(0,n,2)

g=function(x){exp(-x)/(1+x^2)*(x>0)*(x<1)}
f=function(x){exp(-x)/(exp((1-j)/5)-exp(-j/5))}

for (i in 1:n){
  estimates[i,1]=mean(g(runif(m)))
  for ( j in 1:k){
    u=runif(m/k)
    x=-log(exp(-(j-1)/5)-u*(exp((1-j)/5)-exp(-j/5)))
    t2[j]=mean(g(x)/f(x))
  }
  estimates[i,2]=sum(t2)
}
apply(estimates,2,mean)
apply(estimates,2,var)

## -----------------------------------------------------------------------------
set.seed(0)
n=20
alpha=0.05
u=0;sigma=1;
number=100000
Icl=vector();ucl=vector();

for (i in 1:number){
  x=rnorm(n,mean = u,sd=sigma)
  s=sqrt(var(x)*n/(n-1))
  ucl[i]=mean(x)+s/sqrt(n)*qt(1-alpha/2,n-1)
  Icl[i]=mean(x)-s/sqrt(n)*qt(1-alpha/2,n-1)
}
Icl1=mean(Icl);ucl1=mean(ucl)

a=0;
for (i in 1:number){
  if (Icl[i]<0&&ucl[i]>0)a=a+1
  else a=a
}
p=a/number

print(c(Icl1,ucl1))
print(p)

## -----------------------------------------------------------------------------
set.seed(0)
n=20;alpha=0.05;number=1000;
Icl=vector();ucl=vector();

for (i in 1:number){
  x=rchisq(n,df=2)
  s=sqrt(var(x)*n/(n-1))
  ucl[i]=mean(x)+s/sqrt(n)*qt(1-alpha/2,n-1)
  Icl[i]=mean(x)-s/sqrt(n)*qt(1-alpha/2,n-1)
}
Icl1=mean(Icl);ucl1=mean(ucl)

a=0;
for (i in 1:number){
  if (Icl[i]<2&&ucl[i]>2)a=a+1
  else a=a
}
p=a/number

print(c(Icl1,ucl1))
print(p)

## -----------------------------------------------------------------------------
set.seed(12345)

sk = function(x) {
  xbar = mean(x)
  m3 = mean((x - xbar)^3)
  m2 = mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}

# beta(a,a)
pwr_beta = function(a){
 alpha = 0.1
 n = 20
 m = 1e4
 N = length(a)
 pwr = numeric(N)
 cv = qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
 
 for (j in 1:N) { 
  sktests = numeric(m)
  for (i in 1:m) { 
   x = rbeta(n, a[j], a[j])
   sktests[i] = as.integer(abs(sk(x))>= cv)
  }
  pwr[j] = mean(sktests)
 }
 se = sqrt(pwr * (1-pwr) / m) 
 return(list(pwr = pwr,se = se))
}

 a = c(seq(0,1,0.1),seq(1,20,1),seq(20,100,10))
 pwr = pwr_beta(a)$pwr
 # plot the power
 se = pwr_beta(a)$se
 plot(a, pwr, type = "b", xlab = "a", ylab = "pwr", pch=16)
 abline(h = 0.1, lty = 2)
 lines(a, pwr+se, lty = 4)
 lines(a, pwr-se, lty = 4)

## -----------------------------------------------------------------------------
set.seed(0)
m=1000;alpha=0.055;
sigma1=1;sigma2=1.5
number1=numeric(m);number2=numeric(m)
pwr1=numeric(3);pwr2=numeric(3)

count5test=function(x,y){
  x=x-mean(x)
  y=y-mean(y)
  outx=sum(x>max(y))+sum(x<min(y))
  outy=sum(y>max(x))+sum(y<min(x))
  return(as.integer(max(c(outx,outy))>=5))
}

a=c(10,100,2000)
for (i in 1:3){
  n1=a[i];n2=a[i];
  for (j in 1:m){
    x=rnorm(n1,0,sqrt(sigma1))
    y=rnorm(n2,0,sqrt(sigma2))
    s1=var(x)
    s2=var(y)
    F_value=s1/s2
    number1[j]=count5test(x,y)
    number2[j]=as.integer(F_value>qf(1-alpha/2,n1-1,n2-1)||F_value<qf(alpha/2,n1-1,n2-1))
  }
  pwr1[i]=mean(number1)
  pwr2[i]=mean(number2)
}
print(cbind(pwr1,pwr2))

## ----warning=FALSE------------------------------------------------------------
library(MASS)
Mardia<-function(mydata){
  n=nrow(mydata)
  c=ncol(mydata)
  central<-mydata
  for(i in 1:c){
    central[,i]<-mydata[,i]-mean(mydata[,i])
  }
  sigmah<-t(central)%*%central/n
  a<-central%*%solve(sigmah)%*%t(central)
  b<-sum(colSums(a^{3}))/(n*n)
  test<-n*b/6
  chi<-qchisq(0.95,c*(c+1)*(c+2)/6)
  as.integer(test>chi)
}

set.seed(1234)
mu <- c(0,0,0)
sigma <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
m=1000
n<-c(10, 20, 30, 50, 100, 500)
#m: number of replicates; n: sample size
a=numeric(length(n))
for(i in 1:length(n)){
  a[i]=mean(replicate(m, expr={
    mydata <- mvrnorm(n[i],mu,sigma) 
    Mardia(mydata)
  }))
}
print(a)

## -----------------------------------------------------------------------------
library(MASS)
set.seed(7912)
set.seed(7912)
mu1 <- mu2 <- c(0,0,0)
sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
sigma2 <- matrix(c(100,0,0,0,100,0,0,0,100),nrow=3,ncol=3)
sigma=list(sigma1,sigma2)
m=1000
n=50
#m: number of replicates; n: sample size
epsilon <- c(seq(0, .06, .01), seq(.1, 1, .05))
N <- length(epsilon)
pwr <- numeric(N)
for (j in 1:N) { #for each epsilon
  e <- epsilon[j]
  sktests <- numeric(m)
  for (i in 1:m) { #for each replicate
    index=sample(c(1, 2), replace = TRUE, size = n, prob = c(1-e, e))
    mydata<-matrix(0,nrow=n,ncol=3)
    for(t in 1:n){
      if(index[t]==1) mydata[t,]=mvrnorm(1,mu1,sigma1) 
      else mydata[t,]=mvrnorm(1,mu2,sigma2)
    }
    sktests[i] <- Mardia(mydata)
  }
  pwr[j] <- mean(sktests)
}
plot(epsilon, pwr, type = "b",
     xlab = bquote(epsilon), ylim = c(0,1))
abline(h = .05, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(epsilon, pwr+se, lty = 3)
lines(epsilon, pwr-se, lty = 3)

## ----warning=FALSE------------------------------------------------------------
data(law,package = "bootstrap")
n=nrow(law)
x=law$LSAT
y=law$GPA

theta.hat=cor(x,y)
theta.jack=numeric(n)

for(i in 1:n){
  theta.jack[i]=cor(x[-i],y[-i])
}
bias.jack=(n-1)*(mean(theta.jack)-theta.hat)
se.jack=(n-1)*sqrt(var(theta.jack)/n)

round(c(bias.jack=bias.jack,se.jack=se.jack),9)

## ----warning=FALSE------------------------------------------------------------
library(boot)
data(aircondit,package = "boot")
B=600
n=nrow(aircondit)
#Bootstrap
func=function(dat,index){
  x=dat[index,]
  theta=mean(x)
  return(theta)
}
bootstrap_result=boot(data=aircondit,statistic=func,R=B)
boot.ci(bootstrap_result,conf=0.95,type=c("norm","basic","perc","bca"))

## ----warning=FALSE------------------------------------------------------------
##bootstrap
set.seed(0)
library("bootstrap")
lambda_hat=eigen(cov(scor))$values
theta_hat=lambda_hat[1]/sum(lambda_hat)
B=200
n=nrow(scor)

##Jackknif
theta_j=vector()
for (i in 1:n){
  x2=scor[-i,]
  lambda=eigen(cov(x2))$values
  theta_j[i]=lambda[1]/sum(lambda)
}
bias_j=(n-1)*(mean(theta_j)-theta_hat)
se_j=sqrt((n-1)*mean((theta_j-mean(theta_j))^2))
round(c(bias_j=bias_j,se_j=se_j),9)

## ----warning=FALSE------------------------------------------------------------
library(DAAG);attach(ironslag)
data("ironslag",package = "DAAG")
n = length(magnetic) #in DAAG ironslag

e11=e22=e33=e44=matrix(0,n*(n-1)/2,2)
# for n-fold cross validation
# fit models on leave-two-out samples
j=0;

for (k in 1:n) {
  for(i in k+1:n){
    if(i>n)break
    j=j+1
    y=magnetic[-k][-i+1]
    x=chemical[-k][-i+1]
    
    J1=lm(y~x)
    yhat11=J1$coef[1]+J1$coef[2]*chemical[k]
    yhat12=J1$coef[1]+J1$coef[2]*chemical[i]
    e11[j,1]=magnetic[k]-yhat11
    e11[j,2]=magnetic[i]-yhat12
    
    J2=lm(y~x+I(x^2))
    yhat21=J2$coef[1]+J2$coef[2]*chemical[k]+J2$coef[3]*chemical[k]^2
    yhat22=J2$coef[1]+J2$coef[2]*chemical[i]+J2$coef[3]*chemical[i]^2
    e22[j,1]=magnetic[k]-yhat21
    e22[j,2]=magnetic[i]-yhat22
    
    J3=lm(log(y)~x)
    logyhat31=J3$coef[1]+J3$coef[2]*(chemical[k])
    logyhat32=J3$coef[1]+J3$coef[2]*(chemical[i])
    yhat31=exp(logyhat31)
    yhat32=exp(logyhat32)
    e33[j,1]=magnetic[k]-yhat31
    e33[j,2]=magnetic[i]-yhat32
    
    J4=lm(log(y)~log(x))
    logyhat41=J4$coef[1]+J4$coef[2]*log(chemical[k])
    logyhat42=J4$coef[1]+J4$coef[2]*log(chemical[i])
    yhat41=exp(logyhat41)
    yhat42=exp(logyhat42)
    e44[j,1]=magnetic[k]-yhat41
    e44[j,2]=magnetic[i]-yhat42
  }
}
c(mean(e11^2),mean(e22^2),mean(e33^2),mean(e44^2))
lm(magnetic~chemical+I(chemical^2))

## -----------------------------------------------------------------------------
set.seed(0)
n1 = 20;n2 = 30
mu1 = mu2 = 0
sigma1 = sigma2 = 1
m = 1000;
alphahat1=alphahat2=numeric(m)

count5test = function(x, y) {
  X = x - mean(x)
  Y = y - mean(y)
  outx = sum(X > max(Y)) + sum(X < min(Y))
  outy = sum(Y > max(X)) + sum(Y < min(X))
  return(as.integer(max(c(outx, outy)) > 5))
}

count5test_permutation = function(z) {
  n = length(z)
  x = z[1:(n/2)]
  y = z[-(1:(n/2))]
  X = x - mean(x)
  Y = y - mean(y)
  outx = sum(X > max(Y)) + sum(X < min(Y)) 
  outy = sum(Y > max(X)) + sum(Y < min(X))
  return(as.integer(max(c(outx, outy)) > 5))
}

permutation = function(z,R) {
  n = length(z)
  out = numeric(R)
  for (r in 1: R){
    p = sample(1:n ,n ,replace = FALSE)
    out[r] = count5test_permutation(z[p])
  }
  sum(out)/R
}              

for (i in 1:m){
  x = rnorm(n1, mu1, sigma1)
  y = rnorm(n2, mu2, sigma2)
  x = x - mean(x) 
  y = y - mean(y)
  z = c(x,y)
  alphahat1[i]=count5test(x, y)
  alphahat2[i]=permutation(z,1000)
}

alphahat11=mean(alphahat1)
alphahat22=mean(alphahat2)
round(c(count5test=alphahat11,count5test_permutation=alphahat22),6)

## ----warning=FALSE------------------------------------------------------------
library(energy)
library(Ball)
library(RANN)
library(boot)

Tn <- function(z, ix, sizes,k) {
  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0);
  z <- z[ix, ];
  NN <- nn2(data=z, k=k+1)  
  block1 <- NN$nn.idx[1:n1,-1]
  block2 <- NN$nn.idx[(n1+1):n,-1]
  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
  (i1 + i2) / (k * n)
}
eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R,
                   sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}


# 1.Unequal variances and equal expectations
mu1 <- c(0,0,0)
sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
mu2 <- c(0,0,0)
sigma2 <- matrix(c(2,0,0,0,3,0,0,0,4),nrow=3,ncol=3)
n1=n2=20
n <- n1+n2 
N = c(n1,n2)
k=3
R=999
m=100
set.seed(1234)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  mydata1 <- mvrnorm(n1,mu1,sigma1)
  mydata2 <- mvrnorm(n2,mu2,sigma2)
  mydata <- rbind(mydata1,mydata2)
  p.values[i,1] <- eqdist.nn(mydata,N,k)$p.value
  p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
pow


# 2. Unequal variances and unequal expectations
mu1 <- c(0,0,0)
sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
mu2 <- c(0.5,-0.5,0.5)
sigma2 <- matrix(c(2,0,0,0,2,0,0,0,2),nrow=3,ncol=3)
n1=n2=20
n <- n1+n2 
N = c(n1,n2)
k=3
R=999
m=100
set.seed(1234)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  mydata1 <- mvrnorm(n1,mu1,sigma1)
  mydata2 <- mvrnorm(n2,mu2,sigma2)
  mydata <- rbind(mydata1,mydata2)
  p.values[i,1] <- eqdist.nn(mydata,N,k)$p.value
  p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
pow


# 3.Non-normal distributions
n1=n2=20
n <- n1+n2 
N = c(n1,n2)
k=3
R=999
m=100
set.seed(1234)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  mydata1 <- as.matrix(rt(n1,1,2),ncol=1)
  mydata2 <- as.matrix(rt(n2,2,5),ncol=1)
  mydata <- rbind(mydata1,mydata2)
  p.values[i,1] <- eqdist.nn(mydata,N,k)$p.value
  p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
pow

n1=n2=20
n <- n1+n2 
N = c(n1,n2)
k=3
R=999
m=100
set.seed(1234)
p.values <- matrix(NA,m,3)
rbimodel<-function(n,mu1,mu2,sd1,sd2){
  index=sample(1:2,n,replace=TRUE)
  x=numeric(n)
  index1<-which(index==1)
  x[index1]<-rnorm(length(index1), mu1, sd1)
  index2<-which(index==2)
  x[index2]<-rnorm(length(index2), mu2, sd2)
  return(x)
}
for(i in 1:m){
  mydata1 <- as.matrix(rbimodel(n1,0,0,1,2),ncol=1)
  mydata2 <- as.matrix(rbimodel(n2,1,1,4,3),ncol=1)
  mydata <- rbind(mydata1,mydata2)
  p.values[i,1] <- eqdist.nn(mydata,N,k)$p.value
  p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
pow


# 4.Unbalanced samples
mu1 <- c(0,0,0)
sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
mu2 <- c(0.5,-0.5,0.5)
sigma2 <- matrix(c(2,0,0,0,2,0,0,0,2),nrow=3,ncol=3)
n1=10
n2=100
n <- n1+n2 
N = c(n1,n2)
k=3
R=999
m=100
set.seed(1234)
p.values <- matrix(NA,m,3)
for(i in 1:m){
  mydata1 <- mvrnorm(n1,mu1,sigma1)
  mydata2 <- mvrnorm(n2,mu2,sigma2)
  mydata <- rbind(mydata1,mydata2)
  p.values[i,1] <- eqdist.nn(mydata,N,k)$p.value
  p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
pow


## -----------------------------------------------------------------------------
set.seed(0)

lap_f = function(x) 0.5*exp(-abs(x))

rw.Metropolis = function(sigma, x0, N){
  x = numeric(N)
  x[1] = x0
  u = runif(N)
  k = 0
  for (i in 2:N) {
    y = rnorm(1, x[i-1], sigma)
    if (u[i] <= (lap_f(y) / lap_f(x[i-1]))) x[i] = y 
    else {
      x[i] = x[i-1]
      k = k+1
    }
  }
  return(list(x = x, k = k))
}

N = 2000
sigma = c(.05, .5, 2, 16)
x0 = 25
rw1 = rw.Metropolis(sigma[1],x0,N)
rw2 = rw.Metropolis(sigma[2],x0,N)
rw3 = rw.Metropolis(sigma[3],x0,N)
rw4 = rw.Metropolis(sigma[4],x0,N)

#plot
# par(mfrow=c(2,2))  #display 4 graphs together
rw = cbind(rw1$x, rw2$x, rw3$x,  rw4$x)
for (j in 1:4) {
  plot(rw[,j], type="l",
       xlab=bquote(sigma == .(round(sigma[j],3))),
       ylab="X", ylim=range(rw[,j]))
}

#number of candidate points rejected
Rej = cbind(rw1$k, rw2$k, rw3$k, rw4$k)
Acc = round((N-Rej)/N,4)
rownames(Acc) = "Accept rates"
colnames(Acc) = paste("sigma=",sigma)
knitr::kable(Acc)

## ----warning=FALSE------------------------------------------------------------
set.seed(0)
library(GeneralizedHyperbolic)
Gelman.Rubin = function(psi) {
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i-th row of X
  psi = as.matrix(psi)
  n = ncol(psi)
  k = nrow(psi)
  psi.means = rowMeans(psi)        #row means
  B = n * var(psi.means)           #between variance est.
  psi.w = apply(psi, 1, "var")     #within variances
  W = mean(psi.w)                  #within est.
  v.hat = W*(n-1)/n + (B/n)        #upper variance est.
  r.hat = v.hat / W                #G-R statistic
  return(r.hat)
}

laplace.chain = function(sigma, N, X1) {
  #generates a Metropolis chain for Normal(0,1)
  #with Normal(X[t], sigma) proposal distribution
  #and starting value X1
  x = rep(0, N)
  x[1] = X1
  u = runif(N)
  param = c(0, 1, 1)
  for (i in 2:N) {
    xt = x[i-1]
    y = rnorm(1, xt, sigma) #candidate point
    r1 = dskewlap(y, param = param) * dnorm(xt, y, sigma)
    r2 = dskewlap(xt, param = param) * dnorm(y, xt, sigma)
    r = r1 / r2
    if (u[i] <= r) x[i] = y else
      x[i] = xt
  }
  return(x)
}
sigma =1            #parameter of proposal distribution
k = 4               #number of chains to generate
n = 15000           #length of chains
b = 1000            #burn-in length

#choose overdispersed initial values
x0 = c(-10, -5, 5, 10)

#generate the chains
X = matrix(0, nrow=k, ncol=n)
for (i in 1:k)
  X[i, ] = laplace.chain(sigma, n, x0[i])

#compute diagnostic statistics
psi = t(apply(X, 1, cumsum))
for (i in 1:nrow(psi)){
  psi[i,] = psi[i,] / (1:ncol(psi))
}
print(Gelman.Rubin(psi))

#plot psi for the four chains
#par(mfrow=c(2,2))
for (i in 1:k){
  plot(psi[i, (b+1):n], type="l",xlab=i, ylab=bquote(psi))
}


par(mfrow=c(1,1)) #restore default
#plot the sequence of R-hat statistics
rhat = rep(0, n)
for (j in (b+1):n){rhat[j] = Gelman.Rubin(psi[,1:j])}
plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
k=c(4:25,100,500,1000)

S = function(a,k){
  ck = sqrt(a^2*k/(k+1-a^2))
  pt(ck,df=k,lower.tail=FALSE)
}

solve = function(k){
  output = uniroot(function(a){S(a,k)-S(a,k-1)},lower=1,upper=2)
  output$root
}

root = matrix(0,2,length(k))

for (i in 1:length(k)){
  root[2,i]=round(solve(k[i]),4)
}

root[1,] = k
rownames(root) = c('k','A(k)')
root

## ----warning=FALSE------------------------------------------------------------
library(nloptr)
# Mle 
eval_f0 = function(x,x1,n.A=444,n.B=132,nOO=361,nAB=63) {
  
  r1 = 1-sum(x1)
  nAA = n.A*x1[1]^2/(x1[1]^2+2*x1[1]*r1)
  nBB = n.B*x1[2]^2/(x1[2]^2+2*x1[2]*r1)
  r = 1-sum(x)
  return(-2*nAA*log(x[1])-2*nBB*log(x[2])-2*nOO*log(r)-
           (n.A-nAA)*log(2*x[1]*r)-(n.B-nBB)*log(2*x[2]*r)-nAB*log(2*x[1]*x[2]))
}


# constraint
eval_g0 = function(x,x1,n.A=444,n.B=132,nOO=361,nAB=63) {
  return(sum(x)-0.999999)
}

opts = list("algorithm"="NLOPT_LN_COBYLA","xtol_rel"=1.0e-8)

mle = NULL
r = matrix(0,1,2)
r = rbind(r,c(0.25,0.25))        # the beginning value of p0 and q0
j = 2
while (sum(abs(r[j,]-r[j-1,]))>1e-8) {
  res = nloptr( x0=c(0.3,0.1),
                eval_f=eval_f0,
                lb = c(0,0), ub = c(1,1), 
                eval_g_ineq = eval_g0, 
                opts = opts, x1=r[j,],n.A=444,n.B=132,nOO=361,nAB=63 )
  j = j+1
  r = rbind(r,res$solution)
  mle = c(mle,-eval_f0(x=r[j,],x1=r[j-1,]))
}
#the result of EM algorithm
r=cbind(r,mle)
colnames(r)=c("p","q","log_Mle")
r

#the log_max likelihood values
plot(mle,type = 'l')

## ----warning=FALSE------------------------------------------------------------
attach(mtcars)

formulas = list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)
#1 for loops
f1 = vector("list", length(formulas))
for (i in seq_along(formulas)){
  f1[[i]] = lm(formulas[[i]], data = mtcars)
}
f1
#2 lapply
f2 = lapply(formulas, function(x) lm(formula = x, data = mtcars))
f2

## ----warning=FALSE------------------------------------------------------------
set.seed(0)
trials = replicate(
  100,
  t.test(rpois(10, 10), rpois(7, 10)),
  simplify = FALSE
)
# anonymous function:
sapply(trials, function(x) x[["p.value"]])

## extra challenge:
sapply(trials, "[[" ,3)

## ----warning=FALSE------------------------------------------------------------
attach(mtcars)
lapply1=function(data,funct,output_type){
  new=Map(funct,data)
  vapply(new,function(x) x , output_type)
}

##   Example
lapply1(mtcars,median,double(1))

## ----warning=FALSE------------------------------------------------------------
set.seed(0)
library(Rcpp)
library(microbenchmark)

# R
lap_f = function(x) exp(-abs(x))

rwR_Metropolis = function(sigma, x0, N){
  x = numeric(N)
  x[1] = x0
  u = runif(N)
  k = 0

  for (i in 2:N) {
    y = rnorm(1, x[i-1], sigma)
    if (u[i] <= (lap_f(y) / lap_f(x[i-1]))) x[i] = y 
    else {
      x[i] = x[i-1]
      k = k+1
      }
    }
    return(list(x = x, k = k))
  }

dir_cpp = 'C:/Users/admin/Desktop/StatComp20068/src/'
sourceCpp(paste0(dir_cpp,"rwC_Metropolis.cpp"))
x0 = 15
N = 2000
sigma = 2
rwR=rwR_Metropolis(sigma,x0,N)
rwC=rwC_Metropolis(sigma,x0,N)

#par(mfrow=c(1,2))
plot(rwR$x,type = "l",xlab =bquote(sigma==2),ylab = "rwR",ylim = range(rwR$x))
plot(rwC,type = "l",xlab =bquote(sigma==2),ylab = "rwC",ylim = range(rwC))

## -----------------------------------------------------------------------------
set.seed(0)
rwR = rwR_Metropolis(sigma,x0,N)$x[-(1:100)]
rwC = rwC_Metropolis(sigma,x0,N)[-(1:100)]
qqplot(rwR,rwC)
abline(a=0,b=1,col='black')

## -----------------------------------------------------------------------------
set.seed(0)
(time = microbenchmark(rwR=rwR_Metropolis(sigma,x0,N),rwC=rwC_Metropolis(sigma,x0,N)))

