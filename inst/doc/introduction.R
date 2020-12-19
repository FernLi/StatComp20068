## -----------------------------------------------------------------------------
find_k1k2_1=function(n,N,k,a){
  s=0;t=0;
  repeat{
    s=s+choose(k,t)*choose(N-k,n-t)/choose(N,n)
    if(s>a/2){break;}
    t=t+1;
  }
  k1=t;
  
  t=n;m=0;
  repeat{
    m=m+choose(k,t)*choose(N-k,n-t)/choose(N,n)
    if(m>a/2){break;}
    t=t-1;
  }
  k2=t;vec=c(k1,k2);
  return(vec);
}

find_k1k2_2=function(p,n,a){
  k1=0;k2=0;s=0;t=0;k=0;
  repeat{
    s=s+choose(n,k)*p^(k)*(1-p)^(n-k);
    if(s>a/2){break;}
    k=k+1;
  }
  k1=k;
  
  k=n;
  repeat{
    t=t+choose(n,k)*p^(k)*(1-p)^(n-k);
    if(t>a/2){ break;}
    k=k-1;
  }
  k2=k;vec=c(k1,k2);
  return(vec);
}
f1<-function(n,N,k,t0,p){
  p=0;t=0
  repeat{
    p=p+(choose(k,t)*choose(N-k,n-t))/choose(N,n)
    if(t>t0){
      return(p)
      break
    }
    t=t+1
  }
}

f2<-function(p,n,t){
  s=0;k=0
  repeat{
    s=s+choose(n,k)*(p^k)*(1-p)^(n-k)
    if(k>t){
      return(s)
      break
    }
    k=k+1
  }
}


## -----------------------------------------------------------------------------
#小样本下的总体比例检验
N=300;n=50;t=10;theta=0.4;k=N*theta;a=0.05
find_k1k2_1(n,N,k,a)
p1=f1(n,N,k,t);p1
if(p1<a){
  print("拒绝H0")
}else{
  print("不拒绝H0")
}
#大样本下的总体比例检验
N=10000;n=1000;t=35;theta=0.038;k=N*theta
find_k1k2_2(theta,n,a)
p2=f2(theta,n,t);p2
if(p2<a){
  print("拒绝H0")
}else{
  print("不拒绝H0")
}

## -----------------------------------------------------------------------------
Chisquare_test<-function(x,nrow,ncol){
  n=sum(x);
  ni<-rep(0,nrow);nj<-rep(0,ncol)
  for(i in 1:nrow){
    for(j in 1:ncol){
      ni[i]=x[i,j]+ni[i]
    }
  }
  for(j in 1:ncol){
    for(i in 1:nrow){
      nj[j]=x[i,j]+nj[j]
    }
  }
  e<-matrix(nrow=nrow,ncol=ncol);x_square=0;
  for(i in 1:nrow){
    for(j in 1:ncol){
      e[i,j]=ni[i]*nj[j]/n
      x_square=((x[i,j]-e[i,j])^2/e[i,j])+x_square
    }
  }
  y=list("Total value of each row"=ni,"Total value of each column"=nj,"Expected frequency"=e,"Chi-square sum"=x_square)
  return(y)
}

## -----------------------------------------------------------------------------
#Test whether different viewers pay the same attention to the three types of programs  Consistency test
x=matrix(c(83,70,45,91,86,15,41,38,10),nrow=3,ncol=3)
Chisquare_test(x,nrow=3,ncol=3)
chisq.test(x,correct=TRUE)
#Test whether there is an association between blood type and liver disease      Independence test
y<-matrix(c(98,67,13,18,38,41,8,12,289,262,57,30),nrow=4,ncol=3)
Chisquare_test(y,nrow=4,ncol=3)
chisq.test(y,correct=TRUE)

