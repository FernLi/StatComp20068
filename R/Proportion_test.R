#' @title Hypergeometric distribution interval estimation
#' @param n Sample size
#' @param N Total number
#' @param k k=N*theta
#' @param a Significance level
#'
#' @return vec
#' @export
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

#' @title Estimation of Hypergeometric Distribution Interval under Large Sample
#' @param p P value
#' @param n Sample size
#' @param a Significance level
#'
#' @return vec
#' @export
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

#' @title Hypergeometric distribution interval estimation
#' @param n Sample size
#' @param N Total number
#' @param k k=N*theta
#' @param t0 Number of data with specific attributes
#' @param p P value
#'
#' @return p
#' @export
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

#' @title Estimation of Hypergeometric Distribution Interval under Large Sample
#' @param p theta
#' @param n Sample size
#' @param t Number of data with specific attributes
#'
#' @return s
#' @export
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


