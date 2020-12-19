#' @title Chi-square statistics
#'
#' @param x Data matrix X
#' @param nrow Number of rows of data matrix X
#' @param ncol Number of columns of data matrix X
#'
#' @returny The value of the chi-square statistic
#' @export
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