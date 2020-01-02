#include <Rcpp.h>
using namespace Rcpp;

//' @title A random walk Metropolis sampler
//' @description A random walk Metropolis sampler for generating the standard Laplace distribution
//' @param x0: the initial point(double)
//' @param sigma: the standard deviation in the normal distribution(double)
//' @param N: the length of the chain(int)
//' @import Rcpp
//' @useDynLib SC19083
//' @return a random sample of size N
//' @examples
//' \dontrun{
//' N<-100
//' sigma<-1
//' x0<-0
//' Metropolis(sigma,x0,N)
//' }
//' @export
// [[Rcpp::export]]
NumericVector Metropolis(double sigma, double x0, int N){
  NumericVector x(N);
  x[0]=x0;
  for (int i=1;i<N;i++) {
    double e=runif(1)[0];
    double z=rnorm(1,x[i-1],sigma)[0];
    if (e<=exp(abs(x[i-1])-abs(z))) x[i]=z;
    else {
      x[i]=x[i-1];
    }
  }
  NumericVector out=x;
  return (out);
}