#include <Rcpp.h>
using namespace Rcpp;

//' @title A Gibbs sampler using Rcpp for a bivariate density
//' @description A Gibbs sampler using Rcpp
//' @param N the number of samples
//' @param n the number of between-sample random numbers
//' @param a the parameter of bivariate density
//' @param b the parameter of bivariate density
//' @return a random sample of size n
//' @examples
//' \dontrun{
//' set.seed(1211)
//' gibb_C<-gibbsC(120,15,1,1)
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix gibbsC(int N, int n,int a,int b) {
  NumericMatrix mat(N, 2);
  float x = 0, y = 0;
  mat(0,0)=1.0; 
  mat(0,1)=0.5;
  for(int i = 1; i < N; i++) {
    y=mat(i-1,1);
    mat(i,0)=rbinom(1, n, y)[0];
    x=mat(i,0);
    mat(i,1)=rbeta(1, x+a, n-x+b)[0];
  }
  return(mat);
}