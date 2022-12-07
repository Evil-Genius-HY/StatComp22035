# include <Rcpp.h>
using namespace Rcpp;

//' @title PageRank Algorithm
//' @description PageRank Algorithm using Rcpp
//' @param MA Matrix of adjacency
//' @param d damping rate
//' @param MaxIter Maximum number of iteration
//' @param IniVec Initial iteration vector
//' @return Score result of each node
//' @examples
//' \dontrun{
//' NumNode <- 50
//' nrow <- ncol <- NumNode
//' data(MA)
//' IniVec <- runif(NumNode, 0, 1)
//' Score <- PageRank(MA, 0.85, 100000, IniVec)
//' data <- data.frame(
//' ndex = 1:NumNode,  
//' value = Score
//' )
//' ggplot(data, aes(x=index, y=value)) + 
//' geom_bar(stat = "identity", width=0.5)
//' }
//' @export
// [[Rcpp::export]]
NumericVector PageRank(NumericMatrix MA, double d, int MaxIter, NumericVector IniVec){
  
  int N = MA.nrow();
  NumericMatrix MT(N, N);
  NumericVector UpdVec(N);
  
  for(int i = 0; i < N; i++){
    double sum = 0;
    for(int j = 0; j < N; j++){
      sum = sum + MA(j,i);
    }
    if(sum>0){
      for(int k = 0; k < N; k++){
        MT(k,i) = MA(k,i)/sum;
      }
    }
  }
  
  for(int i = 0; i < MaxIter; i++){
    for(int j = 0; j < N; j++){
      UpdVec[j] = 0;
      for(int k = 0; k < N; k++){
        UpdVec[j] = UpdVec[j] + d*MT(j,k)*IniVec[k];
      }
      UpdVec[j] = UpdVec[j] + (1-d)/N;
    }
    double Diff = 0;
    for(int l = 0; l < N; l++){
      Diff = Diff + abs(UpdVec[l]-IniVec[l]);
    }
    if(Diff < 1e-10){
      break;
    }
    for(int m = 0; m < N; m++){
      IniVec[m] = UpdVec[m];
    }
  }
  
  return(UpdVec);
}
//'@useDynLib StatComp22035



//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp
//' @param N length of chain
//' @param mu parameter of bivariate normal distribution
//' @param sigma parameter of bivariate normal distribution
//' @param rho parameter of bivariate normal distribution
//' @param initial initial value
//' @return a random sample of size \code{N}
//' @export
// [[Rcpp::export]]
NumericMatrix gibbsC(NumericVector mu, NumericVector sigma, double rho,
                     NumericVector initial, int N) {
  NumericMatrix XY(N, 2);
  double x, y, m1, m2;
  XY(0, 0) = initial[0];
  XY(0, 1) = initial[1];
  for(int i = 1; i < N; i++) {
    y = XY(i - 1, 1);
    m1 = mu[0] + rho * (y - mu[1]) * sigma[0] / sigma[1];
    XY(i, 0) = rnorm(1, m1, sqrt(1 - rho * rho))[0] * sigma[0];
    x = XY(i, 0);
    m2 = mu[1] + rho * (x - mu[0]) * sigma[1] / sigma[0];
    XY(i, 1) = rnorm(1, m2, sqrt(1 - rho * rho))[0] * sigma[1];
  }
  return(XY);
}
//'@useDynLib StatComp22035