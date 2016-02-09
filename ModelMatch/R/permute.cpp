#include <Rcpp.h>
using namespace Rcpp;

// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
inline int randWrapper(const int n) { return floor(unif_rand()*n); }

// [[Rcpp::export]]
Rcpp::NumericVector randomShuffle(Rcpp::NumericVector a) {

  // clone a into b to leave a alone
  Rcpp::NumericVector b = Rcpp::clone(a);
  std::random_shuffle(b.begin(), b.end(), randWrapper);
  return b;
}

// [[Rcpp::export]]
List permute_within_groups_cpp(List strata){
  List permuted = clone(strata);
  int G = strata.size();

  for(int g = 0; g < G; g++){
    Rcpp::List tmp = strata(g);
    tmp(2) = randomShuffle(tmp(2)); // Third column is treatment
    permuted[g] = tmp;
  }
  return permuted;
}

// [[Rcpp::export]]
double meanC(NumericVector x) {
  int n = x.size();
  double total = 0;

  for(int i = 0; i < n; ++i) {
    total += x[i];
  }
  return total / n;
}


// [[Rcpp::export]]
NumericVector within_group_mean_cpp(List strata, NumericVector prediction, NumericVector response, float shift = 0){
  NumericVector resid = prediction-response;
  int G = strata.size();
  NumericVector means(G);

  for(int g = 0; g < G; g++){
    Rcpp::List tmp = strata(g);
    NumericVector treat_ind = tmp(2); // Third column is treatment
    int ng = treat_ind.size();
    NumericVector index = tmp(0); // First column is indices
    float tt = 0;
    float cc = 0;
    float nt = 0;
    float nc = 0;
    for(int nn = 0; nn < ng; nn++){
      if (treat_ind(nn) == 1){
        tt = tt + resid[index(nn)];
        nt = nt+1;
      } else {
        cc = cc + resid[index(nn)];
        nc = nc+1;
      }
    }
    means(g) = cc/nc - tt/nt + shift;
  }
  return means;
}
