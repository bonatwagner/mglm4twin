#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::sp_mat ef_variability_op(arma::sp_mat sensitivity, List WE, arma::vec k4, arma::vec W) {
  int n_par = WE.size();
  sp_mat V = sp_mat(n_par, n_par);
  for (int i = 0; i <= n_par - 1; i++) {
    vec Wi;
    Wi = as<sp_mat>(WE[i]).diag();
    for (int j = 0; j <= i; j++) {
      vec Wj;
      Wj = as<sp_mat>(WE[j]).diag();
      V(i,j) = -2 * sensitivity(i,j) + sum(k4 % Wi % W % Wj % W);
    }
  }
  return V;
}


//
//
// arma::sp_mat ef_variability_op(arma::sp_mat sensitivity, List W, arma::vec k4) {
//   int n_par = W.size();
//   sp_mat V = sp_mat(n_par, n_par);
//   for (int i = 0; i <= n_par - 1; i++) {
//     vec Wi;
//     Wi = as<sp_mat>(W[i]).diag();
//     for (int j = 0; j <= i; j++) {
//       vec Wj;
//       vec WW;
//       Wj = as<sp_mat>(W[j]).diag();
//       V(i,j) = -2 * sensitivity(i,j) + sum(k4 % Wi % Wj);
//     }
//   }
//   return V;
// }
