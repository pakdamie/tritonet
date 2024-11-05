#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

#include <RcppArmadilloExtensions/sample.h>

arma::mat timesTwo(arma::mat adj, arma::mat x) {
        arma::rowvec disp0 = -0.1 * x.row(0);
        arma::mat disp = (adj.each_row() % (-0.1 * x.row(0)));
        arma::rowvec tryit = adj.row(0) %  (-0.1 * x.row(0));
        return(disp);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
a <- timesTwo(adjacency_matrix_adj, PS_mat)
*/
