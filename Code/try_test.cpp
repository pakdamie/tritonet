#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
//[[Rcpp::export]]
NumericMatrix tryo(int num_patch, 
                   NumericMatrix v) {
        
        
        NumericVector tryit = {2,2,2};
        
        arma::mat trymat(3, 3);
        for (int p = 0; p < 3; p++){
                trymat (p, p) = tryit[p];
        }
        
        arma::mat v2 = Rcpp::as<arma::mat>(v);
        
        NumericMatrix disp_mat = Rcpp::wrap(v2 * trymat);
        
        return (transpose(disp_mat) );
        

}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
tryo(3,a)
*/
