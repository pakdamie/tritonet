#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericMatrix tryo(int num_patch) {
        NumericVector tryit = {1,2,3};
        
        NumericMatrix trymat = ( num_patch);
        for (int p = 0; p < 3; p++){
                trymat (p, p) = tryit[p];
        }
        return trymat;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
tryo(4)
*/
