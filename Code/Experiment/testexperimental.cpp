#include <Rcpp.h>
using namespace Rcpp;


#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix sweep_cpp(NumericMatrix x,
                      NumericVector lastCol) {
        
        int n = x.nrow();
        int m = x.ncol();
        
        NumericMatrix res(n, m);
        int i, j;
        
        for (j = 0; j < m; j++) {
                for (i = 0; i < n; i++) {
                        res(i, j) =  x(i, j) * lastCol[j];
                }
        }
        
        return res;
}


// [[Rcpp::export]]
List discrete_trito_model(NumericMatrix HS, 
                          NumericMatrix HI,
                          NumericMatrix HR,
                          NumericMatrix PS,
                          NumericMatrix PI,
                          NumericMatrix SS,
                          NumericMatrix SI,
                          NumericMatrix adj,
                          NumericVector param
                        ) {
        
        double b_H = param["b_H"];
        double b_P = param["b_P"];
        double b_S = param["b_S"];
        double mu_H = param["mu_H"];
        double theta_P = param["theta_P"];
        double theta_S = param["theta_S"];
        double theta_H = param["theta_H"];
        double f_P = param["f_P"];
        double f_S = param["f_S"];
        double gamma = param["gamma"];
        double c_PP = param["c_PP"];
        double c_SS = param["c_SS"];
        double c_PS = param["c_PS"];
        double c_SP = param["c_SP"];
        double d = param["d"];
        
 NumericVector NH =  HS( 0 , _ ) + HI(0, _) + HR(0,_); 
 NumericVector NP = PS( 0 , _ ) + PI(0, _);
 NumericVector NS = SS( 0 , _ ) + SI(0, _);
 
 NumericVector HP_ratio = PI(0,_)/NH;
 NumericVector HS_ratio = SI(0,_)/NH;
 NumericVector H_ratio = HI(0,_)/NH;
 
            
 NumericVector infections_H = (theta_P * f_P * HP_ratio) + (theta_S * f_S * HS_ratio);
 NumericVector infections_P = (theta_H * f_P) *  H_ratio;
 NumericVector infections_S = (theta_H * f_S) *  H_ratio;

 NumericVector HS_Rates = (b_H * NH) - (infections_H + mu_H) * HS(0,_);
 NumericVector HI_Rates = (infections_H * HS(0,_)) - (gamma + mu_H) * HI(0,_);
 NumericVector HR_Rates = (gamma * HI(0,_)) - (mu_H * HR(0,_));
 
 
 NumericVector disp = rowSums(sweep_cpp(adj, (- d* PS(0,_)))) + 
         transpose(adj) %*% (d * PS_mat[j, ]) ;
        
 
 
 
 
 
 return Rcpp::List::create(Rcpp::Named("HS") = HS_Rates,
                           Rcpp::Named("HI") = HI_Rates,
                           Rcpp::Named("HR") = HR_Rates,
                           Rcpp::Named("d") = disp);
}
        

/*** R
param_ALT<- c(b_H = 1/100,
              b_P = 1/100,
              b_S = 1/100,
              mu_H = 1/100,
              f_P = 0.050, # Biting rate of the p. vector
              f_S = 0.030, # Biting rate of the s.vector
              theta_P = 0.50, # Transmission probability of p. vector
              theta_S = 0.30, theta_H = 0.50,
              gamma = 1/100, c_PP = 0, c_SS = 0, c_PS = 0,c_SP=0,d = 0.01) # Transmission probability of s. vector)


HS = matrix(1, nrow = 3, ncol = 3)
HI = matrix(2, nrow = 3, ncol = 3)
HR = matrix(3, nrow = 3, ncol = 3)

PS = matrix(1, nrow = 3, ncol = 3)
PI = matrix(2, nrow = 3, ncol = 3)

SS = matrix(1, nrow = 3, ncol = 3)
SI = matrix(2, nrow = 3, ncol = 3)

adj = matrix(c(0,1,1,1,0,1,1,1,0), nrow =3, ncol =3)


discrete_trito_model(HS, HI, HR,PS, PI, SS,SI, adj,param_ALT)
rcppFun(PS, 1/5 * P[1,])
*/


