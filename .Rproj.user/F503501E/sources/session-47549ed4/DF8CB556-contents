#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List discrete_trito_model_rcpp_ONEPATCH(
                arma::mat HS, //Susceptible human matrix
                arma::mat HI, //Infectious human matrix
                arma::mat HR, //Recovered human matrix
                arma::mat PS, //Susceptible p.vector matrix
                arma::mat PI, // Infectious p.vector matrix
                arma::mat SS, // Susceptible s.vector matrix
                arma::mat SI, // Infectious s.vector matrix
                arma::mat adj, // Adjacency matrix
                Rcpp::NumericVector param){ //Parameter data.frame 
        
        double b_H = param["b_H"]; // Birth rate of human
        double b_P = param["b_P"]; // Birth rate of p.vector
        double b_S = param["b_S"]; // Birth rate of s.vector
        double mu_H = param["mu_H"]; // Death rate of human
        double theta_P = param["theta_P"]; // Transmission prob of p.vector
        double theta_S = param["theta_S"]; // Transmission prob of s.vector
        double theta_H = param["theta_H"]; // Transmission prob of human
        double f_P = param["f_P"]; //Biterate of p.vector
        double f_S = param["f_S"];  //Bite rate of s.vector
        double gamma = param["gamma"]; //Recovery rate for human
        double c_PP = param["c_PP"]; //Competition coefficient for primary on primary
        double c_SS = param["c_SS"]; //Competition coefficient for secondary on secondary
        double c_PS = param["c_PS"];  //Competition coefficient for primary on secondary
        double c_SP = param["c_SP"]; //Competition coefficient for secondary on primary
        int ntime = param["ntime"] - 1; //How long to simulate for (stupid indexing)
        int disturbance_time = param["disturbance_time"] - 1; //When is disturbance time?
        double delta_T = param["delta_T"]; //Time steps
        double mortality_P = param["mortality_P"]; //Surviviorship of p.vector
        double mortality_S = param["mortality_S"]; //Survivorship of s.vector
        
        //The number of patches is dependent on either the column or row (should 
        //be equal)
        int patch_num = adj.n_rows;
        
        for(int j = 0; j < ntime; j++){
                
                //Total populations of human, primary, and secondary
                arma::rowvec NH = HS.row(j) + HI.row(j) + HR.row(j); 
                arma::rowvec NP = PS.row(j) + PI.row(j);
                arma::rowvec NS = SS.row(j) + SI.row(j);
                
                //Ratios of vector to human and infected human to total human
                arma::rowvec HP_ratio = PI.row(j)/NH;
                arma::rowvec HS_ratio = SI.row(j)/NH;
                arma::rowvec H_ratio =  HI.row(j)/NH;
                
                //Primary and secondary will infect humans
                arma::rowvec infections_H = (theta_P * f_P * HP_ratio) + (theta_S * f_S * HS_ratio);
                
                //Primary and secondary will be infected by human
                arma::rowvec infections_P = (theta_H * f_P) *  H_ratio;
                arma::rowvec infections_S = (theta_H * f_S) *  H_ratio;
                
                //Demographic/infection rates of humans
                arma::rowvec HS_Rates = (b_H * NH) - (infections_H + mu_H) % HS.row(j);
                arma::rowvec HI_Rates = (infections_H % HS.row(j)) - (gamma + mu_H) * HI.row(j);
                arma::rowvec HR_Rates = (gamma * HI.row(j)) - (mu_H * HR.row(j));
                
                //Demographic/infection rates of primary vector
                arma::rowvec PS_Rates = (b_P * NP) - (infections_P + 
                        (c_SP * NS) + (c_PP * NP)) % PS.row(j); 
                
                arma::rowvec PI_Rates = (infections_P % PS.row(j)) - (
                        (c_SP * NS) + (c_PP * NP)) % PI.row(j); 
                
                //Demographic/infection of secondary vector
                arma::rowvec  SS_Rates = (b_S * NS) - (infections_S +
                        (c_PS * NP) + (c_SS * NS)) % SS.row(j);
                
                arma::rowvec  SI_Rates = (infections_S % SS.row(j)) - (
                        (c_PS * NP) + (c_SS * NS)) % SI.row(j);
                
                
              
                //Let's look what the new population size will be
                arma::rowvec total_change_HS = HS.row(j) + (HS_Rates *  delta_T);
                arma::rowvec total_change_HI = HI.row(j) + (HI_Rates *  delta_T);
                arma::rowvec total_change_HR = HR.row(j) + (HR_Rates *  delta_T);
                arma::rowvec total_change_PS = PS.row(j) + (PS_Rates * delta_T) ;
                arma::rowvec total_change_PI = PI.row(j) + (PI_Rates *  delta_T);
                arma::rowvec total_change_SS = SS.row(j) + (SS_Rates *  delta_T);
                arma::rowvec total_change_SI = SI.row(j) + (SI_Rates *  delta_T);
                
                total_change_HS = arma::clamp(total_change_HS, 0, total_change_HS.max());
                total_change_HI = arma::clamp(total_change_HI, 0, total_change_HI.max());
                total_change_HR = arma::clamp(total_change_HR, 0, total_change_HR.max());
                total_change_PS = arma::clamp(total_change_PS, 0, total_change_PS.max());
                total_change_PI = arma::clamp(total_change_PI, 0, total_change_PI.max());
                total_change_SS = arma::clamp(total_change_SS, 0, total_change_SS.max());
                total_change_SI = arma::clamp(total_change_SI, 0, total_change_SI.max());
                
                
                if (j != disturbance_time) {
                        
                        // We then update these row by row
                        HS.row(j + 1) = total_change_HS;
                        HI.row(j + 1) =  total_change_HI;
                        HR.row(j + 1) =  total_change_HR;
                        PS.row(j + 1) =  total_change_PS;
                        PI.row(j + 1) =  total_change_PI;
                        SS.row(j + 1) =  total_change_SS;
                        SI.row(j + 1) =  total_change_SI;
                }
                
                else if (j == disturbance_time) {
                       
        
                        HS.row(j + 1) = total_change_HS;
                        HI.row(j + 1) = total_change_HI;
                        HR.row(j + 1) = total_change_HR;
                        PS.row(j + 1) = total_change_PS * mortality_P;;
                        PI.row(j + 1) = total_change_PI * mortality_P;;
                        SS.row(j + 1) = total_change_SS * mortality_S;;
                        SI.row(j + 1) = total_change_SI * mortality_S;;
                        
                }
        }
        
        
        
        return Rcpp::List::create(Rcpp::Named("HS") = HS,
                                  Rcpp::Named("HI") = HI,
                                  Rcpp::Named("HR") = HR,
                                  Rcpp::Named("PS") = PS,
                                  Rcpp::Named("PI") = PI,
                                  Rcpp::Named("SS") = SS,
                                  Rcpp::Named("SI") = SI);
        
        
}

