#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List discrete_trito_model_rcpp_ONEPATCH(
  arma::mat HS, // Human (Susceptible) Matrix
  arma::mat HI, // Human (Infectious) Matrix
  arma::mat HR, // Human (Recovered) Matrix
  arma::mat PS, // Primary Vector (Susceptible) Matrix
  arma::mat PI, // Primary Vector (Infected) Matrix
  arma::mat MS, // Secondary Vector (Susceptible) Matrix
  arma::mat MI, // Secondary Vector (Susceptible) Matrix
  Rcpp::NumericVector param){ //Parameter data.frame 
        
  double b_H = param["b_H"]; // Birth rate of human
  double b_P = param["b_P"]; // Birth rate of p.vector
  double b_M = param["b_M"]; // Birth rate of s.vector
  double mu_H = param["mu_H"]; // Death rate of human
  double theta_P = param["theta_P"]; // Transmission prob of p.vector
  double theta_M = param["theta_M"]; // Transmission prob of s.vector
  double theta_H = param["theta_H"]; // Transmission prob of human
  double f_P = param["f_P"]; //Biterate of p.vector
  double f_M = param["f_M"];  //Bite rate of s.vector
  double gamma = param["gamma"]; //Recovery rate for human
  double c_PP = param["c_PP"]; //Competition coefficient for primary on primary
  double c_MM = param["c_MM"]; //Competition coefficient for secondary on secondary
  double c_PM = param["c_PM"];  //Competition coefficient for primary on secondary
  double c_MP = param["c_MP"]; //Competition coefficient for secondary on primary
  int ntime = param["ntime"] - 1; //How long to simulate for (stupid indexing)
  int disturbance_time = param["disturbance_time"] - 1; //When is disturbance time?
  double delta_T = param["delta_T"]; //Time steps
  double mortality_P = param["mortality_P"]; //Surviviorship of p.vector
  double mortality_M = param["mortality_M"]; //Survivorship of s.vector
        
  for(int j = 0; j < ntime; j++){
                
    // Total populations of human, primary, and secondary at time j
    // for rate calculations 
    arma::rowvec NH = HS.row(j) + HI.row(j) + HR.row(j); 
    arma::rowvec NP = PS.row(j) + PI.row(j);
    arma::rowvec NM = MS.row(j) + MI.row(j);
    
    // Ratios of vector to human and infected human to total human
    arma::rowvec HP_ratio = PI.row(j)/NH;
    arma::rowvec HM_ratio = MI.row(j)/NH;
    arma::rowvec H_ratio =  HI.row(j)/NH;
    
    // Primary and secondary will infect humans
    arma::rowvec infections_H_fromP = (theta_P * f_P * HP_ratio);
    arma::rowvec infections_H_fromM = (theta_M * f_M * HM_ratio);
    arma::rowvec infections_H = infections_H_fromP + infections_H_fromM;
    
    //Primary and secondary will be infected by human
    arma::rowvec infections_P = (theta_H * f_P) *  H_ratio;
    arma::rowvec infections_M = (theta_H * f_M) *  H_ratio;
    
    //Demographic/infection rates of humans
    arma::rowvec HS_Rates = (b_H * NH) - (infections_H + mu_H) % HS.row(j);
    arma::rowvec HI_Rates = (infections_H % HS.row(j)) - (gamma + mu_H) * HI.row(j);
    arma::rowvec HR_Rates = (gamma * HI.row(j)) - (mu_H * HR.row(j));
    
    //Demographic/infection rates of primary vector
    arma::rowvec PS_Rates = (b_P * NP) - 
      (infections_P + (c_MP * NM) + (c_PP * NP)) % PS.row(j); 
    
    arma::rowvec PI_Rates = (infections_P % PS.row(j)) - 
       ((c_MP * NM) + (c_PP * NP)) % PI.row(j); 
    
    //Demographic/infection of secondary vector
    rma::rowvec  MS_Rates = (b_M * NM) - 
      (infections_M + (c_PM * NP) + (c_MM * NM)) % MS.row(j);
    
    arma::rowvec  MI_Rates = (infections_M % MS.row(j)) - 
      ((c_PM * NP) + (c_MM * NM)) % MI.row(j);
    
    //The total new changes
    arma::rowvec total_change_HS = HS.row(j) + (HS_Rates *  delta_T);
    arma::rowvec total_change_HI = HI.row(j) + (HI_Rates *  delta_T);
    arma::rowvec total_change_HR = HR.row(j) + (HR_Rates *  delta_T);
    arma::rowvec total_change_PS = PS.row(j) + (PS_Rates *  delta_T);
    arma::rowvec total_change_PI = PI.row(j) + (PI_Rates *  delta_T);
    arma::rowvec total_change_MS = MS.row(j) + (MS_Rates *  delta_T);
    arma::rowvec total_change_MI = MI.row(j) + (MI_Rates *  delta_T);
    
    //A bit of a wonky way to ensure that the total change is not 
    //"Clamp each element to the [min_val,max_val] interval;
    //any value lower than min_val will be set to min_val, 
    //and any value higher than max_val will be set to max_val"
    
    //I chose an unbelievable value (1e30 for the maximum)
    
    total_change_HS = arma::clamp(total_change_HS, 0, 1e30);
    total_change_HI = arma::clamp(total_change_HI, 0, 1e30);
    total_change_HR = arma::clamp(total_change_HR, 0, 1e30);
    total_change_PS = arma::clamp(total_change_PS, 0, 1e30);
    total_change_PI = arma::clamp(total_change_PI, 0, 1e30);
    total_change_MS = arma::clamp(total_change_MS, 0, 1e30);
    total_change_MI = arma::clamp(total_change_MI, 0, 1e30);
               
               
    if (j != disturbance_time) {
                       
       // We then update these row by row
       HS.row(j + 1) = total_change_HS;
       HI.row(j + 1) =  total_change_HI;
       HR.row(j + 1) =  total_change_HR;
       PS.row(j + 1) =  total_change_PS;
       PI.row(j + 1) =  total_change_PI;
       MS.row(j + 1) =  total_change_MS;
       MI.row(j + 1) =  total_change_MI;
     }
               
      else if (j == disturbance_time) {
           HS.row(j + 1) = total_change_HS;
           HI.row(j + 1) = total_change_HI;
           HR.row(j + 1) = total_change_HR;
           PS.row(j + 1) = total_change_PS * mortality_P;
           PI.row(j + 1) = total_change_PI * mortality_P;
           MS.row(j + 1) = total_change_MS * mortality_M;
           MI.row(j + 1) = total_change_MI * mortality_M;
                        
                }
        }
        
        
        
        return Rcpp::List::create(
          Rcpp::Named("HS") = HS,
          Rcpp::Named("HI") = HI,
          Rcpp::Named("HR") = HR,
          Rcpp::Named("PS") = PS,
          Rcpp::Named("PI") = PI,
          Rcpp::Named("MS") = MS,
          Rcpp::Named("MI") = MI);
        
        
}

