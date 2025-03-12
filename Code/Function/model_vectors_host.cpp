#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List model_vectors_host(
    arma::mat HS, // Human (Susceptible) Matrix
    arma::mat HI, // Human (Infectious) Matrix
    arma::mat HR, // Human (Recovered) Matrix
    arma::mat PS, // Primary Vector (Susceptible) Matrix
    arma::mat PI, // Primary Vector (Infected) Matrix
    arma::mat MS, // Secondary Vector (Susceptible) Matrix
    arma::mat MI, // Secondary Vector (Infected) Matrix
    Rcpp::NumericVector param // Parameter vector
) {
  
  // Parameters
  double b_H = param["b_H"]; // Human birth rate
  double b_P = param["b_P"]; // Primary vector birth rate
  double b_M = param["b_M"]; // Secondary vector birth rate
  double mu_H = param["mu_H"]; // Human death rate
  double theta_P = param["theta_P"]; // Primary vector transmission probability
  double theta_M = param["theta_M"]; // Secondary vector transmission probability
  double theta_H = param["theta_H"]; // Human transmission probability
  double f_P = param["f_P"]; // Primary vector bite rate
  double f_M = param["f_M"]; // Secondary vector bite rate
  double gamma = param["gamma"]; // Human recovery rate
  double c_PP = param["c_PP"]; // Competition coefficient: primary on primary
  double c_MM = param["c_MM"]; // Competition coefficient: secondary on secondary
  double c_PM = param["c_PM"]; // Competition coefficient: primary on secondary
  double c_MP = param["c_MP"]; // Competition coefficient: secondary on primary
  double mu_V = param["mu_V"]; //Vector mortality rate
  int ntime = param["ntime"] - 1; // Total simulation time
  int disturbance_time = param["disturbance_time"] - 1; // Disturbance time
  int infection_time = param["infection_time"] - 1;
  double delta_T = param["delta_T"]; // Time step
  double mortality_P = param["mortality_P"]; // Primary vector mortality
  double mortality_M = param["mortality_M"]; // Secondary vector mortality
  

  //Time anything
  for (int j = 0; j < ntime; j++) {
    
    // Total populations for rate calculations
    arma::rowvec NH = HS.row(j) + HI.row(j) + HR.row(j);
    arma::rowvec NP = PS.row(j) + PI.row(j);
    arma::rowvec NM = MS.row(j) + MI.row(j);
   
   if(j == infection_time){
     HI.row(j) = 1;
   }
   else{
    HI.row(j) = HI.row(j);
   }
  
    // Infection ratios
    arma::rowvec HP_ratio = PI.row(j) / NH;
    arma::rowvec HM_ratio = MI.row(j) / NH;
    arma::rowvec H_ratio =  HI.row(j) / NH;
    
    // Human infection rates
    arma::rowvec infections_H_fromP = theta_P * f_P * HP_ratio;
    arma::rowvec infections_H_fromM = theta_M * f_M * HM_ratio;
    arma::rowvec infections_H = infections_H_fromP + infections_H_fromM;
    
    // Vector infection rates
    arma::rowvec infections_P = theta_H * f_P * H_ratio;
    arma::rowvec infections_M = theta_H * f_M * H_ratio;
    
    // Human demographic/infection rates
    arma::rowvec HS_Rates = (b_H * NH) - (infections_H + mu_H) % HS.row(j);
    arma::rowvec HI_Rates = (infections_H % HS.row(j)) - (gamma + mu_H) * HI.row(j);
    arma::rowvec HR_Rates = (gamma * HI.row(j)) - (mu_H * HR.row(j));
    
    // Primary vector rates
    arma::rowvec PS_Rates = (b_P * NP) * (1-(c_MP * NM + c_PP * NP)) - (mu_V + infections_P) % PS.row(j);
    arma::rowvec PI_Rates = (infections_P % PS.row(j)) - (mu_V * PI.row(j));
    
    // Secondary vector rates
    arma::rowvec MS_Rates = (b_M * NM) *(1-(c_PM * NP + c_MM * NM)) - (mu_V + infections_M) % MS.row(j);
    arma::rowvec MI_Rates = (infections_M % MS.row(j)) - (mu_V * MI.row(j));
    
    // Update changes
    arma::rowvec total_change_HS = HS.row(j) + HS_Rates * delta_T;
    arma::rowvec total_change_HI = HI.row(j) + HI_Rates * delta_T;
    arma::rowvec total_change_HR = HR.row(j) + HR_Rates * delta_T;
    arma::rowvec total_change_PS = PS.row(j) + PS_Rates * delta_T;
    arma::rowvec total_change_PI = PI.row(j) + PI_Rates * delta_T;
    arma::rowvec total_change_MS = MS.row(j) + MS_Rates * delta_T;
    arma::rowvec total_change_MI = MI.row(j) + MI_Rates * delta_T;
    
    // Clamp values to avoid negatives
    total_change_HS = arma::clamp(total_change_HS, 0, 1e30);
    total_change_HI = arma::clamp(total_change_HI, 0, 1e30);
    total_change_HR = arma::clamp(total_change_HR, 0, 1e30);
    total_change_PS = arma::clamp(total_change_PS, 0, 1e30);
    total_change_PI = arma::clamp(total_change_PI, 0, 1e30);
    total_change_MS = arma::clamp(total_change_MS, 0, 1e30);
    total_change_MI = arma::clamp(total_change_MI, 0, 1e30);
    
    // Apply changes, with disturbance adjustments
    if (j != disturbance_time) {
      HS.row(j + 1) = total_change_HS;
      HI.row(j + 1) = total_change_HI;
      HR.row(j + 1) = total_change_HR;
      PS.row(j + 1) = total_change_PS;
      PI.row(j + 1) = total_change_PI;
      MS.row(j + 1) = total_change_MS;
      MI.row(j + 1) = total_change_MI;
      
    } else {
      
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
    Rcpp::Named("MI") = MI
  );
}
