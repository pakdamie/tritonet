#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List discrete_trito_model_rcpp(
    arma::mat HS, //Susceptible human matrix
    arma::mat HI, //Infectious human matrix
    arma::mat HR, //Recovered human matrix
    arma::mat PS, //Susceptible p.vector matrix
    arma::mat PI, // Infectious p.vector matrix
    arma::mat MS, // Susceptible s.vector matrix
    arma::mat MI, // Infectious s.vector matrix
    arma::mat adj, // Adjacency matrix
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
  double d = param["d"]; //Dispersal rate
  int ntime = param["ntime"] - 1; //How long to simulate for (stupid indexing)
  int disturbance_time = param["disturbance_time"] - 1; //When is disturbance time?
  double delta_T = param["delta_T"]; //Time steps
  double prop = param["prop"]; //Proportion of patches to be affected
  double mortality_P = param["mortality_P"]; //Surviviorship of p.vector
  double mortality_M = param["mortality_M"]; //Survivorship of s.vector
  
  //The number of patches is dependent on either the column or row (should 
  //be equal)
  int patch_num = adj.n_rows;
  
  double sample_size = floor(prop * patch_num); // Floor the coverage
  
  // Generate the sequence from 1 to patch_num
  arma::vec patches = arma::linspace<arma::vec>(1, patch_num,patch_num)-1;            
  // Sample without replacement
  
  arma::vec sampled_coverage = Rcpp::RcppArmadillo::sample(patches, sample_size, false);

  
  
  
  for(int j = 0; j < ntime; j++){
    
    //Total populations of human, primary, and secondary
    arma::rowvec NH = HS.row(j) + HI.row(j) + HR.row(j); 
    arma::rowvec NP = PS.row(j) + PI.row(j);
    arma::rowvec NM = MS.row(j) + MI.row(j);
    
    //Ratios of vector to human and infected human to total human
    arma::rowvec HP_ratio = PI.row(j)/NH;
    arma::rowvec HM_ratio = MI.row(j)/NH;
    arma::rowvec H_ratio =  HI.row(j)/NH;
    
    //Primary and secondary will infect humans
    arma::rowvec infections_H = (theta_P * f_P * HP_ratio) + (theta_M * f_M * HM_ratio);
    
    //Primary and secondary will be infected by human
    arma::rowvec infections_P = (theta_H * f_P) *  H_ratio;
    arma::rowvec infections_M = (theta_H * f_M) *  H_ratio;
    
    //Demographic/infection rates of humans
    arma::rowvec HS_Rates = (b_H * NH) - (infections_H + mu_H) % HS.row(j);
    arma::rowvec HI_Rates = (infections_H % HS.row(j)) - (gamma + mu_H) * HI.row(j);
    arma::rowvec HR_Rates = (gamma * HI.row(j)) - (mu_H * HR.row(j));
    
    //Demographic/infection rates of primary vector
    arma::rowvec PS_Rates = (b_P * NP) - (infections_P + 
      (c_MP * NM) + (c_PP * NP)) % PS.row(j); 
    
    arma::rowvec PI_Rates = (infections_P % PS.row(j)) - (
      (c_MP * NM) + (c_PP * NP)) % PI.row(j); 
    
    //Demographic/infection of secondary vector
    arma::rowvec  MS_Rates = (b_M * NM) - (infections_M +
      (c_PM * NP) + (c_MM * NM)) % MS.row(j);
    
    arma::rowvec  MI_Rates = (infections_M % MS.row(j)) - (
      (c_PM * NP) + (c_MM * NM)) % MI.row(j);
    
    
    // Dispersal of the primary and secondary vectors
    arma::vec dispersal_PS = arma::sum((adj.each_row() % (-d * PS.row(j))),1)
      + adj.t() *  (d * PS.row(j).t());
    
    arma::vec  dispersal_PI = arma::sum((adj.each_row() % (-d * PI.row(j))),1) +
    adj.t() *  (d * PI.row(j).t());
    
    arma::vec  dispersal_MS = arma::sum((adj.each_row() % (-d * MS.row(j))),1) +
      adj.t() *  (d * MS.row(j).t());
    
    arma::vec  dispersal_MI = arma::sum((adj.each_row() % (-d * MI.row(j))),1) +
      adj.t() *  (d * MI.row(j).t());
    
    //Let's look what the new population size will be
    arma::rowvec total_change_HS = HS.row(j) + (HS_Rates *  delta_T);
    arma::rowvec total_change_HI = HI.row(j) + (HI_Rates *  delta_T);
    arma::rowvec total_change_HR = HR.row(j) + (HR_Rates *  delta_T);
    arma::rowvec total_change_PS = PS.row(j) + (PS_Rates + dispersal_PS.t()) * delta_T;
    arma::rowvec total_change_PI = PI.row(j) + (PI_Rates + dispersal_PI.t()) *  delta_T;
    arma::rowvec total_change_MS = MS.row(j) + (MS_Rates + dispersal_MS.t()) *  delta_T;
    arma::rowvec total_change_MI = MI.row(j) + (MI_Rates + dispersal_MI.t()) *  delta_T;
    
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
      arma::uvec sampled_index = arma::conv_to<arma::uvec>::from(sampled_coverage);
      
      HS.row(j + 1) = total_change_HS;
      HI.row(j + 1) = total_change_HI;
      HR.row(j + 1) = total_change_HR;
      
      total_change_PS.cols(sampled_index) = total_change_PS.cols(sampled_index) * mortality_P;
      total_change_PI.cols(sampled_index) = total_change_PI.cols(sampled_index) * mortality_P;
      total_change_MS.cols(sampled_index) = total_change_MS.cols(sampled_index) * mortality_M;
      total_change_MI.cols(sampled_index) = total_change_MI.cols(sampled_index) * mortality_M;
      
      
      PS.row(j + 1) = total_change_PS;
      PI.row(j + 1) = total_change_PI;
      MS.row(j + 1) = total_change_MS;
      MI.row(j + 1) = total_change_MI;
      
    }
  }
  
  
  
  return Rcpp::List::create(Rcpp::Named("HS") = HS,
                            Rcpp::Named("HI") = HI,
                            Rcpp::Named("HR") = HR,
                            Rcpp::Named("PS") = PS,
                            Rcpp::Named("PI") = PI,
                            Rcpp::Named("MS") = MS,
                            Rcpp::Named("MI") = MI,
                            Rcpp::Named("Target_Patch") = sampled_coverage);
  
  
}

