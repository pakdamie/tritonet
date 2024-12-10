#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat calculate_RE_onepatch(
    Rcpp::List x, 
    Rcpp::NumericVector param, 
    int secondInclude){
  
  double theta_H = param["theta_H"]; //transmission probability of human
  double theta_P = param["theta_P"]; //transmission probbility of p.vector
  double theta_M = param["theta_M"]; //transmission probbility of s.vector
  double f_P = param["f_P"]; // Biting rate of p.vector
  double f_M = param["f_M"];  // Biting rate of s.vector
  double gamma = param["gamma"];  // Recovery rate of human
  double mu_H = param["mu_H"]; // Human mortality rate
  double c_MP = param["c_MP"]; //Competition of secondary on primary
  double c_PM = param["c_PM"]; //Competition of primary on secondary
  double c_MM = param["c_MM"];  //Competition of secondary on secondary
  double c_PP = param["c_PP"]; //Competition of primary on primary
  double ntime = param["ntime"]; //Length of simulation
        
  arma::mat HS = x[0]; //Susceptible humans
  arma::mat HI = x[1]; //Infected humans
  arma::mat HR = x[2];  //Recovered humans
  
  arma::mat PS = x[3]; //Susceptible Primary
  arma::mat PI = x[4]; //Infectious Primary
  
  arma::mat MS = x[5]; //Susceptible Secondary
  arma::mat MI = x[6]; //Infectious Secondary
  
  arma::mat NH = HS + HI + HR; //Total humans
  arma::mat NP = PS + PI; //Total primary
  arma::mat NM = MS + MI; //Total secondary
  
  // This is a matrix that we fill with the RE and the current
  // populations at time j (row)
  
  arma::mat mat_RE = arma::mat(NM.n_rows, 10, arma::fill::zeros);
        
  //I don't need to really run it that long
  for (int j = 0; j <= ntime - 100; ++j){

    //primary or secondary infecting humans
    double p_inf_h = theta_P * f_P * HS[j]/NH[j];
    double m_inf_h = theta_M * f_M * HS[j]/NH[j];
  
    //human infection primary or secondary
    double h_inf_p = theta_H * f_P * PS[j]/NH[j];
    double h_inf_m = theta_H * f_M * MS[j]/NH[j];
    
    //the F-matrix     
    arma::Mat F_mat = arma::mat(3,3,arma::fill::zeros) ;
      
    //Remember C++ indexing starts at 0 instead of 1
    F_mat(0, 1) =  p_inf_h; //in R (1, 2)
    F_mat(0, 2) =  m_inf_h; //in R (1, 3)
    F_mat(1, 0) =  h_inf_p; //in R(2, 1)
    F_mat(2, 0) =  h_inf_m; //in R(3, 1)
       
    //the V-matrix 
    arma::mat V_mat = arma::mat(3,3,arma::fill::zeros) ;
       
    double h_out = gamma + mu_H;
    double p_out_diag = c_PP * (PS[j] + (2) * PI[j]) + (c_MP *(NM[j]));
    double p_out_nondiag = c_MP * PI[j];
    double m_out_diag = c_MM * (MS[j] + (2) * MI[j]) + (c_PM *(NP[j]));
    double m_out_nondiag = c_PM * MI[j];

    V_mat(0,0) = h_out; // In R (1,1)
    V_mat(1,1) = p_out_diag ; // In R (2,2)
    V_mat(2,2) = m_out_diag; // In R (3,3)
    V_mat(1,2) = p_out_nondiag; // In R (2,3)
    V_mat(2,1) = m_out_nondiag; // In R (3,2)
       
    if (secondInclude == 0){
      arma::mat F_small_mat = F_mat.submat(0, 0, 1, 1);
      arma::mat V_small_mat = V_mat.submat(0, 0, 1, 1);
      arma::mat V_small_inv_mat = inv(V_small_mat);
      arma::mat FV_small_mat = F_small_mat * V_small_inv_mat;
      arma::cx_vec eigval = arma::eig_gen(FV_small_mat);
      mat_RE(j,0) =  arma::abs(eigval).max();
      mat_RE(j,1) = NP[j];
      mat_RE(j,2) = NM[j];
      mat_RE(j,3) = HS[j];
      mat_RE(j,4) = PS[j];
      mat_RE(j,5) = MS[j];
      mat_RE(j,6) = HI[j];
      mat_RE(j,7) = PI[j];
      mat_RE(j,8) = MI[j];
      mat_RE(j,9) = j;
               
     } else {
      arma::mat V_inv_mat = inv(V_mat);    
      arma::mat FV_mat = F_mat * V_inv_mat;
      arma::cx_vec eigval = arma::eig_gen(FV_mat);
      mat_RE(j,0) =  arma::abs(eigval).max();
      mat_RE(j,1) = NP[j];
      mat_RE(j,2) = NM[j];
      mat_RE(j,3) = HS[j];
      mat_RE(j,4) = PS[j];
      mat_RE(j,5) = MS[j];
      mat_RE(j,6) = (theta_H * f_P * PS[j]/NH[j] + 
                    theta_H * f_M * MS[j]/NH[j]) * (HI[j]); 
      mat_RE(j,7) = p_inf_h * (PI[j]);
      mat_RE(j,8) = m_inf_h * (MI[j]);
      mat_RE(j,9) = j;
      }
  }
      
        return(mat_RE);
        
}
