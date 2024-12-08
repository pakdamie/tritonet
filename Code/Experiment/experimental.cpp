#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List discrete_trito_model(
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
  double theta_H = param["theta_H"]; // Transmission prob of human
  double theta_P = param["theta_P"]; // Transmission prob of p.vector
  double theta_S = param["theta_S"]; // Transmission prob of s.vector
  double f_P = param["f_P"]; //Biterate of p.vector
  double f_S = param["f_S"];  //Bite rate of s.vector
  double gamma = param["gamma"]; //Recovery rate for human
  double c_PP = param["c_PP"]; //Competition coefficient for primary on primary
  double c_SS = param["c_SS"]; //Competition coefficient for secondary on secondary
  double c_PS = param["c_PS"];  //Competition coefficient for primary on secondary
  double c_SP = param["c_SP"]; //Competition coefficient for secondary on primary
  double d = param["d"]; //Dispersal rate
  int ntime = param["ntime"]; //How long to simulate for
  int disturbance_time = param["disturbance_time"]; //When is disturbance time?
  double delta_T = param["delta_T"]; //Time steps
  double prop = param["prop"]; //Proportion of patches to be affected
  double mortality_P = param["mortality_P"]; //Surviviorship of p.vector
  double  mortality_S = param["mortality_S"]; //Survivorship of s.vector
        
        //The number of patches is dependent on either the column or row (should 
        //be equal)
        int patch_num = adj.n_rows;
       
       double j = 0;
                
                //Total populations of human, primary, and secondary
                arma::rowvec NH =  HS.row(j) + HI.row(j) + HR.row(j); 
                arma::rowvec NP = PS.row(j) + PI.row(j);
                arma::rowvec NS = SS.row(j) + SI.row(j);
                
                //Ratios of vector to human and infected human to total human
                arma::rowvec HP_ratio = PI.row(j)/NH;
                arma::rowvec HS_ratio = SI.row(j)/NH;
                arma::rowvec H_ratio =  HI.row(j)/NH;
                
                //Primary and secondary will infect humans
                arma::rowvec infections_H = (theta_P * f_P * HP_ratio) + (theta_S * f_S * HS_ratio);
                
                //Primary and secondary will be infected by human
                arma::rowvec infections_P = (theta_H * f_P) *  H_ratio ;
                arma::rowvec infections_S = (theta_H * f_S) *  H_ratio ;
                
        
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
        
        arma::rowvec  SI_Rates = (infections_S % SS.row(j)) - ((c_PS * NP) +
                (c_SS * NS)) % SI.row(j);
        
        
        // Dispersal of the primary and secondary vectors
        arma::vec dispersal_PS = arma::sum((adj.each_row() % (-d * PS.row(j))),1)
                + adj.t() *  (d * PS.row(j).t());
        
        arma::vec  dispersal_PI = arma::sum((adj.each_row() % (-d * PI.row(j))),1) +
        adj.t() *  (d * PI.row(j).t());
        
        arma::vec  dispersal_SS = arma::sum((adj.each_row() % (-d * SS.row(j))),1) +
                adj.t() *  (d * SS.row(j).t());
        
        arma::vec  dispersal_SI = arma::sum((adj.each_row() % (-d * SI.row(j))),1) +
                adj.t() *  (d * SI.row(j).t());
        
        arma::rowvec total_change_HS = HS.row(j) + (HS_Rates *  delta_T);
        arma::rowvec total_change_HI = HI.row(j) + (HI_Rates *  delta_T);
        arma::rowvec total_change_HR = HR.row(j) + (HR_Rates *  delta_T);
        arma::rowvec total_change_PS = PS.row(j) + (PS_Rates + dispersal_PS.t()) * delta_T;
        arma::rowvec total_change_PI = PI.row(j) + (PI_Rates + dispersal_PI.t()) *  delta_T;
        arma::rowvec total_change_SS = SS.row(j) + (SS_Rates + dispersal_SS.t()) *  delta_T;
        arma::rowvec total_change_SI = SI.row(j) + (SI_Rates + dispersal_SI.t()) *  delta_T;
        
        
        total_change_HS = arma::clamp(total_change_HS, 0, total_change_HS.max());
        total_change_HI = arma::clamp(total_change_HI, 0, total_change_HI.max());
        total_change_HR = arma::clamp(total_change_HR, 0, total_change_HR.max());
        total_change_PS = arma::clamp(total_change_PS, 0, total_change_PS.max());
        total_change_PI = arma::clamp(total_change_PI, 0, total_change_PI.max());
        total_change_SS = arma::clamp(total_change_SS, 0, total_change_SS.max());
        total_change_SI = arma::clamp(total_change_SI, 0, total_change_SI.max());
        
        
        if (j != disturbance_time) {
                
                // We then update these row by row
                HS.row(j + 1) = total_change_HS ;
                HI.row(j + 1) =  total_change_HI;
                HR.row(j + 1) =  total_change_HR;
                PS.row(j + 1) =  total_change_PS ;
                PI.row(j + 1) =  total_change_PI ;
                SS.row(j + 1) =  total_change_SS ;
                SI.row(j + 1) =  total_change_SI ;
        }
        
        else if (j == disturbance_time) {
                
                double sample_size = floor(prop * patch_num); // Floor the coverage
                
                // Generate the sequence from 1 to patch_num
                arma::vec patches = arma::linspace<arma::vec>(1, patch_num,patch_num)-1;            
                // Sample without replacement
                
                arma::vec sampled_coverage = Rcpp::RcppArmadillo::sample(patches, sample_size, false);
                arma::uvec sampled_index= arma::conv_to<arma::uvec>::from(sampled_coverage);
                
                
                HS.row(j + 1) = total_change_HS;
                HI.row(j + 1) = total_change_HI;
                HR.row(j + 1) = total_change_HR;
                
                total_change_PS.cols(sampled_index) = total_change_PS.cols(sampled_index) * mortality_P;
                total_change_PI.cols(sampled_index) = total_change_PI.cols(sampled_index) * mortality_P;
                total_change_SS.cols(sampled_index) = total_change_SS.cols(sampled_index) * mortality_S;
                total_change_SI.cols(sampled_index) = total_change_SI.cols(sampled_index) * mortality_S;
                
                
                PS.row(j + 1) = total_change_PS;
                PI.row(j + 1) = total_change_PI;
                SS.row(j + 1) = total_change_SS;
                SI.row(j + 1) = total_change_SI;
                
        
}
        
        return Rcpp::List::create(Rcpp::Named("HS") = NH,
                                  Rcpp::Named("HI") = NP,
                                  Rcpp::Named("HR") = NS,
                                  Rcpp::Named("PS") = HP_ratio,
                                  Rcpp::Named("PI") = HS_ratio,
                                  Rcpp::Named("SS") = HS_Rates);
        
        
}





/*** R
HS = matrix(1,nrow = 3,ncol =3)
HI = matrix(1,nrow = 3,ncol =3)
HR = matrix(1,nrow = 3,ncol =3)
PS = matrix(1,nrow = 3,ncol =3)
PI = matrix(1,nrow = 3,ncol =3)
SS = matrix(1,nrow = 3,ncol =3)
SI = matrix(1,nrow = 3,ncol =3)
adj = matrix(c(0,1,1,1,0,1,0,1,1),nrow = 3, ncol =3)

param_ALT<- c(b_H =1,
              b_P = 1/100,
              b_S = 1/100,
              mu_H = 1/100,
              f_P = 0.050, # Biting rate of the p. vector
              f_S = 0.030, # Biting rate of the s.vector
              theta_P = 0.50, # Transmission probability of p. vector
              theta_S = 0.30, theta_H = 0.50,
              gamma = 1/100, c_PP = 0, c_SS = 0, c_PS = 0,c_SP=0,d = 0.01,
              ntime = 1000, disturbance_time=500, delta_T=0.1,prop=0.5, mortality_P=0.8,
              mortality_S = 0.5) # Transmission probability of s. vector)



discrete_trito_model(HS,HI ,
                     HR ,PS ,PI , SS ,SI ,adj,param_ALT)    

*/
