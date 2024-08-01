#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
//[[Rcpp::export]]

List trito_metapop(
                double t,  
                NumericVector y, 
                int patch_num,
                NumericVector parms,
                NumericMatrix disp_mat_P) {
        
        //Demographic parameters
        NumericVector b_H = rep(parms["b_H"],patch_num);    //Human birth rate
        NumericVector b_P = rep(parms["b_P"],patch_num);   //P.vector birth rate
        NumericVector b_S = rep(parms["b_S"],patch_num);   //S. vector birth rate
        NumericVector mu_H = rep(parms["mu_H"],patch_num); //Human death rate
        NumericVector mu_P = rep(parms["mu_P"],patch_num); //P. vector death rate
        NumericVector mu_S = rep(parms["mu_S"],patch_num); // S. vector death rate
        
        //Force of infection parameters
        NumericVector a_P = rep(parms["a_P"],patch_num); //biting rate of the p. vector
        NumericVector a_S = rep(parms["a_S"],patch_num);  //biting rate of the s.vector
        
        NumericVector phi_P =  rep(parms["phi_P"],patch_num);  //transmission probability of p. vector
        NumericVector phi_S =  rep(parms["phi_S"],patch_num);   //transmission probability of s. vector
        NumericVector phi_H =  rep(parms["phi_H"],patch_num);   //transmission probability of human
        
        // Recovery rate
        NumericVector gamma = rep(parms["gamma"],patch_num);   //recovery rate of infected human
        
        //#competition coefficient
        NumericVector c_PS =  rep(parms["c_PS"],patch_num);   //competition effect of p.vector on s.vector
        NumericVector c_SP =  rep(parms["c_SP"],patch_num);  //competition effect of s.vector on p.vector
        NumericVector disp_max = rep(parms["disp_max"],patch_num);
        
        // FOI
        NumericVector FOI_P = a_P * phi_P; // #FOI for a primary vector
        NumericVector FOI_S = a_S * phi_S;  //#FOI for a secondary vector
        NumericVector FOI_H_P = a_P * phi_H; //#FOI for a human to a primary vector
        NumericVector FOI_H_S = a_S * phi_H; //#FOI for a human to a secondary vector
        
        // Human host population
        NumericVector HS = y[Range(0,patch_num - 1)];
        NumericVector HI = y[Range(patch_num, (2 * patch_num) - 1)];
        NumericVector HR = y[Range((2 * patch_num),((3 * patch_num) - 1))];
        
        //Primary vectors
        NumericVector PS = y[Range((3 * patch_num),((4 * patch_num) - 1))];
        NumericVector PI = y[Range((4 * patch_num),((5 * patch_num) - 1))];
        
        //Secondary vectors
        NumericVector SS = y[Range((5 * patch_num),((6 * patch_num) - 1))];
        NumericVector SI = y[Range((6 * patch_num),((7 * patch_num) - 1))];
        

        Rcpp::NumericVector dy(patch_num * 7);
        
        //Human host population
        NumericVector dHS;
        NumericVector dHI;
        NumericVector dHR;
        
        //Vector population
        NumericVector dPS;
        NumericVector dPI;
        NumericVector dSS;
        NumericVector dSI;
        
        // Population size
        NumericVector N_H = HS + HI + HR; //human host population 
        NumericVector N_P = PS + PI; // p. vector population
        NumericVector N_S = SS + SI; // s. vector population
        
        
        //Disp_matrix adjusted for density-dependence
        arma::vec dispersal_vec = 1.0/((N_S + N_P * 1e10));
        
        
        //Distance_matrix
        
        
        
        arma::mat dispersal_mat(patch_num, patch_num);
        for (int p = 0; p < patch_num; p++){
                dispersal_mat (p, p) = dispersal_vec[p];
        }
        
        arma::mat disp_mat_P_convert = Rcpp::as<arma::mat>(disp_mat_P);
        arma::mat new_disp_matrix = disp_mat_P_convert * dispersal_mat;
        
        NumericMatrix disp_mat = Rcpp::wrap(new_disp_matrix);

   
        //human host population
        
        //susceptible host human population
        dHS = (b_H * N_H) - (FOI_P * HS * (PI/N_P)) - (FOI_S *HS*(SI/N_S)) - (mu_H *HS);
        
        //infected host human population
        dHI = (FOI_P * HS * (PI/N_P)) + (FOI_S * HS * (SI/N_S))- (gamma*HI) - (mu_H*HI);
        
        //recovered host human population
        dHR = (gamma * HI) - (mu_H * HR);
        
        //primary vector population
        
        //susceptible primary vector population
        dPS =  (b_P*N_P) - (FOI_H_P*PS*(HI/N_H)) - (mu_P*PS) - (c_SP*(PS)*(N_S)) +
                (transpose(disp_mat) * PS) -  (disp_mat * PS);
                
      
        //infected primary vector population
        dPI = (FOI_H_P * PS * (HI/N_H)) - (mu_P*PI) - (c_SP*(PI)*(N_S)) + 
                (transpose(disp_mat) * PI) -  (disp_mat * PI);
                
        
        //secondary vector population
        dSS =  (b_S * N_S) - (FOI_H_S *SS* (HI/N_H)) - (mu_S *SS) -  (c_PS*(SS)*(N_P)) +
                (transpose(disp_mat) * SS) -  (disp_mat * SS);
                
        dSI = (FOI_H_S * SS * (HI/N_H)) - (mu_S*SI) -  c_PS*(SI)*(N_P) + 
                (transpose(disp_mat) * SI) -  (disp_mat * SI);      
                
                
        dy[Range(0,patch_num-1)] = dHS ;
        dy[Range(patch_num, (2*patch_num)-1)] = dHI;
        dy[Range((2*patch_num),((3*patch_num)-1))] = dHR;
                
         dy[Range((3*patch_num),((4*patch_num)-1))] = dPS;
         dy[Range((4*patch_num),((5*patch_num)-1))] =    dPI ;
         dy[Range((5*patch_num),((6*patch_num)-1))] =  dSS ;
         dy[Range((6*patch_num),((7*patch_num)-1))] = dSI;
                
                
                

        return List::create(dy);
}

