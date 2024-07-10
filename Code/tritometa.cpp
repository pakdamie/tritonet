#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//[[Rcpp::export]]

List trito_metapop(
                double t,  
                NumericVector y, 
                int patch_num,
                NumericVector parms,
                NumericMatrix disp_mat) {
        
        //Demographic parameters
        NumericVector b_H = rep(parms["b_H"],patch_num); //Human birth rate
        NumericVector b_P = rep(parms["b_P"],patch_num); //P.vector birth rate
        NumericVector b_S = rep(parms["b_S"],patch_num);  //S. vector birth rate
        NumericVector mu_H = rep(parms["mu_H"],patch_num);  //Human death rate
        NumericVector mu_P = rep(parms["mu_P"],patch_num); //P. vector death rate
        NumericVector mu_S = rep(parms["mu_S"],patch_num); // S. vector death rate
        
        //Force of infection parameters
        NumericVector a_P = rep(parms["a_P"],patch_num); //biting rate of the p. vector
        NumericVector a_S =  rep(parms["a_S"],patch_num);  //biting rate of the s.vector
        
        NumericVector phi_P =   rep(parms["phi_P"],patch_num);  //transmission probability of p. vector
        NumericVector phi_S =    rep(parms["phi_S"],patch_num);   //transmission probability of s. vector
        NumericVector phi_H =   rep(parms["phi_H"],patch_num);   //transmission probability of human
        
        // Recovery rate
        NumericVector gamma = rep(parms["gamma"],patch_num);   //recovery rate of infected human
        
        //#competition coefficient
        NumericVector c_PS =  rep(parms["c_PS"],patch_num);   //competitition effect of p.vector on s.vector
        NumericVector c_SP =  rep(parms["c_SP"],patch_num);  //competition effect of s.vector on p.vector
        
        // FOI
        NumericVector FOI_P = a_P * phi_P; // #FOI for a primary vector
        NumericVector FOI_S = a_S * phi_S;  //#FOI for a secondary vector
        NumericVector FOI_H_P = a_P * phi_H; //#FOI for a human to a primary vector
        NumericVector FOI_H_S = a_S * phi_S; //#FOI for a human to a secondary vector
        
        // Human host population
        
        
        NumericVector HS = y[Range(0,patch_num-1)];
        NumericVector HI = y[Range(patch_num, (2*patch_num)-1)];
        NumericVector HR = y[Range((2*patch_num),((3*patch_num)-1))];
        
        NumericVector PS = y[Range((3*patch_num),((4*patch_num)-1))];
        NumericVector PI =  y[Range((4*patch_num),((5*patch_num)-1))];
        NumericVector SS =  y[Range((5*patch_num),((6*patch_num)-1))];
        NumericVector SI =  y[Range((6*patch_num),((7*patch_num)-1))];
        

        Rcpp::NumericVector dy(patch_num * 7);
        
        NumericVector dHS ;
        NumericVector dHI ;
        NumericVector dHR ;
        
        NumericVector dPS;
        NumericVector dPI ;
        NumericVector dSS ;
        NumericVector dSI ;
        
        

        // Population size
        NumericVector N_H = HS + HI + HR; //human host population 
        NumericVector N_P = PS + PI; // p. vector population
        NumericVector N_S = SS + SI; // s. vector population
                
                
        //human host population
        
        //susceptible host human population
        dHS = (b_H * N_H) - (FOI_P * HS * (PI/N_P)) - (FOI_S *HS*(SI/N_S)) - (mu_H *HS);
        
        //infected host human population
        dHI = (FOI_P*HS *(PI/N_P)) + (FOI_S *HS*(SI/N_S))- (gamma*HI) -(mu_H*HI);
        
        //recovered host human population
        dHR = (gamma*HI)-(mu_H*HR);
        
        //primary vector population
        
        //susceptible primary vector population
        dPS =  (b_P * N_P) - (FOI_H_P *PS* (HS/N_H)) - (mu_P *PS) - (c_SP*(PS)*(N_S)) +
                (transpose(disp_mat) * PS) -  (disp_mat * PS);
                
      
        //infected primary vector population
        dPI = (FOI_H_P * PS * (HS/N_H)) - mu_P*PI - c_SP*(PI)*(N_S) + 
                ((disp_mat) * PI) -  (disp_mat * PI);
                
        
        //secondary vector population
        dSS =  (b_S * N_S) - (FOI_H_S *SS* (HS/N_H)) - mu_S *SS -  c_PS*(SS)*(N_P)+
                (transpose(disp_mat) * SS) -  (disp_mat * SS);
                
        dSI = (FOI_H_S * SS * (HS/N_H)) - mu_S*SI -  c_PS*(SI)*(N_P) + 
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


/*** R
patch_num = 3
initial_y <- c(HS = rep(10,patch_num),
                 HI = rep(10,patch_num),
                 HR = rep(0,patch_num),
                 PS = rep(0,patch_num),
                 PI = rep(100,patch_num),
                 SS = rep(100,patch_num),
                 SI = rep(100,patch_num))



parameters_n <- c(
        b_H = 10, #Human birth rate
        b_P = 0, #P.vector birth rate
        b_S = 0, #S. vector birth rate
        mu_H = 0, #Human death rate
        mu_P = 0, #P. vector death rate
        mu_S = 0, #S. vector death rate
        
        a_P = 0.10, #biting rate of the p. vector
        a_S = 0.10, #biting rate of the s.vector
        
        phi_P = 0.90, #transmission probability of p. vector
        phi_S = 0.10, #transmission probability of s. vector
        phi_H  = 0.40, #transmission probability of human
        
        # Recovery rate
        gamma = 1/7,  #recovery rate of infected human
        
        #competition coefficient
        c_PS = 0, #competitition effect of p.vector on s.vector
        c_SP = 0  #competitition effect of s.vector on p.vector
)

df.contact <- matrix(c(0,0,0,
                       0,0,0,
                       0,0,0), ncol =3)

results <- deSolve::lsoda(
        times = 1:50,
        y = initial_y ,
        func = trito_metapop,
        parms = parameters_n,
        patch_num = 3,
        disp_mat =df.contact 
        
)




*/