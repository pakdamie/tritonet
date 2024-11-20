if (secondInclude == 0){
    arma::mat F_small_mat = F_mat.submat(0, 0, 1, 1);
    arma::mat V_small_mat = V_mat.submat(0, 0, 1, 1);
    arma::mat V_small_inv_mat = inv(V_small_mat);
    arma::mat FV_small_mat = F_small_mat * V_small_inv_mat;
    arma::cx_vec eigval = arma::eig_gen(FV_small_mat);
    mat_RE(j) =  arma::abs(eigval).max();
}else{
    arma::mat V_inv_mat = inv(V_mat);    
    arma::mat FV_mat = F_mat * V_inv_mat;
    arma::cx_vec eigval = arma::eig_gen(FV_mat);
    mat_RE(j) =  arma::abs(eigval).max();

NumericVector calculate_RE_Patch_cpp(List x, NumericVector param,
                                     int disturbance_time) {

    double theta_P = param["theta_P"];
    double theta_S = param["theta_S"];
    double theta_H = param["theta_H"];
    double gamma = parameters["gamma"];
    double mu_H = parameters["mu_H"];
    double f_P = parameters["f_P"];
    double f_S = parameters["f_S"];
    double c_SP = parameters["c_SP"];
    double c_PS = parameters["c_PS"];
    double c_SS = parameters["c_SS"];
    double c_PP = parameters["c_PP"];
                
    arma::mat HS = x[0];
    arma::mat HI = x[1];
    arma::mat HR = x[2];
    
    arma::mat PS = x[3];
    arma::mat PI = x[4];
    
    arma::mat SS = x[5];
    arma::mat SI = x[6];
    
    arma::mat NH = HS + HI + HR;
    arma::mat NP = PS + PI;
    arma::mat NS = SS + SI;
                
                
    IntegerVector interest <- seq(disturbance_time - 10, disturbance_time + 100);
                
    List <- NULL; 
    
    for(int t = 0; t < size(interest); t++){
        
        arma::rowvec HS_int = HS.row(t);
        arma::rowvec HI_int = HI.row(t);
        arma::rowvec HR_int = HR.row(t);
        
        arma::rowvec PS_int = PS.row(t);
        arma::rowvec PI_int = PI.row(t);

        arma::rowvec SS_int = SS.row(t);
        arma::rowvec SI_int = SI.row(t);
        
        arma::rowvec NH_int = NH.row(t);
        arma::rowvec NP_int = NP.row(t);
        arma::rowvec NS_int = NS.row(t);
        
    for (int p = 0; p < patch_num; p++)){
     
     PH = theta_P * f_P * (HS_int[k]/NH_int[k]);
     PS = theta_S * f_S * (HS_int[k]/NH_int[k]);
     HP = theta_P * f_P * (PS_int[k]/NH_int[k]);
     HS = theta_S * f_S * (SS_int[k]/NH_int[k]);
     
     arma::mat F_mat = Matrix(c(
                     0, PH, PS,
                     HP, 0, 0,
                     HS, 0, 0), 
                     byrow = TRUE, ncol = 3,sparse = TRUE); 
     
     H_rate = 1/(gamma + mu_H); 
     P_rate = 1/((c_SP * NS_int[k]) + (c_PP * NP_int[k]))
     S_rate = 1/((c_PS * NP_int[k]) + (c_SS * NS_int[k]))
     
     arma::mat V_mat <- Matrix(H_rate 0, 0,
                           0, P_rate, 0,
                           0, 0, S_rate,
                           byrow = TRUE, ncol = 3,sparse = TRUE)
         
         
   arma::mat FV <- F_mat %*% V_mat;
         
   patch_R0_time[[k]] <- cbind(RE = max(eigen(FV)$values),
                                     
                                     patch_num = k,
                                     time = interest_time,
                                     primary_contrib = sum(FV[1,2]),
                                     secondary_contribut = sum(FV[1,3]))
                        }
                        full_time[[t]] <- do.call(rbind, patch_R0_time)
                                
                                
                }
                
                
                full_time_f <- do.call(rbind.data.frame, full_time)
                        maxRE <- full_time_f[which.max(full_time_f$RE),]
                
                max_CV <- max(by(full_time_f, full_time_f$time, function(x)  CV_RE = sd(x$RE)/mean(x$RE)))
                        maxRE$CV <- max_CV
                
                
                return(list(full_time_f ,maxRE))
                        
        }
        
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
