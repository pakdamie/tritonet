#' Calculate the analytical R_effective
#'
#' @param param 
#'
#' @return
#' @export
#'
#' @examples
Calculate_analytical_REff <- function(List_x, param){
   
   theta_H = param["theta_H"] 
   theta_P = param["theta_P"] 
   theta_M = param["theta_M"] 
   f_P = param["f_P"]
   f_M =param["f_M"]
   gamma = param["gamma"]
   mu_H = param["mu_H"]
   c_MP = param["c_MP"]
   c_PM = param["c_PM"]
   c_MM = param["c_MM"]
   c_PP = param["c_PP"]
   ntime = param["ntime"]
    
   HS = List_x[[1]]
   HI = List_x[[2]]
   HR = List_x[[3]]
   
   PS = List_x[[4]]
   PI = List_x[[5]]
   
   MS = List_x[[6]]
   MI = List_x[[7]]
   
   NH = HS + HI + HR
   NP = PS + PI
   NM = MS + MI
       
   Psi_P =  (c_MP * NP) + c_MM * ( MS + (2 * MI))
   Psi_M =  (c_PM * NM) + c_PP * (PS + (2 * PI))
   Psi_C =  (c_MP * PI) * (c_PM * MI)
   wait_time = sqrt(((Psi_P * Psi_M) - Psi_C) * (gamma + mu_H)) * NH
   
   HtoH_P = sqrt(Psi_M * ((theta_H * f_P * PS) * (theta_P * f_P * HS)))/wait_time
   HtoH_M = sqrt(Psi_P * ((theta_H * f_M * MS) * (theta_M * f_M * HS)))/wait_time
   HtoH = HtoH_P + HtoH_M 
   
   MtoP = sqrt((c_PM * MI) * ((theta_H * f_P * PS ) + (theta_M * f_M * HS)))/wait_time
   PtoM = sqrt((c_MP * PI) * ((theta_H * f_M * MS) + (theta_P * f_P * HS)))/wait_time
   
   VtoV = -MtoP - PtoM
   
   #This is the total effective reproductive number
   RE = HtoH + VtoV
   
   RE_DF <- cbind.data.frame(time = 
                    seq(1,ntime),
                    RE = RE,
                    HtoH_P, HtoH_M, 
                    MtoP, PtoM)
        
   return(RE_DF)
}