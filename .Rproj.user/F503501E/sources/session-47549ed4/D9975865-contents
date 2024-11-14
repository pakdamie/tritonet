Calculate_R_effective (adjacency_matrix, inf_states, other_states,
                       F_matrix_test, V_matrix_test, desolve_init,
                       tryit)

NP <-  (calculate_vector_abundance_patch(tmp)[[1]] + calculate_vector_abundance_patch(tmp)[[2]])[,2:87]
colnames(NP) <- paste("NP", 1:86, sep="")

NS <-  (calculate_vector_abundance_patch(tmp)[[3]] + calculate_vector_abundance_patch(tmp)[[4]])[,2:87]
colnames(NS) <- paste("NS", 1:86, sep="")

NV <- NP  + NS
colnames(NV) <- paste("NV", 1:86, sep="")

NH<- (calculate_host_abundance_patch(tmp)[[1]] + (calculate_host_abundance_patch(tmp))[[2]] + 
         (calculate_host_abundance_patch(tmp))[[3]])[,2:87]
colnames(NH) <- paste("NH", 1:86, sep="")

                                                
tryit<- cbind(tmp, NP, NS, NV, NH)

