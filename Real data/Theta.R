#### ADMM ??????lamda ????
f_theta=function(Y_b,Z_b,p_b,N_b,Z_matrix_b,A_b,Lamda1,Lamda2,kappa,gamma,iter_admm,eps,penalty_index,Z_inv_b,ll1,ll2,temp_intial_theta,Lamda3,ll3){
Temp=array(0,c(p_b,N_b,ll1,ll2,ll3))
hat_y=array(0,c(N_b,ll1,ll2,ll3))
for(lk in 1:ll3) {
for(lj in 1:ll2){
   temp_theta=temp_intial_theta[,,lk]
         for(li in 1:ll1) {
          #res_admm=f_admm(Y_b,Z_b,temp_theta,p_b,N_b,as.matrix(Z_matrix_b),as.matrix(A_b),Lamda1[li],Lamda2[lj],kappa,gamma,iter_admm,eps,as.matrix(Z_inv_b),penalty_index)
           res_admm=admmcpp(Y_b,Z_b,temp_theta,as.matrix(Z_matrix_b),as.matrix(A_b),Lamda1[li],Lamda2[lj],kappa,gamma,iter_admm,eps,as.matrix(Z_inv_b))
          #f_admm=function(Y_b,Z_b,the,len_the,n_m,Z_matrix,A,lamda1,lamda2,kappa,gamma,iter_admm,eps,Z_inv,penalty_index)
        # Y_b,Z_b,the,len_the,n_m,Z_matrix,A,lamda1,lamda2,kappa,gamma,iter_admm,eps,Z_inv
            temp_theta=res_admm$theta_b
        Temp[,,li,lj,lk]=matrix(temp_theta,nrow=p_b,byrow=F)
        hat_y[,li,lj,lk]=as.numeric(Z_matrix_b%*%temp_theta)
         }
}
}
res=list(Theta=Temp,Hat_y=hat_y)
return(res)
}