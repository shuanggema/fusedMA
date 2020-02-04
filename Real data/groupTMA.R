f_group_T=function(M_T,n_T,Y,Z,Z_matrix_T,Z_inv_T,A_T,d2,h,N,B,temp_sam_z,temp_theta,Lamda1,Lamda2,kappa,gamma,iter_admm,eps,penalty_index,iter_fw,Intial_theta_T){

Hat_Y_B=array(0,c(N,B,ll1,ll2,ll3))
Theta_B=array(0,c(d2,N,B,ll1,ll2,ll3))
Sample_z_B=matrix(0,d2,B)
Sample_y_B=matrix(0,N,B)
hat_theta=array(NA,c(d2,N,ll1,ll2,ll3))
hat_alpha=array(NA,c(d2,N,ll1,ll2,ll3))
hat_size=array(NA,c(ll1,ll2,ll3))
bic=array(NA,c(ll1,ll2,ll3))
hat_label=array(NA,c(N,N,ll1,ll2,ll3))
hat_weight=array(NA,c(B,ll1,ll2,ll3))
gc()
CB=1:B
res_theta=foreach(b=CB)%dopar%f_theta(Y,Z[,temp_sam_z[,b]],h,N,Z_matrix_T[,,b],A_T,Lamda1,Lamda2,kappa,gamma,iter_admm,eps,penalty_index,Z_inv_T[,,b],ll1,ll2,array(Intial_theta_T[,,b,],c(h,N,ll3)),Lamda3,ll3)
for(b in 1:B){
   index_z=temp_sam_z[,b]
  # res_theta=f_theta(Y,Z[,index_z],h,N,Z_matrix_T[,,b],A_T,Lamda1,Lamda2,kappa,gamma,iter_admm,eps,penalty_index,Z_inv_T[,,b],ll1,ll2,Intial_theta_T[,,b])
  theta_b=res_theta[[b]]$Theta
  Theta_B[index_z,,b,,,]=theta_b
  Sample_z_B[index_z,b]=rep(1,h)
  Sample_y_B[,b]=rep(1,N)
  Hat_Y_B[,b,,,]=res_theta[[b]]$Hat_y
  gc()
 }


 for(lk in 1:ll3){
for(lj in 1:ll2){

for(li in 1:ll1){
      res_final=f_final(M_T,n_T,Y,Hat_Y_B[,,li,lj,lk],Theta_B[,,,li,lj,lk],Sample_z_B,Sample_y_B)
      gc()
hat_theta[,,li,lj,lk]=res_final$hat_theta
hat_alpha[,1:res_final$hat_size,li,lj,lk]=res_final$hat_alpha
hat_size[li,lj,lk]=res_final$hat_size
bic[li,lj,lk]=res_final$bic
hat_label[,1:res_final$hat_size,li,lj,lk]=res_final$hat_label
hat_weight[,li,lj,lk]=res_final$hat_weight

}
}
}
res=list(hat_theta=hat_theta,hat_alpha=hat_alpha,hat_size=hat_size,bic=bic,hat_label=hat_label,hat_weight=hat_weight)
return(res)
gc()
}