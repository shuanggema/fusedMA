### ??machine ???sumodel 
f_group=function(M,Y,Z,Z_matrix,Z_inv,A,d2,n,B,temp_sam_z,temp_sam_y,Lamda1,Lamda2,kappa,gamma,iter_admm,eps,penalty_index,iter_fw,Intial_theta,Lamda3){

Hat_Y_B=array(0,c(N,B,ll1,ll2,ll3))
Theta_B=array(0,c(d2,N,B,ll1,ll2,ll3))
Sample_y_B=matrix(0,N,B)
Sample_z_B=matrix(0,d2,B)

hat_theta=array(NA,c(d2,N,ll1,ll2,ll3))
hat_alpha=array(NA,c(d2,N,ll1,ll2,ll3))
hat_size=array(NA,c(ll1,ll2,ll3))
bic=array(NA,c(ll1,ll2,ll3))
hat_label=array(NA,c(N,N,ll1,ll2,ll3))
hat_weight=array(NA,c(B,ll1,ll2,ll3))
gc()
if(M==1){CM=rep(M,B)}else{CM=rep(1:M,rep(M,B))}
CB=rep(1:B,M)
res_theta=foreach(m=CM,b=CB)%dopar%f_theta(Y[temp_sam_y[,m]],Z[temp_sam_y[,m],temp_sam_z[,b]],h,n,Z_matrix[,,b,m],A,Lamda1,Lamda2,kappa,gamma,iter_admm,eps,penalty_index,Z_inv[,,b,m],ll1,ll2,array(Intial_theta[,,b,m,],c(h,n,ll3)),Lamda3,ll3)
### ??machine ???sumodel 

for(m in 1:M){
   index_y=temp_sam_y[,m]
   Y_b=Y[index_y]
   
   for(b in 1:B){
   index_z=temp_sam_z[,b]
 #res_theta=f_theta(Y_b,Z[index_y,index_z],h,n,Z_matrix[,,b,m],A,Lamda1,Lamda2,kappa,gamma,iter_admm,eps,penalty_index,Z_inv[,,b,m],ll1,ll2,Intial_theta[,,b,m,],Lamda3,ll3)
  theta_b=res_theta[[b+B*(m-1)]]$Theta
  Theta_B[index_z,index_y,b,,,]=theta_b
  Sample_z_B[index_z,b]=rep(1,h)
  Sample_y_B[index_y,b]=rep(1,n)
  Hat_Y_B[index_y,b,,,]=res_theta[[b+B*(m-1)]]$Hat_y
  gc()
 }
}
gc()
### ???
for(lk in 1:ll3){
for(lj in 1:ll2){

for(li in 1:ll1){
      res_final=f_final(Y,Hat_Y_B[,,li,lj,lk],Theta_B[,,,li,lj,lk],Sample_z_B,Sample_y_B)
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
 gc()
res=list(hat_theta=hat_theta,hat_alpha=hat_alpha,hat_size=hat_size,bic=bic,hat_label=hat_label,hat_weight=hat_weight)
return(res)
}