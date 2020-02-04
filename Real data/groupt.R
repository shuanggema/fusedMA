f_group_t=function(Y,Z,p,N,Z_matrix_t,A_t,Lamda1,Lamda2,kappa,gamma,iter_admm,eps,penalty_index,Z_inv_t,ll1,ll2,Intial_theta_t,Lamda3,ll3){
res_theta=f_theta(Y,Z,p,N,Z_matrix_t,A_t,Lamda1,Lamda2,kappa,gamma,iter_admm,eps,penalty_index,Z_inv_t,ll1,ll2,Intial_theta_t,Lamda3,ll3)
hat_theta=array(NA,c(d2,N,ll1,ll2,ll3))
hat_alpha=array(NA,c(d2,N,ll1,ll2,ll3))
hat_size=array(NA,c(ll1,ll2,ll3))
bic=array(NA,c(ll1,ll2,ll3))
hat_label=array(NA,c(N,N,ll1,ll2,ll3))
gc()
for(lk in 1:ll3){
for(lj in 1:ll2){
  for(li in 1:ll1){
   theta=res_theta$Theta[,,li,lj,lk]
   theta[abs(theta)<=sparse]=0
   hat_y=diag(Z%*%theta)
   seq=1:N
   K=1
   Group_W=matrix(0,N,N)
   Label_W=matrix(0,N,N)
while(length(seq)>0){
     i=seq[1]
     ind=seq[-1]
     id=c(i)
     if(length(ind)>0){
       for(ii in 1:length(ind)){
       temp_id=f_trun(theta[,ind[ii]],theta[,i])
       if(temp_id==T){id=c(id,ind[ii])}
       }
     }
     Group_W[1:length(id),K]=id
     Label_W[id,K]=rep(1,length(id))
     seq=seq[!seq%in%id]
     K=K+1
  }
  K=K-1
  Label_W=Label_W[,1:K]
  if(K==1){W_label_len=f_label_len(Label_W)}
  if(K!=1){W_label_len=apply(Label_W,2,f_label_len)}
  hat_alpha_W=theta%*%Label_W/matrix(rep(W_label_len,d2),ncol=K,byrow=T)
  hat_alpha_W[which(is.na(hat_alpha_W[,1])==TRUE),]=0
  gc()

hat_theta[,,li,lj,lk]=theta
hat_alpha[,1:K,li,lj,lk]=hat_alpha_W
hat_size[li,lj,lk]=K
bic[li,lj,lk]=log(mean((Y-hat_y)^2))+log(N*(sum(hat_alpha_W!=0)))*log(N)*(sum(hat_alpha_W!=0))/N 
hat_label[,1:K,li,lj,lk]=Label_W
gc()
   
}
}
}   
res=list(hat_theta=hat_theta,hat_alpha=hat_alpha,hat_size=hat_size,bic=bic,hat_label=hat_label)
return(res)
}


