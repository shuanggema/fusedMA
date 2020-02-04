#### ?????????????theta ?????????group
f_final=function(M,n,Y,Hat_Y_B,Theta_B,Sample_z_B,Sample_y_B)
{
 ### ???????hat_Y_B, ??weight
 res_w=f_weight(Y,Hat_Y_B,iter_fw,eps,B)
 hat_weight=res_w$w
 Matrix_w1=matrix(rep(hat_weight,d2),ncol=B,byrow=T)
 WTheta=matrix(0,nrow=d2,ncol=N)
 W_hat_y=rep(0,N)
###  ??weight????Theta_B????????
 for(j in 1:N){  
    ## theta
    Temp_TW=apply(Theta_B[,j,]*Matrix_w1*Sample_z_B,1,sum)
    #Temp_WS=apply(Matrix_w1*Sample_z_B,1,sum) 
    #WTheta[,j]=Temp_TW/Temp_WS
    WTheta[,j]=Temp_TW
    ## y
    #W_hat_y[j]=sum(Hat_Y_B[j,]*res_w$w*Sample_y_B[j,])/sum(res_w$w*Sample_y_B[j,])
#    W_hat_y[j]=sum(Hat_Y_B[j,]*res_w$w*Sample_y_B[j,])
  } 
WTheta[which(is.na(WTheta[,1])==TRUE),]=0
WTheta[abs(WTheta)<=sparse]=0
W_hat_y=diag(Z%*%WTheta)
seq=1:N
K=1
Group_W=matrix(0,N,N)
Label_W=matrix(0,N,N)
####????
while(length(seq)>0){
   i=seq[1]
   ind=seq[-1]
   id=c(i)
   if(length(ind)>0){
     for(ii in 1:length(ind)){
     temp_id=f_trun(WTheta[,ind[ii]],WTheta[,i])
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

W_hat_y[which(is.na(W_hat_y)==TRUE)]=0

  
if(K==1){W_label_len=f_label_len(Label_W)}
if(K!=1){W_label_len=apply(Label_W,2,f_label_len)}
  
hat_alpha=WTheta%*%Label_W/matrix(rep(W_label_len,d2),ncol=K,byrow=T)
hat_alpha[which(is.na(hat_alpha[,1])==TRUE),]=0
bic=rep(0,M)
for(i in 1:M){
bic[i]=log(mean((Y[n*(i-1)+1:n]-W_hat_y[n*(i-1)+1:n])^2))+log(n*(sum(hat_alpha!=0)))*log(n)*(sum(hat_alpha!=0))/n}
bic=sum(bic)



gc()
res=list(hat_weight=hat_weight,hat_size=K,hat_theta=WTheta,hat_label=Label_W,hat_alpha=hat_alpha,hat_Y=W_hat_y,bic=bic)
return(res)
 }
