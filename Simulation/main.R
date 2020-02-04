rm(list=ls(all=TRUE))
#setwd("D:/subgroup/R programming/20190124")
#Rcpp::sourceCpp("admmcpp.cpp")
set.seed(99)
nCores <- Sys.getenv("TOTALPROCS")
print(nCores)
library(foreach)
library(doMC)
registerDoMC(nCores)
library(doRNG)
registerDoRNG(1)
Nsim=50
### the total sample size
N=100
### the total dimension
p=100
### the sparsity
s=4
### size of group
KK=2
### the probabilty of different group
pro=rep(1/KK,KK)
## the dimension of subject-specific covariates
d2=100
alpha0=2
#T_alpha=runif(s,1.8,2)
T_alpha=rep(alpha0,s)
theta0=matrix(c(c(-T_alpha,rep(0,d2-s)),c(T_alpha,rep(0,d2-s))),d2,length(pro),byrow=F)
rho2=0
## to determin the number of the covariate in each model
h=10
## the size of bootstrap or the size of sub model
B=d2/h
## the sample size in each bootstrap
n_a=n=100
##### the size of local machine
M_a=M=N/n
######## the tuning parameter
ll1=10
ll2=8
ll3=1
Lamda1=seq(0.000001,2,length.out=ll1)
Lamda2=c(0.000001,0.00001,0.0001,seq(0.001,0.25,length.out=ll2-3))
Lamda3=seq(2,2,length.out=ll3)
## for intial value

#### number of boostrap  for estimating intial value
inti_B=20

###### the iteration of admm
iter_admm=300
## the iteration of greedy algorithm
iter_fw=100
## the pre tuning for ADMM
eps=0.01
### special for grouping
eps_group=2.5
sparse=0.001
gamma=3
#####rho
kappa=1
penalty_index=3
M_T=1
n_T=N

#library("lpSolve")
library("MASS")
library("stats")  
#library("glmnet") 
library("ncvreg")
library("Matrix")
library("Rcpp")
library("RcppArmadillo")
#library("mixtools")
library("admscad")
#library(microbenchmark)


source("data.R")
source("fun.R")
source("intial.R")
source("admm.R")
source("weight.R")
source("group.R")
source("Theta.R")
source("final.R")
source("groupt.R")
source("groupTMA.R")
source("tranf.R")





T_label=matrix(0,N,Nsim)
Sam_Z=array(NA,c(h,B,Nsim))
Sam_Y=array(NA,c(n,M,Nsim))
## distributed






#################### Total model averaging##########
##########estimated values
Hat_theta_T=array(NA,c(d2,N,Nsim))
Hat_alpha_T=array(NA,c(d2,N,Nsim))
Hat_size_T=matrix(NA,1,Nsim)
Hat_cluster_T=array(NA,c(N,N,Nsim))
Hat_weight_T=array(NA,c(B,Nsim))
Hat_label_T=matrix(NA,N,Nsim)
####label
Accuracy_T=matrix(NA,1,Nsim)
####theta
BIAS_T=array(NA,c(d2,KK,Nsim))
MEDIAN_T=array(NA,c(d2,KK,Nsim))
SD_T=array(NA,c(d2,KK,Nsim))
mean_size_T=matrix(NA,1,Nsim)
sd_size_T=matrix(NA,1,Nsim)


OPT_T=OPT=OPT_t=OPT_a=matrix(NA,Nsim,3)
Tim_a=Tim=Tim_t=Tim_T=rep(NA,Nsim)


temp_sam_z=matrix(0,h,B)
temp_sam_y=matrix(0,n,M)


######## Total distributed
#temp_hat_theta_T=array(NA,c(d2,N,ll1,ll2,ll3,Nsim))
#temp_hat_alpha_T=array(NA,c(d2,N,ll1,ll2,ll3,Nsim))
#temp_hat_size_T=array(NA,c(ll1,ll2,ll3,Nsim))
#temp_bic_T=array(NA,c(ll1,ll2,ll3,Nsim))
#temp_hat_label_T=array(NA,c(N,N,ll1,ll2,ll3,Nsim))
#temp_hat_weight_T=array(NA,c(B,ll1,ll2,ll3,Nsim))

Intial_theta_T=array(0,c(h,N,B,ll3))
gc()


BIAS_G=array(NA,c(d2,KK,Nsim))
MEDIAN_G=array(NA,c(d2,KK,Nsim))
SD_G=array(NA,c(d2,KK,Nsim))
Hat_alpha_G=array(NA,c(d2,KK,Nsim))
Hat_weight_G=array(NA,c(B,Nsim))

Hat_alpha_c=array(0,c(d2,KK,Nsim))
Hat_mix_c=matrix(0,nrow=KK,ncol=Nsim)
Hat_label_c=matrix(NA,N,Nsim)
library("fmrs")

library("mixtools")
library("HDclassif")
library("RSKC")
###true
Hat_alpha_tr=array(0,c(d2,KK,Nsim))
Hat_mix_tr=matrix(0,nrow=KK,ncol=Nsim)
Hat_label_tr=matrix(NA,N,Nsim)

## all
Hat_alpha_all=array(0,c(d2,KK,Nsim))
Hat_mix_all=matrix(0,nrow=KK,ncol=Nsim)
Hat_label_all=matrix(NA,N,Nsim)

## partly
Hat_alpha_pa=array(0,c(d2,KK,Nsim))
Hat_mix_pa=matrix(0,nrow=KK,ncol=Nsim)
Hat_label_pa=matrix(NA,N,Nsim)

###  cluster

Hat_alpha_cl=array(0,c(d2,KK,Nsim))
Hat_mix_cl=matrix(0,nrow=KK,ncol=Nsim)
Hat_label_cl=matrix(NA,N,Nsim)

### sparse kmeans

Hat_alpha_skm=array(0,c(d2,KK,Nsim))
Hat_mix_skm=matrix(0,nrow=KK,ncol=Nsim)
Hat_label_skm=matrix(NA,N,Nsim)

### kmeans

Hat_alpha_km=array(0,c(d2,KK,Nsim))
Hat_mix_km=matrix(0,nrow=KK,ncol=Nsim)
Hat_label_km=matrix(NA,N,Nsim)

Hat_alpha_WT=array(NA,c(d2,N,Nsim))
DATA_X=array(NA,c(d2,N,Nsim))
DATA_Y=matrix(NA,N,Nsim)
Temp_hat_theta_T=array(NA,c(d2,N,ll1,ll2,ll3,Nsim))
Temp_hat_alpha_T=array(NA,c(d2,N,ll1,ll2,ll3,Nsim))
Temp_hat_size_T=array(NA,c(ll1,ll2,ll3,Nsim))
Temp_bic_T=array(NA,c(ll1,ll2,ll3,Nsim))
Temp_hat_label_T=array(NA,c(N,N,ll1,ll2,ll3,Nsim))
Temp_hat_weight_T=array(NA,c(B,ll1,ll2,ll3,Nsim))
Temp_Intial_theta_T=array(0,c(h,N,B,ll3,Nsim))

for (i in 1:Nsim){
#T_group=sample(1:KK,N,replace=T,prob=pro)
#T_label[,i]=T_group
#data=gdat(N,rho2,theta0,T_group)
#Y=data$Y
#Z=data$Z

## total Chen
data=fmrs.gendata1(nObs =N, nComp=KK, nCov =d2,
                     coeff =c(0,theta0[,1],0,theta0[,2]), dispersion =rep(0.5,KK),
                     mixProp =pro, rho = rho2,
                     disFamily = "norm")
T_group=data$T_group
T_label[,i]=T_group
Y=data$y
Z=data$x
DATA_X[,,i]=Z
DATA_Y[,i]=Y

temp_sam_y=matrix(1:N, nrow=n,byrow=F)
#temp_sam_z=matrix(1:p,nrow=h,byrow=F)
temp_sam_z=matrix(0,nrow=h,ncol=B)
ind_sam=sample((s+1):d2,d2-s,replace=F)

temp_sam_z=matrix(c(1:s,ind_sam),nrow=h,byrow=F)

sink("sam100.R",append=TRUE)
print(i)
print(temp_sam_z)
sink()
Sam_Z[,,i]=temp_sam_z
Sam_Y[,,i]=temp_sam_y

res_tr=regmixEM(Y,Z[,temp_sam_z[1:s,1]],k=KK,addintercept=FALSE)

###true
Hat_alpha_tr[temp_sam_z[1:s,1],,i]=res_tr$beta
Hat_mix_tr[,i]=res_tr$lambda
Hat_label_tr[,i]=apply(res_tr$posterior,1,which.max)

## all
res_all=regmixEM(Y,Z[,temp_sam_z[,1]],k=KK,addintercept=FALSE)
Hat_alpha_all[temp_sam_z[,1],,i]=res_all$beta
Hat_mix_all[,i]=res_all$lambda
Hat_label_all[,i]=apply(res_all$posterior,1,which.max)

## partly
res_pa=regmixEM(Y,Z[,c(temp_sam_z[1:2,2],temp_sam_z[3:(h),1])],k=KK,addintercept=FALSE)
Hat_alpha_pa[c(temp_sam_z[1:2,2],temp_sam_z[3:(h),1]),,i]=res_pa$beta
Hat_mix_pa[,i]=res_pa$lambda
Hat_label_pa[,i]=apply(res_pa$posterior,1,which.max)


##  cluster
res_cl=hddc(cbind(Y,Z),K=KK)

Hat_mix_cl[,i]=res_cl$prop
Hat_label_cl[,i]=res_cl$class
for(imm in 1:KK){
index_cl=which(res_cl$class==imm)
Z_cl=Z[index_cl,]
res_Gcl=cv.ncvreg(Z_cl,Y[index_cl])
lam_Gcl=res_Gcl$fit
out_Gcl=lam_Gcl$beta[,res_Gcl$min]
Hat_alpha_cl[,,i]=out_Gcl[-1]
}




## sparse kmeans

res_skm=RSKC(cbind(Y,Z),KK, alpha = 0,L1=1)

Hat_label_skm[,i]=res_skm$labels
for(imm in 1:KK){
index_skm=which(res_skm$labels==imm)
Hat_mix_skm[imm,i]=length(index_skm)/N
Z_skm=Z[index_skm,]
res_Gskm=cv.ncvreg(Z_skm,Y[index_skm])
lam_Gskm=res_Gskm$fit
out_Gskm=lam_Gskm$beta[,res_Gskm$min]
Hat_alpha_skm[,,i]=out_Gskm[-1]
}


##kmeans

res_km=kmeans(Y,KK)
Hat_label_km[,i]=res_km$cluster
for(imm in 1:KK){
index_km=which(res_km$cluster==imm)
Hat_mix_km[imm,i]=length(index_km)/N
Z_km=Z[index_km,]
res_Gkm=cv.ncvreg(Z_km,Y[index_km])
lam_Gkm=res_Gkm$fit
out_Gkm=lam_Gkm$beta[,res_Gkm$min]
Hat_alpha_km[,,i]=out_Gkm[-1]
}

## total Chen
res.mle <- fmrs.mle(y = Y, x =Z,delta = data$delta, nComp =KK, disFamily = "norm",
                   initCoeff = rnorm(KK*p+KK),
                   initDispersion = rep(1,KK),
                   initmixProp = rep(1/KK, KK),nIterNR = 200)
res.lam <- fmrs.tunsel(y = Y, x =Z,delta = data$delta,nComp =KK, disFamily = "norm",
                      initCoeff = c(coefficients(res.mle)),
                      initDispersion = dispersion(res.mle),
                      initmixProp = mixProp(res.mle),
                      penFamily = "adplasso",nIterNR = 200)
#show(res.lam)

res.var <- fmrs.varsel(y = Y, x =Z, delta = data$delta,
                      nComp = ncomp(res.mle), disFamily = "norm",
                      initCoeff=c(coefficients(res.mle)),
                      initDispersion = dispersion(res.mle),
                      initmixProp = mixProp(res.mle),
                      penFamily = "adplasso",
                      lambPen = slot(res.lam, "lambPen"),nIterNR = 200)
Hat_alpha_c[,1:KK,i]=round(coefficients(res.var)[-1,],3)
Hat_mix_c[,i]=mixProp(res.var)
Hat_label_c[,i]=apply(weights(res.var),1,which.max)


gc()

Z_matrix_T=array(0,c(N,N*h,B))
Z_inv_T=array(0,c(N*h,N*h,B))
Delta_t=Matrix(0,nrow=N,ncol=(N*(N-1)/2),sparse=T)
e_t=diag(N)
for(ix in 1:(N-1)){
ind=(ix-1)*N-ix*(ix-1)/2
Delta_t[,(ind+1):(ind+N-ix)]=e_t[,ix]-e_t[,(ix+1):N]
for(is in 1:B){
Z_matrix_T[ix,((ix-1)*h+1):((ix-1)*h+h),is]=Z[ix,temp_sam_z[,is]]
gc()
}
}
gc()
Delta_t=t(Delta_t)
gc()
#Delta_t=matrix(Delta_t,ncol=N,byrow=F)
A_T=kronecker(as.matrix(Delta_t),diag(1,h))
Comp_G=matrix(0,N,B)
Theta_G=array(0,c(d2,KK,B))
for(is in 1:B){
Z_matrix_T[N,((N-1)*h+1):((N-1)*h+h),is]=Z[N,temp_sam_z[,is]]
Z_inv_T[,,is]=solve(t(Z_matrix_T[,,is])%*%Z_matrix_T[,,is]+kappa*t(A_T)%*%A_T+kappa*diag(N*h))
for(imm in 1:KK){
index_G=which(T_group==imm)
Z_G=Z[index_G,temp_sam_z[,is]]
res_G=cv.ncvreg(Z_G,Y[index_G])
lam_G=res_G$fit
out_G=lam_G$beta[,res_G$min]
Theta_G[temp_sam_z[,is],imm,is]=out_G[2:(h+1)]
Comp_G[index_G,is]=Z_G%*%out_G[2:(h+1)]
}



 for(ll in 1:ll3){
    Intial_theta_T[,,is,ll]= f_intial(Z[,temp_sam_z[,is]],Y,Lamda3[ll],inti_B)$the}
gc()
}

Weight_G=f_weight(Y,Comp_G,iter_fw,eps,B)$w
Hat_weight_G[,i]=Weight_G
for(imm in 1:KK){
Hat_alpha_G[,imm,i]=Theta_G[,imm,]%*%Weight_G
}
gc() 

timestart_T<-Sys.time()
res_group_T=f_group_T(M_T,n_T,Y,Z,Z_matrix_T,Z_inv_T,A_T,d2,h,N,B,temp_sam_z,temp_theta,Lamda1,Lamda2,kappa,gamma,iter_admm,eps,penalty_index,iter_fw,Intial_theta_T)
timeend_T<-Sys.time() 
Tim_T[i]=timeend_T-timestart_T

temp_hat_theta_T=res_group_T$hat_theta 
temp_hat_alpha_T=res_group_T$hat_alpha 
temp_hat_size_T=res_group_T$hat_size 
temp_bic_T=res_group_T$bic 
temp_hat_label_T=res_group_T$hat_label
temp_hat_weight_T=res_group_T$hat_weight 
temp_opt_T=which(temp_bic_T==min(temp_bic_T),arr.ind=T)

opt_T=matrix(1,nrow=1,ncol=3)
if(ll2==1&ll3==1){opt_T[1,1]=temp_opt_T}
if(ll2==1&ll3!=1){opt_T[1,]=c(temp_opt_T[1,1],1,temp_opt_T[1,2])}
if(ll2!=1&ll3==1){opt_T[1,]=c(temp_opt_T[1,1],temp_opt_T[1,2],1)}
OPT_T[i,]=opt_T[1,]
hat_alpha_T=temp_hat_alpha_T[,,opt_T[1,1],opt_T[1,2],opt_T[1,3]]
hat_size_T=temp_hat_size_T[opt_T[1,1],opt_T[1,2],opt_T[1,3]]
hat_cluster_T=temp_hat_label_T[,,opt_T[1,1],opt_T[1,2],opt_T[1,3]]
Hat_weight_T[,i]=temp_hat_weight_T[,opt_T[1,1],opt_T[1,2],opt_T[1,3]]
res_tranf_T=f_tranf_c(KK,hat_size_T,hat_alpha_T,theta0,hat_cluster_T,N)
Hat_theta_T[,,i]=temp_hat_theta_T[,,opt_T[1,1],opt_T[1,2],opt_T[1,3]]
Hat_alpha_T[,,i]=res_tranf_T$Alpha
Hat_size_T[1,i]=hat_size_T
Hat_cluster_T[,,i]=res_tranf_T$Cluster
Hat_label_T[,i]=res_tranf_T$Label
gc()
Temp_hat_theta_T[,,,,,i]=temp_hat_theta_T
Temp_hat_alpha_T[,,,,,i]=temp_hat_alpha_T
Temp_hat_size_T[,,,i]=temp_hat_size_T
Temp_bic_T[,,,i]=temp_bic_T
Temp_hat_label_T[,,,,,i]=temp_hat_label_T
Temp_hat_weight_T[,,,,i]=temp_hat_weight_T
Temp_Intial_theta_T[,,,,i]= Intial_theta_T

save_a=Hat_label_T[,1:i]
save_b=T_label[,1:i]
save_c=Hat_alpha_T[,1:max(Hat_size_T[1,i]),1:i]
save_d=Hat_weight_T[,1:i]
save_e=Hat_alpha_G[,,1:i]
save_f=Hat_weight_G[,1:i]

save(theta0,save_a,save_b,save_c,save_d,save_e,save_f,file="/gpfs/ysm/project/bh563/subgroup/2/resT2019.RData")


gc()







gc()




Accuracy_T[i]=mean(Hat_label_T[,i]==T_label[,i])
### reslut of theta
for(ix in 1:KK){
  for(iy in 1:s){
    
  
    BIAS_T[iy,ix,i]=(Hat_alpha_T[iy,ix,i]-theta0[iy,ix])
    BIAS_G[iy,ix,i]=(Hat_alpha_G[iy,ix,i]-theta0[iy,ix])
    
    if(i>1){
      
      MEDIAN_T[iy,ix,i]=median(Hat_alpha_T[iy,ix,1:i])
      SD_T[iy,ix,i]=sd(Hat_alpha_T[iy,ix,1:i])
      MEDIAN_G[iy,ix,i]=median(Hat_alpha_G[iy,ix,1:i])
      SD_G[iy,ix,i]=sd(Hat_alpha_G[iy,ix,1:i])
    }else{
    MEDIAN_T[iy,ix,i]=mean(Hat_alpha_T[iy,ix,1:i])
    SD_T[iy,ix,i]=0
    MEDIAN_G[iy,ix,i]=mean(Hat_alpha_G[iy,ix,1:i])
    SD_G[iy,ix,i]=0
    }
    
    cat(c("total average","\n"))
    cat(c("i","group","dim","bias","median","SD","\n"))
    cat(c(i,ix,iy,BIAS_T[iy,ix,i], MEDIAN_T[iy,ix,i],SD_T[iy,ix,i],"\n"))
  }
  
}

if(i>1){
  cat(c("i","size","\n"))
 
  cat(c("total average","\n"))
  cat(c("i","mean","median","sd","\n"))
  cat(c(i,mean(Hat_size_T[1,1:i]),median(Hat_size_T[1,1:i]),sd(Hat_size_T[1,1:i]),"\n"))
  cat(c("i","Accuracy","\n"))
  cat(c(i,mean(Accuracy_T[1:i]),"\n"))

} 


sink("res100.R",append=TRUE)
Temp_Bias_T=matrix(0,s,KK)
Temp_Bias_G=matrix(0,s,KK)

for(ix in 1:KK){
  for(iy in 1:s){
    
    Temp_Bias_T[iy,ix]=mean(BIAS_T[iy,ix,1:i])
    Temp_Bias_G[iy,ix]=mean(BIAS_G[iy,ix,1:i])
  }
} 
cat(c(i,'\n'))

cat(c(round(Temp_Bias_T,4),'\n'))

cat(c(round(SD_T[1:s,1:KK,i],4),'\n'))
cat(c(round(Temp_Bias_G,4),'\n'))

cat(c(round(SD_G[1:s,1:KK,i],4),'\n'))


cat(c("meansize","mediansize","accuray",'\n'))
cat(c(mean(Hat_size_T[1,1:i]),'\n')) 
cat(c(median(Hat_size_T[1,1:i]),'\n')) 
cat(c(mean(Accuracy_T[1:i]),'\n'))

sink()

gc()

save.image(file="/gpfs/ysm/scratch60/bh563/subgroup/2/res2019.RData")
}





