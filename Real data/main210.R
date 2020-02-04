#rm(list=ls(all=TRUE))
#setwd("D:/subgroup/R programming/R program summary/Real data/selection/210")
#Rcpp::sourceCpp("admmcpp.cpp")
rm(list=ls(all=TRUE))
#setwd("D:/subgroup/Realdata/1")
set.seed(99)
KK=2
## to determin the number of the covariate in each model
h=10
## the size of bootstrap or the size of sub model
#B=d2/h
B=25
hh=50
data=read.table("BCRA_imageFeature.txt",header=T)
rownames(data)=data[,1]
data=data[,-1]
BCRA_Y=read.csv("BCRA_clinical_data_y.csv")
select_Y=read.csv("select_Y.csv")
#BCRA_Y=read.csv("BCRA_y.csv")
###????not available ????0?y
Y1=select_Y$LYMPH_NODE_EXAMINED_COUNT
Y2=select_Y$LYMPH_NODES_EXAMINED_HE_COUNT
Ind=which((Y1!="[Not Available]"&Y2!="[Not Available]"&Y2!=0)==T)
select_Y=select_Y[Ind,]
sel_Y=as.numeric(as.character(Y2))/as.numeric(as.character(Y1))
sel_Y=log(sel_Y/(1-sel_Y))
sel_ind=which((is.infinite(sel_Y)==F&is.na(sel_Y)==F)==T)
select_Y=select_Y[sel_ind,]
 ### ?image?????? 
newColname=sapply(colnames(data), function(a){gsub("\\.","-", substr(a, 1,12))})
ind=pmatch(newColname,select_Y[,2])
  ## ??data ????
ind1=which(is.na(ind)==F)
data=data[,ind1]
 ## ??NA
dat=apply(data, 1, function(a){a[is.na(a)] = median(as.numeric(a), na.rm = T); a})
## column and row rotation
  ### all column are 0
dat=dat[,-c(291,293)]
  ## ??????0?
med=which(apply(dat,2,median)==0)
dat=dat[,-c(med)]
  ## normalization
dat=apply(dat,2,scale)
  ## ?????????
Ind=unique(which(is.na(dat)==T,arr.ind=T)[,2])
dat=dat[,-Ind]
  ## ?????
which(colnames(dat)=="Neighbors_SecondClosestDistance_Adjacent")
which(colnames(dat)=="AreaShape_EulerNumber")
which(colnames(dat)=="Neighbors_FirstClosestDistance_Adjacent")
dat=dat[,-c(14,34,104)]
  ## ??y
ind2=ind[is.na(ind)==F]
select_Y=select_Y[ind2,]
  ## y???
Y1=select_Y$LYMPH_NODE_EXAMINED_COUNT
Y2=select_Y$LYMPH_NODES_EXAMINED_HE_COUNT
sel_Y=as.numeric(as.character(Y2))/as.numeric(as.character(Y1))
sel_Y=log(sel_Y/(1-sel_Y))
sel_ind=which(is.infinite(sel_Y)==F)
dat=dat[sel_ind,]
sel_Y=sel_Y[sel_ind]
set.seed(99)

#nCores <- Sys.getenv("TOTALPROCS")
#nCores <-Sys.getenv("NUMBER_OF_PROCESSORS")
nCores <- 3

print(nCores)
library(foreach)
library(doMC)
registerDoMC(nCores)
library(doRNG)
registerDoRNG(1)
Nsim=1
### the total sample size
N=nrow(dat)

### the total dimension
p=ncol(dat)
### the sparsity
s=10
### size of group

### the probabilty of different group
pro=rep(1/KK,KK)
## the dimension of subject-specific covariates
d2=p
alpha0=1
#T_alpha=runif(s,0.8,1)
T_alpha=rep(alpha0,s)
theta0=matrix(c(c(-T_alpha,rep(0,d2-s)),c(T_alpha,rep(0,d2-s))),d2,length(pro),byrow=F)
rho2=0

## the sample size in each bootstrap
n_a=n=N
##### the size of local machine
M_a=M=N/n
######## the tuning parameter
ll1=10
ll2=8
ll3=1
Lamda1=seq(0.000001,1.7,length.out=ll1)
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
library("flexmix")
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
Temp_hat_theta_T=array(NA,c(d2,N,ll1,ll2,ll3,Nsim))
Temp_hat_alpha_T=array(NA,c(d2,N,ll1,ll2,ll3,Nsim))
Temp_hat_size_T=array(NA,c(ll1,ll2,ll3,Nsim))
Temp_bic_T=array(NA,c(ll1,ll2,ll3,Nsim))
Temp_hat_label_T=array(NA,c(N,N,ll1,ll2,ll3,Nsim))
Temp_hat_weight_T=array(NA,c(B,ll1,ll2,ll3,Nsim))
Temp_Intial_theta_T=array(0,c(h,N,B,ll3,Nsim))

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
BIC_tr=BIC_all=BIC_pa=BIC_cl=BIC_skm=BIC_km=BIC_c=BIC_T=rep(0,Nsim)
MSE_tr=MSE_all=MSE_pa=MSE_cl=MSE_skm=MSE_km=MSE_c=MSE_T=rep(0,Nsim)
ID=matrix(0,Nsim,p)
T_SAM=matrix(0,139,Nsim)
SEL_Y=sel_Y
DAT=dat
for (i in 1:Nsim){
#T_group=sample(1:KK,N,replace=T,prob=pro)
#T_label[,i]=T_group
#data=gdat(N,rho2,theta0,T_group)
#Y=data$Y
#Z=data$Z

## total Chen
#data=fmrs.gendata1(nObs =N, nComp=KK, nCov =d2,
#                     coeff =c(0,theta0[,1],0,theta0[,2]), dispersion =rep(0.5,KK),
#                     mixProp =pro, rho = rho2,
#                     disFamily = "norm")
#T_group=data$T_group
#T_label[,i]=T_group

#Y=data$y
#Z=data$x 
test_sam=1:139
T_SAM[,i]=test_sam
sel_Y=SEL_Y[test_sam]
dat=DAT[test_sam,]
T_group=rep(2,N)
for(xx in 1:N)
{T_group[xx]=as.numeric(select_Y$RACE[xx]=="WHITE")}
T_group=T_group+1
T_group_test=T_group[-c(1:N)]
T_group=T_group[c(1:N)]
#Y_test=sel_Y[-c(1:N)]
Y=sel_Y[c(1:N)]
#Z_test=dat[-c(1:N),]
Z=(dat[c(1:N),])

DATA_X[,,i]=Z
DATA_Y[,i]=Y

## screening
par=rep(0,p)
for(j in 1:p){res_par=flexmix(Y~Z[,j]-1,k=KK)
              par_res=parameters(res_par)
              par[j]=sum(abs(par_res[1,]))
     }

id=order(par,decreasing=T)
ID[i,]=id
Z=Z[,id]

#Z_test=Z_test[,id]
temp_sam_y=matrix(1:N, nrow=n,byrow=F)
#temp_sam_z=matrix(1:p,nrow=h,byrow=F)
temp_sam_z=matrix(0,nrow=h,ncol=B)
#ind_sam=sample((s+1):d2,d2-s,replace=F)
#temp_sam_z=matrix(c(1:s,ind_sam,sample(1:d2,B*h-d2)),nrow=h,byrow=F)

temp_sam_z=matrix(c(1:d2,sample(1:d2,B*h-d2)),nrow=h,byrow=F)

sink("sam100.R",append=TRUE)
print(i)
print(temp_sam_z)
sink()
Sam_Z[,,i]=temp_sam_z
Sam_Y[,,i]=temp_sam_y
Z_test=Z
Y_test=Y
###true




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
#Comp_G=matrix(0,N,B)
#Theta_G=array(0,c(d2,KK,B))
for(is in 1:B){
Z_matrix_T[N,((N-1)*h+1):((N-1)*h+h),is]=Z[N,temp_sam_z[,is]]
Z_inv_T[,,is]=solve(t(Z_matrix_T[,,is])%*%Z_matrix_T[,,is]+kappa*t(A_T)%*%A_T+kappa*diag(N*h))




 for(ll in 1:ll3){
    Intial_theta_T[,,is,ll]= f_intial(Z[,temp_sam_z[,is]],Y,Lamda3[ll],inti_B)$the}
gc()
}

#Weight_G=f_weight(Y,Comp_G,iter_fw,eps,B)$w
#Hat_weight_G[,i]=Weight_G
#for(imm in 1:KK){
#Hat_alpha_G[,imm,i]=Theta_G[,imm,]%*%Weight_G
#}
gc() 
 #theta0=Hat_alpha_G[,,i]
timestart_T<-Sys.time()
res_group_T=f_group_T(M_T,n_T,Y,Z,Z_matrix_T,Z_inv_T,A_T,d2,h,N,B,temp_sam_z,temp_theta,Lamda1,Lamda2,kappa,gamma,iter_admm,eps,penalty_index,iter_fw,Intial_theta_T)
timeend_T<-Sys.time() 
Tim_T[i]=timeend_T-timestart_T
#save.image(file="/gpfs/ysm/scratch60/bh563/realdata/h20/2/res10.RData")

temp_hat_theta_T=res_group_T$hat_theta 
temp_hat_alpha_T=res_group_T$hat_alpha 
temp_hat_size_T=res_group_T$hat_size 
temp_bic_T=res_group_T$bic 
temp_hat_label_T=res_group_T$hat_label
temp_hat_weight_T=res_group_T$hat_weight 
temp_opt_T=which(temp_bic_T==min(temp_bic_T),arr.ind=T)

Temp_hat_theta_T[,,,,,i]=temp_hat_theta_T
Temp_hat_alpha_T[,,,,,i]=temp_hat_alpha_T
Temp_hat_size_T[,,,i]=temp_hat_size_T
Temp_bic_T[,,,i]=temp_bic_T
Temp_hat_label_T[,,,,,i]=temp_hat_label_T
Temp_hat_weight_T[,,,,i]=temp_hat_weight_T
Temp_Intial_theta_T[,,,,i]= Intial_theta_T
#save.image(file="/gpfs/ysm/scratch60/bh563/subgroup/1/restest2019.RData")

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
BIC_T[i]=Temp_bic_T[opt_T[1,1],opt_T[1,2],opt_T[1,3],i]
MSE_T[i]=sum((apply(abs(Y_test-(Z_test%*%hat_alpha_T[,1:hat_size_T])),1,min))^2/98)
gc()

save_a=Hat_label_T[,1:i]
save_b=T_label[,1:i]
save_c=Hat_alpha_T[,1:max(Hat_size_T[1,i]),1:i]
save_d=Hat_weight_T[,1:i]

#save.image(file="res210.RData")
}





