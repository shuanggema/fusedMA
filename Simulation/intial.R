###??

f_intial=function(Y,Z_matrix,A,N,p,lamda3,inti_size,inti_B){
    c1=array(0,c(p,N,ll3))
    diff=rep(0,ll3)
    for(i in 1:ll3){
    hat=solve(t(Z_matrix)%*%Z_matrix+lamda3[i]*t(A)%*%A)%*%t(Z_matrix)%*%Y
    c1[,,i]=matrix(hat,nrow=p,byrow=F)
    diff[i]=sum((Y-Z_matrix%*%hat)^2)
    }
    j=which.min(diff)
    res=list(the=c1[,,j])
    return(res)
    gc()
   }
   

##the EM algorithm to fit the mixture of linear regression
#mixlinone estimates the mixture regression parameters by MLE based on ONE
#initial value
mixlinone<-function(x,y,bet,sig,pr,m=2){
run=0; n=length(y);
# X=cbind(rep(1,n),x);
X<-x;
if(length(sig)>1 ){ #the case when the variance is unequal

r=matrix(rep(0,m*n),nrow=n);pk=r;lh=0;
for(j in seq(m))
{r[,j]=y-X%*%bet[j,];lh=lh+pr[j]*dnorm(r[,j],0,sig[j]);}
lh=sum(log(lh));
#E-steps
repeat
{ prest=c(bet,sig,pr);run=run+1;plh=lh;
for(j in seq(m))
{ pk[,j]=pr[j]*pmax(10^(-300),dnorm(r[,j],0,sig[j]))}
pk=pk/matrix(rep(apply(pk,1,sum),m),nrow=n);
#M-step
np=apply(pk,2,sum);pr=np/n;lh=0;
for(j in seq(m))
{w=diag(pk[,j]);
bet[j,]=ginv(t(X)%*%w%*%X)%*%t(X)%*%w%*%y;
r[,j]= y-X%*%bet[j,]; sig[j]=sqrt(t(pk[,j])%*%(r[,j]^2)/np[j]);
lh=lh+pr[j]*dnorm(r[,j],0,sig[j]);}
lh=sum(log(lh));dif=lh-plh;
if(dif<10^(-5)|run>500){break}}}
else{ #the case when the variance is equal
r=matrix(rep(0,m*n),nrow=n);pk=r; lh=0
for(j in seq(m))
{r[,j]=y-X%*%bet[j,];lh=lh+pr[j]*dnorm(r[,j],0,sig);}
lh=sum(log(lh));

#E-steps
repeat
{ prest=c(bet,sig,pr);run=run+1;plh=lh;
for(j in seq(m))
{ pk[,j]=pr[j]* pmax(10^(-300),dnorm(r[,j],0,sig)) }
pk=pk/matrix(rep(apply(pk,1,sum),m),nrow=n);
#M-step
np=apply(pk,2,sum);pr=np/n;
for(j in seq(m))
{ w=diag(pk[,j]);
bet[j,]=ginv(t(X)%*%w%*%X)%*%t(X)%*%w%*%y;
r[,j]= y-X%*%bet[j,]; }
sig=sqrt(sum(pk*(r^2))/n);lh=0;
for(j in seq(m))
{lh=lh+pr[j]*dnorm(r[,j],0,sig);}
lh=sum(log(lh));
dif=lh-plh;
if(dif<10^(-5)|run>500){break}}
sig=sig*rep(1,m)}
est=list(theta= matrix(c(bet,sig,pr),nrow=m),likelihood=lh,run=run,diflh=dif,pk=pk)
est}
##mixlin based on 20 initial values
mixlin <-function(x,y,k=2,numini=20)

{ n=length(y); x=matrix(x,nrow=n)
a=dim(x);p=a[2]; #n1=2*p;
n1=0.5*n
bet= matrix(rep(0,k*p),nrow=k);sig=0;
for(j in seq(k))
{ind=sample(1:n,n1); 
# X=cbind(rep(1,n1),x[ind,]);
X=x[ind,];
bet[j,]=ginv(t(X)%*%X)%*%t(X) %*%y[ind];
sig=sig+sum((y[ind] -X%*%bet[j,])^2);}
pr=rep(1/k,k);sig=sig/n1/k;
est=mixlinone(x,y,bet,sig,pr,k);lh=est$likelihood;
obj=rep(0,numini); obj[1]=lh;
for(i in seq(numini-1))
{bet= matrix(rep(0,k*p),nrow=k);sig=0;
for(j in seq(k))
{ind=sample(1:n,n1); 
# X=cbind(rep(1,n1),x[ind,]);
X=x[ind,];
bet[j,]=ginv(t(X)%*%X)%*%t(X) %*%y[ind];
sig=sig+sum((y[ind] -X%*%bet[j,])^2);}
pr=rep(1/k,k);sig=sig/n1/k;pest=est;plh=lh;
est=mixlinone(x,y,bet,sig,pr,k);lh=est$likelihood;obj[i+1]=lh;
if(lh<plh){est=pest;lh=plh;}}
est=list(theta=est$theta,likelihood=est$likelihood,run=est$run,diflh=est$dif,objlh=obj,pk=est$pk)
est}

 #res2=cv.ncvreg(Z[index1,],Y[index1])
  #      lam2=res2$fit
   #     out2=lam2$beta[,res2$min]

 f_intial=function(Z,Y,lamda3,inti_B){
  res=mixlin(Z,Y,k=lamda3,numini=inti_B)
  Label=apply(res$pk,1,which.max)
  #res=regmixEM(Y,Z,addintercept =F)
  #Label=apply(res$posterior,1,which.max)
  nc=ncol(Z)
  nr=nrow(Z)
  C1=matrix(0,nc,nr)
  for(ini in 1:lamda3){
  index=which(Label==ini)
  if(length(index)>=1){
  res2=cv.ncvreg(Z[index,],Y[index])
  lam2=res2$fit
  out2=lam2$beta[,res2$min]
  C1[,index]=out2[2:(nc+1)]}
  }
  res=list(the=C1)
  return(res)
    gc() 
  }
       