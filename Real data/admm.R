### admm 
## Y_b the observation for response n_m*1 vector
## Z_b the observation for covariate matrix n_m*h matrix
## the intial estimates
## n_m # sample size
## h dimension
## Z_matrix n [number of row] (n_m*h)[numer of colum] matrix
## Z_inv n_m*h n_m*h matrix
## A (n_m*(n_m-1)/2)*h) (n_m*h)  matrix
## len_the=h, lamda1=1,lamda2=0.1,kappa=3,gamma=2,iter_admm=100,eps=0.001, penalty_index=3
ST=function(z,t){ 
      res=pmax(0,z-t)-pmax(0,-z-t)
      return(res)}
S=function(z,t){
  normdelta=sqrt(sum(z^2))
  if((1-t/normdelta)>0){
    a=(1-t/normdelta)*z}else{
    a=rep(0,length(z))}
    return(a)
    }
 ############ scad penalty
 ## zeta vetor :h*1
  f_scad=function(lamda,kappa,zeta,gamma){
     normdelta=sqrt(sum(zeta^2))
     if(normdelta<=(lamda+lamda/kappa)){
     a=S(zeta,lamda/kappa)}else if(normdelta>gamma*lamda){a=zeta}else{
     a=S(zeta,gamma*lamda/((gamma-1)*kappa))/(1-1/((gamma-1)*kappa))}
     return(a)
     } 
      ##########lasso penalty
f_lasso=function(lamda,kappa,zeta){a=S(zeta,lamda/kappa);return(a)}
#########  mcp penalty
f_mcp=function(lamda,kappa,zeta,gamma){
  normdelta=sqrt(sum(zeta^2))
  if(normdelta<=gamma*lamda){
    a=S(zeta,lamda/kappa)/(1-1/(gamma*kappa))}else{
      a<-zeta}
  return(a)
}
  
f_admm=function(Y_b,Z_b,the,len_the,n_m,Z_matrix,A,lamda1,lamda2,kappa,gamma,iter_admm,eps,Z_inv,penalty_index)
              {
              the1=c(the)
              The1=matrix(the1,ncol=n_m,byrow=F)
              Eta1=matrix(A%*%the1,nrow=len_the,byrow=F)
              V1=matrix(0,len_the,ncol(Eta1))
              v1=c(V1)
              eta1=c(Eta1)
              l=0
              diff=11
              h_dual=1
              u1=c(the)
              vebar1=rep(0,n_m*len_the)
              if(penalty_index==3){
              while(diff>eps)
                  {
                   l=l+1
                   the1=Z_inv%*%(t(Z_matrix)%*%Y_b+kappa*t(A)%*%((eta1)-(v1)/kappa)+kappa*(u1-vebar1/kappa))
                    Xi=the1+vebar1/kappa
                    u1=ST(Xi,lamda2/kappa)  
                    vebar1=vebar1+kappa*(the1-u1)
                    Eta1=matrix(A%*%as.numeric(the1),nrow=len_the,byrow=F)
                    Zeta=Eta1+V1/kappa
                   for(j in 1:ncol(Zeta)) 
                    {
                      Eta1[,j]=f_scad(lamda1,kappa,Zeta[,j],gamma)
                      } 
                    eta1=c(Eta1)
                    Atheta=matrix(A%*%as.numeric(the1),nrow=len_the,byrow=F)
                    V1=V1+kappa*(Atheta-Eta1)
                    v1=c(V1)
                    
                    diff=max(sqrt(sum((Atheta-Eta1)^2)),sqrt(sum((the1-u1)^2)))
                    #cat(c(l,diff,"\n"))
                    if(l>=iter_admm){cat(c("the max iteration is reaching","\n"));break}
                      }
                      }
              else if(penalty_index==2){
              while(diff>eps)
                  {
                   l=l+1
                   the1=Z_inv%*%(t(Z_matrix)%*%Y_b+kappa*t(A)%*%((eta1)-(v1)/kappa)+kappa*(u1-vebar1/kappa))
                    Xi=the1+vebar1/kappa
                    u1=ST(Xi,lamda2/kappa)  
                    vebar1=vebar1+kappa*(the1-u1)
                    Eta1=matrix(A%*%as.numeric(the1),nrow=len_the,byrow=F)
                    Zeta=Eta1+V1/kappa
                   for(j in 1:ncol(Zeta)) 
                    {
                      Eta1[,j]=f_mcp(lamda1,kappa,Zeta[,j],gamma)
                      } 
                    eta1=c(Eta1)
                    Atheta=matrix(A%*%as.numeric(the1),nrow=len_the,byrow=F)
                    V1=V1+kappa*(Atheta-Eta1)
                    v1=c(V1)
                    
                    diff=max(sqrt(sum((Atheta-Eta1)^2)),sqrt(sum((the1-u1)^2)))
                    #cat(c(l,diff,"\n"))
                    if(l>=iter_admm){cat(c("the max iteration is reaching","\n"));break}
                      }
              }
                    
               else {
              while(diff>eps)
                  {
                   l=l+1
                   the1=Z_inv%*%(t(Z_matrix)%*%Y_b+kappa*t(A)%*%((eta1)-(v1)/kappa)+kappa*(u1-vebar1/kappa))
                    Xi=the1+vebar1/kappa
                    u1=ST(Xi,lamda2/kappa)  
                    vebar1=vebar1+kappa*(the1-u1)
                    Eta1=matrix(A%*%as.numeric(the1),nrow=len_the,byrow=F)
                    Zeta=Eta1+V1/kappa
                   for(j in 1:ncol(Zeta)) 
                    {
                      Eta1[,j]=f_lasso(lamda1,kappa,Zeta[,j])
                      } 
                    eta1=c(Eta1)
                    Atheta=matrix(A%*%as.numeric(the1),nrow=len_the,byrow=F)
                    V1=V1+kappa*(Atheta-Eta1)
                    v1=c(V1)
                    diff=max(sqrt(sum((Atheta-Eta1)^2)),sqrt(sum((the1-u1)^2)))
                   # cat(c(l,diff,"\n"))
                    if(l>=iter_admm){cat(c("the max iteration is reaching","\n"));break}
                      }
              }
                   res=list(theta_b=u1,eta_b=Eta1,Diff=diff,IT=l)
                   return(res)
 gc()
               ## output u1: n_m*h-vector, Eta1:h (n_m(n_m-1)/2)-matirx,  diff, a scalar
                    } 