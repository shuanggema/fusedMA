f_log=function(b,a){all(b==a)}
f_norm=function(x){norm(x,"2")}

###????????
f_trun=function(b,a){norm(b-a,"2")<=eps_group}
f_which=function(b,a){which(b==a)}
f_label_len=function(a){
  a=a[a!=0]
  return(length(a))
}
f_label_min=function(a){
  a=a[a!=0]
  return(min(a))
}
S=function(z,t){
  normdelta=sqrt(sum(z^2))
  if((1-t/normdelta)>0){
    a=(1-t/normdelta)*z}else{
    a=rep(0,length(z))}
    return(a)
    }
 ############ scad penalty
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
############ scad penalty
f_scad=function(lamda,kappa,zeta,gamma){
  normdelta<-sqrt(sum(zeta^2))
  if(normdelta<=(lamda+lamda/kappa)){
    a=S(zeta,lamda/kappa)}else if(normdelta>gamma*lamda){
      a=zeta}else{
        a=S(zeta,gamma*lamda/((gamma-1)*kappa))/(1-1/((gamma-1)*kappa))}
  return(a)
}
  
   ST=function(z,t){ 
      res=pmax(0,z-t)-pmax(0,-z-t)
      return(res)}
  