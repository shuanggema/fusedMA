 ### ??????? ?weight
 CR_derivate=function(Y_M,Hat_Y_M,w){
    t(Hat_Y_M)%*%(Hat_Y_M%*%w-Y_M)
  }
 
 
  f_weight=function(Y_M,Hat_Y_M,iter_fw,eps,B){
    diff=88
    i=1
    w0=rep(0,B)
    W=matrix(0,B,iter_fw)
    while(i<=iter_fw&&diff>=eps){
      g=CR_derivate(Y_M,Hat_Y_M,w0)
      JK=which.min(g)
      e=rep(0,B)
      e[JK]=1
      w1=w0+2*(1+i)^{-1}*(e-w0)
      W[,i]=w1
      i=i+1
      diff=max(abs(w1-w0))
      w0=w1
    }
    ans=list(w=w1,W=W)
    return(ans)
  }
  