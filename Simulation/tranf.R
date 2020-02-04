           
f_tranf_c=function(K,Size,Alpha,newtheta0,Cluster,n){
      ind=seq(0,0,length=K)
      Label=rep(0,n)
      if(Size>1){
       for(iy in 1:K){
       diff_alpha=apply(Alpha[,1:Size]-newtheta0[1:s,iy],2,f_norm)
       com=cbind(diff_alpha,1:Size)
       com=com[sort.list(com[,1]),]
       ind[iy]=com[1,2]
     }
     }
     if(Size==1){ind=1}
     ind=unique(ind)
     alpha1=Alpha[,ind]
     alpha2=Alpha[,-ind]
     label1=Cluster[,ind]
     label2=Cluster[,-ind]
     Alpha=cbind(alpha1,alpha2)
     Cluster=cbind(label1,label2)
    ######relabel according true theta0 
     for(iy in 1:Size){
     index=which(Cluster[,iy]==1)
     Label[index]=iy
     }
    res=list(Alpha=Alpha,Cluster=Cluster,Label=Label) 
    return(res)
}