gdat=function(N,rho2,theta0,T_group)
{
  
  
  ## 2 The generation of subject-specific Z 
  Ma2=matrix(rep(1:d2,d2),ncol=d2,byrow=F)
  CovM2=rho2^(abs(Ma2-t(Ma2)))
  Z=mvrnorm(N,rep(0,d2),CovM2)
  eps=rnorm(N,0,0.5)
  ## generate Y
  Y=diag(Z%*%theta0[,T_group])+eps
  
  ##data
  dat=list(Y=Y,Z=Z)
  return(dat)
}


fmrs.gendata1=function(nObs,nComp,nCov,coeff,dispersion,mixProp,rho,umax,disFamily = "lnorm")



{

  if(missing(disFamily)) disFamily = "lnorm"



  if(sum(mixProp) != 1)

    stop("The sum of mixing proportions must be 1.")

  if(sum(dispersion <= 0) != 0)

    stop("Dispersion parameters cannot be zero or negative.")

  if(rho > 1 | rho < -1)

    stop("The correlation cannot be less than -1 or greater thatn 1.")



  mu <- rep(0, nCov)

  Sigma <- diag(nCov)



  for(i in 1:nCov){

    for(j in 1:nCov){

      Sigma[i,j] <- rho^abs(i-j)

    }}



  X <- matrix(rnorm(nCov * nObs), nObs)

  X <- scale(X, TRUE, FALSE)

  X <- X %*% svd(X, nu = 0)$v

  X <- scale(X, FALSE, TRUE)

  eS <- eigen(Sigma, symmetric = TRUE)

  ev <- eS$values

  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), nCov) %*% t(X)

  nm <- names(mu)

  if (is.null(nm) && !is.null(dn <- dimnames(Sigma)))

    nm <- dn[[1L]]

  dimnames(X) <- list(nm, NULL)

  if (nObs == 1)

    cX = drop(X)

  else cX = t(X)

  cX <- scale(cX)

  colnames(cX) <-  paste("X", 1:nCov,sep = ".")



  coef0 <- matrix(coeff, nrow = nComp, ncol = nCov+1, byrow = TRUE)

  mixProp0 <- cumsum(mixProp)



  yobs <-c()

  c <- rep()

  dlt <- c()

  u <- c()

  tobs <- c()



  if(disFamily == "lnorm"){

    for(i in 1:nObs){

      epss <- rnorm(1)

      u1 <- runif(1)

      k = length(which(mixProp0<=u1)) + 1

      u[i] = k

      yobs[i] <- coef0[k,1] + coef0[k,-1] %*% cX[i,] + dispersion[k] * epss



      c[i] <- log(runif(1, 0, umax))

      tobs[i] <- exp(min(yobs[i],c[i]))

      dlt[i] <- (yobs[i] < c[i])*1

    }

  }else if(disFamily=="norm"){

    for(i in 1:nObs){

      epss <- rnorm(1)

      u1 <- runif(1)

      k = length(which(mixProp0<=u1)) + 1

      u[i] = k

      yobs[i] <- coef0[k,1] + coef0[k,-1] %*% cX[i,] + dispersion[k] * epss

      tobs[i] <- yobs[i]

      dlt[i] <- 1

    }

  }else if(disFamily=="weibull"){

    for(i in 1:nObs){

      ext <- log(rexp(1))

      u1 <- runif(1)

      k = length(which(mixProp0<=u1)) + 1

      yobs[i] <- coef0[k,1] + coef0[k,-1] %*% cX[i,] + dispersion[k] * ext



      c[i]<- log(runif(1, 0, umax))

      tobs[i] <- exp(min(yobs[i],c[i]))

      dlt[i] <- (yobs[i] < c[i])*1

    }

  } else{

    stop("The family of sub-distributions are not specified correctly.")

  }

  return(list(y = tobs, delta = dlt, x = cX, disFamily = disFamily,T_group=u))

}