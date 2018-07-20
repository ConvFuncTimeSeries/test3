#' Genreates Two-Regime (VAR) Models
#'
#' Genreates two-regime multivariate vector auto-regressive models.
#' @param nob number of observations.
#' @param thr threshold value.
#' @param phi1 VAR coefficient matrix of regime 1.
#' @param phi2 VAR coefficient matrix of regime 2.
#' @param sigma1 innovational covariance matrix of regime 1.
#' @param sigma2 innovational covariance matrix of regime 2.
#' @param c1 constant vector of regime 1.
#' @param c2 constatn vector of regime 2.
#' @param delay two elements (i,d) with "i" being the component index and "d" the delay for threshold variable.
#' @param ini burn-in period.
#' @return mTAR.sim returns a list with following components:
#' \item{series}{a time series following the two-regime multivariate VAR model.}
#' \item{at}{innovation of the time series.}
#' \item{threshold}{threshold value.}
#' \item{delay}{two elements (i,d) with "i" being the component index and "d" the delay for threshold variable.}
#' \item{n1}{number of observations in regime 1.}
#' \item{n2}{number of observations in regime 2.}
#' @examples
#' phi1=matrix(c(0.5,0.7,0.3,0.2),2,2)
#' phi2=matrix(c(0.4,0.6,0.5,-0.5),2,2)
#' sigma1=matrix(c(1,0,0,1),2,2)
#' sigma2=matrix(c(1,0,0,1),2,2)
#' c1=c(0,0)
#' c2=c(0,0)
#' delay=c(1,1)
#' y=mTAR.sim(100,0,phi1,phi2,sigma1,sigma2,c1,c2,delay,ini=500)
#' @export
"mTAR.sim" <- function(nob,thr,phi1,phi2,sigma1,sigma2=NULL,c1=NULL,c2=NULL,delay=c(1,1),ini=500){
if(!is.matrix(phi1))phi1=as.matrix(phi1)
if(!is.matrix(phi2))phi2=as.matrix(phi2)
if(!is.matrix(sigma1))sigma1=as.matrix(sigma1)
if(is.null(sigma2))sigma2=sigma1
if(!is.matrix(sigma2))sigma2=as.matrix(sigma2)
k1 <- nrow(phi1)
k2 <- nrow(phi2)
k3 <- nrow(sigma1)
k4 <- nrow(sigma2)
k <- min(k1,k2,k3,k4)
if(is.null(c1))c1=rep(0,k)
if(is.null(c2))c2=rep(0,k)
c1 <- c1[1:k]
c2 <- c2[1:k]
p1 <- ncol(phi1)%/%k
p2 <- ncol(phi2)%/%k
###
"mtxroot" <- function(sigma){
if(!is.matrix(sigma))sigma=as.matrix(sigma)
sigma <- (t(sigma)+sigma)/2
m1 <- eigen(sigma)
P <- m1$vectors
L <- diag(sqrt(m1$values))
sigmah <- P%*%L%*%t(P)
sigmah
}
###
s1 <- mtxroot(sigma1[1:k,1:k])
s2 <- mtxroot(sigma2[1:k,1:k])
nT <- nob+ini
et <- matrix(rnorm(k*nT),nT,k)
a1 <- et%*%s1
a2 <- et%*%s2
###
d <- delay[2]
if(d <= 0)d <- 1
p <- max(p1,p2,d)
### set starting values
zt <- matrix(1,p,1)%*%matrix(c1,1,k)+a1[1:p,]
resi <- matrix(0,nT,k)
ist = p+1
for (i in ist:ini){
 wk = rep(0,k)
 if(zt[(i-d),delay[1]] <= thr){
 resi[i,] <- a1[i,]
 wk <- wk + a1[i,] + c1
 if(p1 > 0){
  for (j in 1:p1){
   idx <- (j-1)*k
   phi <- phi1[,(idx+1):(idx+k)]
   wk <- wk + phi%*%matrix(zt[i-j,],k,1)
   }
  }
 }else{ resi[i,] <- a2[i,]
     wk <- wk+a2[i,] + c2
    if(p2 > 0){
     for (j in 1:p2){
      idx <- (j-1)*k
      phi <- phi2[,(idx+1):(idx+k)]
      wk <- wk + phi%*%matrix(zt[i-j,],k,1)
      }
     }
   }
 zt <- rbind(zt,c(wk))
 }
#### generate data used
n1 = 0
for (i in (ini+1):nT){
 wk = rep(0,k)
 if(zt[(i-d),delay[1]] <= thr){
 n1 = n1+1
 resi[i,] <- a1[i,]
 wk <- wk + a1[i,] + c1
 if(p1 > 0){
  for (j in 1:p1){
   idx <- (j-1)*k
   phi <- phi1[,(idx+1):(idx+k)]
   wk <- wk + phi%*%matrix(zt[i-j,],k,1)
   }
  }
 }else{ resi[i,] <- a2[i,]
     wk <- wk+a2[i,] + c2
    if(p2 > 0){
     for (j in 1:p2){
      idx <- (j-1)*k
      phi <- phi2[,(idx+1):(idx+k)]
      wk <- wk + phi%*%matrix(zt[i-j,],k,1)
      }
     }
   }
 zt <- rbind(zt,c(wk))
 }
mTAR.sim <- list(series = zt[(ini+1):nT,],at = resi[(ini+1):nT,],threshold=thr,
delay=delay, n1=n1, n2 = (nob-n1))
}


#' Estimation of a Multivariate Two-Regime SETAR Model
#'
#' Estimation of a multivariate two-regime SETAR model, including threshold.
#' The procedure of Li and Tong (2016) is used to search for the threshold.
#' @param y a (nT-by-k) data matrix of multivariate time series, where nT is the sample size and k is the dimension.
#' @param p1 AR-order of regime 1.
#' @param p2 AR-order of regime 2.
#' @param d delay for threshold variable, default value is 1.
#' @param thr threshold variable. Estimation is needed if thr = NULL.
#' @param delay two elements (i,d) with "i" being the component and "d" the delay for threshold variable.
#' @param thrV vector of threshold variable. If it is not null, thrV must have the same sample size of that of y.
#' @param Trim lower and upper quantiles for possible threshold value.
#' @param k0 the maximum number of threshold values to be evaluated.
#' @param include.mean logical values indicating whether constant terms are included.
#' @param score the choice of criterion used in selection threshold, namely (AIC, det(RSS)).
#' @references
#' Li, D., and Tong. H. (2016) Nested sub-sample search algorithm for estimation of threshold models. \emph{Statisitca Sinica}, 1543-1554.
#' @return mTAR returns a listh with the following components:
#' \item{data}{the data matrix, y.}
#' \item{beta}{a (p*k+1)-by-(2k) matrices. The first k columns show the estimation results in regime 1, and the second k columns show these in regime 2.}
#' \item{arorder}{AR orders of regimes 1 and 2.}
#' \item{sigma}{estimated innovational covariance matrices of regimes 1 and 2.}
#' \item{residuals}{estimated innovations.}
#' \item{nobs}{numbers of observations in regimes 1 and 2.}
#' \item{model1, model2}{estimated models of regimes 1 and 2.}
#' \item{thr}{threshold value.}
#' \item{delay}{two elements (i,d) with "i" being the component and "d" the delay for threshold variable.}
#' \item{thrV}{vector of threshold variable.}
#' \item{D}{D.}
#' \item{RSS}{RSS.}
#' \item{information}{overall information criteria.}
#' \item{cnst}{logical values indicating whether the constant terms are included in regimes 1 and 2.}
#' \item{sresi}{standardized residuals.} 
#' @examples
#' phi1=matrix(c(0.5,0.7,0.3,0.2),2,2)
#' phi2=matrix(c(0.4,0.6,0.5,-0.5),2,2)
#' sigma1=matrix(c(1,0,0,1),2,2)
#' sigma2=matrix(c(1,0,0,1),2,2)
#' c1=c(0,0)
#' c2=c(0,0)
#' delay=c(1,1)
#' y=mTAR.sim(100,0,phi1,phi2,sigma1,sigma2,c1,c2,delay,ini=500)
#' est=mTAR(y$series,1,1,0,y$series,delay,Trim,300,include.mean,"AIC")
#' @export 
"mTAR" <- function(y,p1,p2,thr=NULL,thrV=NULL,delay=c(1,1),Trim=c(0.15,0.85),k0=300,include.mean=TRUE,score="AIC"){
if(!is.matrix(y))y <- as.matrix(y)
p = max(p1,p2)
if(p < 1)p=1
d = delay[2]
ist = max(p,d)+1
nT <- nrow(y)
ky <- ncol(y)
nobe <- nT-ist+1
if(length(thrV) < nT){
         thrV <- y[(ist-d):(nT-d),delay[1]]
	 }else{
          thrV <- thrV[ist:nT]
	  }
beta <- matrix(0,(p*ky+1),ky*2)
### set up regression framework
yt <- y[ist:nT,]
x1 <- NULL
for (i in 1:p){
 x1=cbind(x1,y[(ist-i):(nT-i),])
 }
X <- x1
if(include.mean){X=cbind(rep(1,nobe),X)}
### search for threshold if needed
if(is.null(thr)){
 D <- NeSS(yt,x1,thrV=thrV,Trim=Trim,k0=k0,include.mean=include.mean,score=score)
##cat("Threshold candidates: ",D,"\n")
 ll = length(D)
 RSS=NULL
 ic = 2
#### If not AIC, generalized variance is used.
 if(score=="AIC")ic <- 1
 for (i in 1:ll){
  thr = D[i]
  idx = c(1:nobe)[thrV <= thr]
  m1a <- MlmNeSS(yt,X,subset=idx,SD=FALSE,include.mean=include.mean)
  m1b <- MlmNeSS(yt,X,subset=-idx,SD=FALSE,include.mean=include.mean)
   RSS=c(RSS,sum(m1a$score[ic],m1b$score[ic]))
  }
#cat("RSS: ",RSS,"\n")
 j=which.min(RSS)
 thr <- D[j]
 cat("Estimated Threshold: ",thr,"\n")
 }else{
 cat("Threshold: ",thr,"\n")
 RSS <- NULL
}
#### The above ends the search of threshold ###
### Final estimation results
###cat("Coefficients are in vec(beta), where t(beta)=[phi0,phi1,...]","\n")
  resi = matrix(0,nobe,ncol(yt))
  sresi <- matrix(0,nobe,ncol(yt))
  idx = c(1:nobe)[thrV <= thr]
  n1 <- length(idx)
  n2 <- nrow(yt)-n1
  k <- ncol(yt)
  X = x1
  if(include.mean){X=cbind(rep(1,nobe),X)}
  m1a <- MlmNeSS(yt,X,subset=idx,SD=TRUE,include.mean=include.mean)
  np <- ncol(X)
  resi[idx,] <- m1a$residuals
  infc = m1a$information
  beta[1:nrow(m1a$beta),1:ky] <- m1a$beta
  cat("Regime 1 with sample size: ",n1,"\n")
  coef1 <- m1a$coef
  coef1 <- cbind(coef1,coef1[,1]/coef1[,2])
  colnames(coef1) <- c("est","s.e.","t-ratio")
  for (j in 1:k){
   jdx=(j-1)*np
   cat("Model for the ",j,"-th component (including constant, if any): ","\n")
   print(round(coef1[(jdx+1):(jdx+np),],4))
   }
 cat("sigma: ", "\n")
 print(m1a$sigma)
 cat("Information(aic,bix,hq): ",infc,"\n")
 m1 <- eigen(m1a$sigma)
 P <- m1$vectors
 Di <- diag(1/sqrt(m1$values))
 siv <- P%*%Di%*%t(P)
 sresi[idx,] <- m1a$residuals%*%siv
 ###
  m1b <- MlmNeSS(yt,X,subset=-idx,SD=TRUE,include.mean=include.mean)
  resi[-idx,] <- m1b$residuals
  beta[1:nrow(m1b$beta),(ky+1):(2*ky)] <- m1b$beta
  cat(" ","\n")
  cat("Regime 2 with sample size: ",n2,"\n")
  coef2 <- m1b$coef
  coef2 <- cbind(coef2,coef2[,1]/coef2[,2])
  colnames(coef2) <- c("est","s.e.","t-ratio")
  for (j in 1:k){
   jdx=(j-1)*np
   cat("Model for the ",j,"-th component (including constant, if any): ","\n")
   print(round(coef2[(jdx+1):(jdx+np),],4))
   }  
  cat("sigma: ","\n")
  print(m1b$sigma)
  cat("Information(aic,bic,hq): ",m1b$information,"\n")
  m2 <- eigen(m1b$sigma)
  P <- m2$vectors
  Di <- diag(1/sqrt(m2$values))
  siv <- P%*%Di%*%t(P)
  sresi[-idx,] <- m1b$residuals%*%siv
#
 fitted.values <- yt-resi
 cat(" ","\n")
 cat("Overall pooled estimate of sigma: ","\n")
 sigma = (n1*m1a$sigma+n2*m1b$sigma)/(n1+n2)
 print(sigma)
 infc = m1a$information+m1b$information
 cat("Overall information criteria(aic,bic,hq): ",infc,"\n")
 Sigma <- cbind(m1a$sigma,m1b$sigma)
mTAR <- list(data=y,beta=beta,arorder=c(p1,p2),sigma=Sigma,residuals = resi,nobs=c(n1,n2),
model1 = m1a, model2 = m1b,thr=thr, delay=delay, thrV=thrV, D=D, RSS=RSS, information=infc,
cnst=rep(include.mean,2),sresi=sresi)
}



#' Estimation of Multivariate TAR Models
#'
#' Estimation of mutlivariate TAR models with given thresholds. It can handle multiple regimes.
#' @param y vector time series
#' @param arorder AR order of each regime. The number of regime is length of arorder.
#' @param thr threshould value(s). There are k-1 threshold for a k-regime model.
#' @param delay two elements (i,d) with "i" being the component and "d" the delay for threshold variable.
#' @param thrV external threhold variable if any. If thrV is not null, it must have the same number of observations as y-series.
#' @param include.mean logical values indicating whether constant terms are included. Default is TRUE for all.
#' @param output a logical value indicating four output. Default is TRUE.
#' @return mTAR.est returns a list with the following components:
#' \item{data}{the data matrix, y.}
#' \item{k}{the dimension of y.}
#' \item{arorder}{AR orders of regimes 1 and 2.}
#' \item{beta}{a (p*k+1)-by-(2k) matrices. The first k columns show the estimation results in regime 1, and the second k columns show these in regime 2.}
#' \item{sigma}{estimated innovational covariance matrices of regimes 1 and 2.}
#' \item{thr}{threshold value.}
#' \item{residuals}{estimated innovations.}
#' \item{sresi}{standardized residuals.}
#' \item{nobs}{numbers of observations in different regimes.}
#' \item{cnst}{logical values indicating whether the constant terms are included in different regimes.}
#' \item{AIC}{AIC value.}
#' \item{delay}{two elements (i,d) with "i" being the component and "d" the delay for threshold variable.}
#' \item{thrV}{values of threshold variable.}
#' @examples
#' phi1=matrix(c(0.5,0.7,0.3,0.2),2,2)
#' phi2=matrix(c(0.4,0.6,0.5,-0.5),2,2)
#' sigma1=matrix(c(1,0,0,1),2,2)
#' sigma2=matrix(c(1,0,0,1),2,2)
#' c1=c(0,0)
#' c2=c(0,0)
#' delay=c(1,1)
#' y=mTAR.sim(100,0,phi1,phi2,sigma1,sigma2,c1,c2,delay,ini=500)
#' est=mTAR.est(y$series,c(1,1),0,delay)
#' @export
"mTAR.est" <- function(y,arorder=c(1,1),thr=c(0),delay=c(1,1),thrV=NULL,include.mean=c(TRUE,TRUE),output=TRUE){
if(!is.matrix(y)) y <- as.matrix(y)
p <- max(arorder)
if(p < 1)p <-1
k <- length(arorder)
d <- delay[2]
ist <- max(p,d)+1
nT <- nrow(y)
ky <- ncol(y)
Y <- y[ist:nT,]
X <- NULL
for (i in 1:p){
    X <- cbind(X,y[(ist-i):(nT-i),])
    }
X <- as.matrix(X)
nobe <- nT-ist+1
if(length(thrV) < nT){
 thrV <- y[(ist-d):(nT-d),delay[1]]
 if(output) cat("Estimation of a multivariate SETAR model: ","\n")
 }else{thrV <- thrV[ist:nT]}
### make sure that the threholds are in increasing order!!!
thr <- sort(thr)
k1 <- length(thr)
### house keeping
 beta <- NULL
 sigma <- NULL
 nobs <- NULL
 resi <- matrix(0,nobe,ky)
 sresi <- matrix(0,nobe,ky)
 npar <- NULL
 AIC <- 99999
#
if(k1 > (k-1)){
 cat("Only the first: ",k1," thresholds are used.","\n")
 THr <- min(thrV)-1
 thr <- thr[1:(k-1)]
 THr <- c(THr,thr,(max(thrV)+1))
 }
if(k1 < (k-1)){cat("Too few threshold values are given!!!","\n")}
if(length(include.mean) != k)include.mean=rep("TRUE",k)
if(k1 == (k-1)){
 THr <- c((min(thrV[1:nobe])-1),thr,c(max(thrV[1:nobe])+1))
 if(output) cat("Thresholds used: ",thr,"\n")
 for (j in 1:k){
  idx <- c(1:nobe)[thrV <= THr[j]]
  jdx <- c(1:nobe)[thrV > THr[j+1]]
  jrm <- c(idx,jdx)
  n1 <- nobe-length(jrm)
  nobs <- c(nobs,n1)
  kdx <- c(1:nobe)[-jrm]
  x1 <- X
  if(include.mean[j]){x1 <- cbind(rep(1,nobe),x1)}
  np <- ncol(x1)*ky
  npar <- c(npar,np)
  beta <- matrix(0,(p*ky+1),k*ky)
  sig <- diag(rep(1,ky))
  if(n1 <= ncol(x1)){
   cat("Regime ",j," does not have sufficient observations!","\n")
   sigma <- cbind(sigma,sig)
   }
   else{
   mm <- MlmNeSS(Y,x1,subset=kdx,SD=TRUE,include.mean=include.mean[j])
   jst = (j-1)*k
   beta1 <- mm$beta
   beta[1:nrow(beta1),(jst+1):(jst+k)] <- beta1
   Sigma <- mm$sigma
   if(output){
    cat("Estimation of regime: ",j," with sample size: ",n1,"\n")   
    cat("Coefficient estimates:t(phi0,phi1,...) ","\n")
    colnames(mm$coef) <- c("est","se")
    print(mm$coef)
    cat("Sigma: ","\n")
    print(Sigma)
    }
   sigma <- cbind(sigma,Sigma)
   resi[kdx,] <- mm$residuals
   m1 <- eigen(Sigma)
   P <- m1$vectors
   Di <- diag(1/sqrt(m1$values))
   Sighinv <- P%*%Di%*%t(P)
   sresi[kdx,] <- mm$residuals %*% Sighinv
   }
  }
  AIC <- 0
  for (j in 1:k){
   sig <- sigma[,((j-1)*k+1):(j*k)]
   d1 <- det(sig)
   n1 <- nobs[j]
   np <- npar[j]
   if(!is.null(sig)){s2 <- sig^2*(n1-np)/n1
               AIC <- AIC + log(d1)*n1 + 2*np
	       }
	}
  if(output) cat("Overal AIC: ",AIC,"\n")
}
mTAR.est <- list(data=y,k = k, arorder=arorder,beta=beta,sigma=sigma,thr=thr,residuals=resi,
sresi=sresi,nobs=nobs,cnst=include.mean, AIC=AIC, delay=delay,thrV=thrV)
}



#' Refine a Fitted 2-Regime Multivariate TAR Models 
#'
#' Refine a fitted 2-regime multivariate TAR models using "thres" as threshold for t-ratios.
#' @param m1 A fitted mTAR object.
#' @param thres The threshold value.
#' @return ref.mTAR returns a list with following components:
#' \item{data}{data matrix, y.}
#' \item{arorder}{AR orders of regimes 1 and 2.}
#' \item{sigma}{estimated innovational covariance matrices of regimes 1 and 2.}
#' \item{beta}{a (p*k+1)-by-(2k) matrices. The first k columns show the estimation results in regime 1, and the second k columns shows these in regime 2.}
#' \item{residuals}{estimated innovations.}
#' \item{sresi}{standard residuals.}
#' \item{criteria}{overall information criteria.}
#' @examples
#' phi1=matrix(c(0.5,0.7,0.3,0.2),2,2)
#' phi2=matrix(c(0.4,0.6,0.5,-0.5),2,2)
#' sigma1=matrix(c(1,0,0,1),2,2)
#' sigma2=matrix(c(1,0,0,1),2,2)
#' c1=c(0,0)
#' c2=c(0,0)
#' delay=c(1,1)
#' y=mTAR.sim(100,0,phi1,phi2,sigma1,sigma2,c1,c2,delay,ini=500)
#' est=mTAR.est(y$series,c(1,1),0,delay)
#' ref.mTAR(est,0)
#' @export
"ref.mTAR" <- function(m1,thres=1.0){
y <- m1$data
arorder <- m1$arorder
sigma <- m1$sigma
nobs <- m1$nobs
thr <- m1$thr
delay <- m1$delay
thrV <- m1$thrV
include.mean <- m1$cnst
mA <- m1$model1
mB <- m1$model2
vName <- colnames(y)
if(is.null(vName)){
        vName <- paste("V",1:ncol(y))
        colnames(y) <- vName
     }
###
p = max(arorder)
d = delay[2]
ist = max(p,d)+1
nT <- nrow(y)
ky <- ncol(y)
nobe <- nT-ist+1
### set up regression framework
yt <- y[ist:nT,]
x1 <- NULL
xname <- NULL
vNameL <- paste(vName,"L",sep="")
for (i in 1:p){
 x1=cbind(x1,y[(ist-i):(nT-i),])
 xname <- c(xname,paste(vNameL,i,sep=""))
 }
colnames(x1) <- xname
X <- x1
if(include.mean[1]){X <- cbind(rep(1,nobe),X)
            xname <- c("cnst",xname)
          }
colnames(X) <- xname
####
### Final estimation results
###cat("Coefficients are in vec(beta), where t(beta)=[phi0,phi1,...]","\n")
 npar <- NULL
  np <- ncol(X)
  beta <- matrix(0,np,ky*2)
  resi = matrix(0,nobe,ky)
  sresi <- matrix(0,nobe,ky)
  idx = c(1:nobe)[thrV <= thr]
  n1 <- length(idx)
  n2 <- nrow(yt)-n1
### Regime 1, component-by-component
  cat("Regime 1 with sample size: ",n1,"\n")
  coef <- mA$coef
  RES <- NULL
  for (ij in 1:ky){
   tra <- rep(0,np)
   i1 <- (ij-1)*np
   jj <- c(1:np)[coef[(i1+1):(i1+np),2] > 0]
   tra[jj] <- coef[(i1+jj),1]/coef[(i1+jj),2]
###   cat("tra: ",tra,"\n")
   kdx <- c(1:np)[abs(tra) >= thres]
  npar <- c(npar,length(kdx))
##
  cat("Significant variables: ",kdx,"\n")
  if(length(kdx) > 0){ 
   xreg <- as.matrix(X[idx,kdx])
   colnames(xreg) <- colnames(X)[kdx]
   m1a <- lm(yt[idx,ij]~-1+xreg)
   m1aS <- summary(m1a)
   beta[kdx,ij] <- m1aS$coefficients[,1]
   resi[idx,ij] <- m1a$residuals
   cat("Component: ",ij,"\n")
   print(m1aS$coefficients)
   }else{
     resi[idx,ij] <- yt[idx,ij]
   }
  RES <- cbind(RES,c(resi[idx,ij]))
  }
###
 sigma1 <- crossprod(RES,RES)/n1
 cat("sigma: ", "\n")
 print(sigma1)
# cat("Information(aic,bix,hq): ",infc,"\n")
 m1 <- eigen(sigma1)
 P <- m1$vectors
 Di <- diag(1/sqrt(m1$values))
 siv <- P%*%Di%*%t(P)
 sresi[idx,] <- RES%*%siv
 ###  Regime 2 ####
  cat("Regime 2 with sample size: ",n2,"\n")
  coef <- mB$coef
  RES <- NULL
  for (ij in 1:ky){
   tra <- rep(0,np)
   i1 <- (ij-1)*np
   jj <- c(1:np)[coef[(i1+1):(i1+np),2] > 0]
   tra[jj] <- coef[(i1+jj),1]/coef[(i1+jj),2]
###   cat("tra: ",tra,"\n")
   kdx <- c(1:np)[abs(tra) >= thres]
##
  cat("Significant variables: ",kdx,"\n")
  npar <- c(npar,length(kdx))
  if(length(kdx) > 0){
   xreg <- as.matrix(X[-idx,kdx])
   colnames(xreg) <- colnames(X)[kdx]
   m1a <- lm(yt[-idx,ij]~-1+xreg)
   m1aS <- summary(m1a)
   beta[kdx,ij+ky] <- m1aS$coefficients[,1]
   resi[-idx,ij] <- m1a$residuals
   cat("Component: ",ij,"\n")
   print(m1aS$coefficients)
   }else{
    resi[-idx,ij] <- yt[-idx,ij]
    }
  RES <- cbind(RES,c(resi[-idx,ij]))
  }
###
 sigma2 <- crossprod(RES,RES)/n2
 cat("sigma: ", "\n")
 print(sigma2)
# cat("Information(aic,bix,hq): ",infc,"\n")
 m1 <- eigen(sigma2)
 P <- m1$vectors
 Di <- diag(1/sqrt(m1$values))
 siv <- P%*%Di%*%t(P)
 sresi[-idx,] <- RES%*%siv
#
 fitted.values <- yt-resi
 cat(" ","\n")
 cat("Overall pooled estimate of sigma: ","\n")
 sigma = (n1*sigma1+n2*sigma2)/(n1+n2)
 print(sigma)
 d1 <- log(det(sigma))
 nnp <- sum(npar)
 TT <- (n1+n2)
 aic <- TT*d1+2*nnp
 bic <- TT*d1+log(TT)*nnp
 hq <- TT*d1+2*log(log(TT))*nnp
 cri <- c(aic,bic,hq)
### infc = m1a$information+m1b$information
cat("Overall information criteria(aic,bic,hq): ",cri,"\n")
 Sigma <- cbind(sigma1,sigma2)
ref.mTAR <- list(data=y,arorder=arorder,sigma=Sigma,beta=beta,residuals=resi,sresi=sresi,criteria=cri)
}


#' Prediction of a Fitted Multivariate TAR Model
#' @param model multivariate TAR model.
#' @param orig forecast origin.
#' @param h forecast horizon.
#' @param iternations number of iterations.
#' @param ci confidence level.
#' @param output a logical value for output.
#' @return mTAR.pred returns a list with components:
#' \item{model}{the multivariate TAR model.}
#' \item{pred}{prediction.}
#' \item{Ysim}{fitted y.}
#' @examples
#' phi1=matrix(c(0.5,0.7,0.3,0.2),2,2)
#' phi2=matrix(c(0.4,0.6,0.5,-0.5),2,2)
#' sigma1=matrix(c(1,0,0,1),2,2)
#' sigma2=matrix(c(1,0,0,1),2,2)
#' c1=c(0,0)
#' c2=c(0,0)
#' delay=c(1,1)
#' y=mTAR.sim(100,0,phi1,phi2,sigma1,sigma2,c1,c2,delay,ini=500)
#' est=mTAR.est(y$series,c(1,1),0,delay)
#' pred=mTAR.pred(est,100,1,300,0.90,TRUE)
#' @export
"mTAR.pred" <- function(model,orig,h=1,iterations=3000,ci=0.95,output=TRUE){
## ci: probability of pointwise confidence interval coefficient
### only works for self-exciting mTAR models.
###
y <- model$data
arorder <- model$arorder
beta <- model$beta
sigma <- model$sigma
thr <- model$thr
include.mean <- model$cnst
delay <- model$delay
p <- max(arorder)
k <- length(arorder)
d <- delay[2]
nT <- nrow(y)
ky <- ncol(y)
if(orig < 1)orig <- nT
if(orig > nT) orig <- nT
if(h < 1) h <- 1
### compute the predictions. Use simulation if the threshold is predicted
Sigh <- NULL
for (j in 1:k){
 sig <- sigma[,((j-1)*ky+1):(j*ky)]
 m1 <- eigen(sig)
 P <- m1$vectors
 Di <- diag(sqrt(m1$values))
 sigh <- P%*%Di%*%t(P)
 Sigh <- cbind(Sigh,sigh)
 }
Ysim <- array(0,dim=c(h,ky,iterations))
for (it in 1:iterations){
  yp <- y[1:orig,]
  et <- matrix(rnorm(h*ky),h,ky)
  for (ii in 1:h){
   t <- orig+ii
   thd <- yp[(t-d),delay[1]]
   JJ <- 1
   for (j in 1:(k-1)){
     if(thd > thr[j]){JJ <- j+1}
     }
   Jst <- (JJ-1)*ky
   at <- matrix(et[ii,],1,ky)%*%Sigh[,(Jst+1):(Jst+ky)]
   x <- NULL
   if(include.mean[JJ]) x <- 1
   pJ <- arorder[JJ]
   phi <- beta[,(Jst+1):(Jst+ky)]
   for (i in 1:pJ){
    x <- c(x,yp[(t-i),])
    }
   yhat <- matrix(x,1,length(x))%*%phi[1:length(x),]
   yhat <- yhat+at
   yp <- rbind(yp,yhat)
   Ysim[ii,,it] <- yhat
   }
 }
### summary
 pred <- NULL
 upp <- NULL
 low <- NULL
 pr <- (1-ci)/2
 pro <- c(pr, 1-pr)
 for (ii in 1:h){
  fst <- NULL
  lowb <- NULL
  uppb <- NULL
  for (j in 1:ky){
   ave <- mean(Ysim[ii,j,])
   quti <- quantile(Ysim[ii,j,],prob=pro)
   fst <- c(fst,ave)
   lowb <- c(lowb,quti[1])
   uppb <- c(uppb,quti[2])
   }
   pred <- rbind(pred,fst)
   low <- rbind(low,lowb)
   upp <- rbind(upp,uppb)
 }
if(output){
 colnames(pred) <- colnames(y)
 cat("Forecast origin: ",orig,"\n")
 cat("Predictions: 1-step to ",h,"-step","\n")
 print(pred)
 cat("Lower bounds of ",ci*100," % confident intervals","\n")
 print(low)
 cat("Upper bounds: ","\n")
 print(upp)
}
mTAR.pred <- list(data=y,pred = pred,Ysim=Ysim)
}


#' This program performs multivariate linear regression analysis.
#' @export
"MlmNeSS" <- function(y,z,subset=NULL,SD=FALSE,include.mean=TRUE){
# z: design matrix, including 1 as its first column if constant is needed.
# y: dependent variables
# subset: locators for the rows used in estimation
## Model is y = z%*%beta+error
#### include.mean is ONLY used in counting parameters. z-matrix should include column of 1 if constant
#### is needed.
#### SD: switch for computing standard errors of the estimates
####
zc = as.matrix(z)
yc = as.matrix(y)
if(!is.null(subset)){
 zc <- zc[subset,]
 yc <- yc[subset,]
 }
n=nrow(zc)
p=ncol(zc)
k <- ncol(y)
coef=NULL
infc = rep(9999,3)
if(n <= p){
  cat("too few data points in a regime: ","\n")
   beta = NULL
   res = yc
   sig = NULL
  }
  else{
   ztz=t(zc)%*%zc/n
   zty=t(zc)%*%yc/n
   ztzinv=solve(ztz)
   beta=ztzinv%*%zty
   res=yc-zc%*%beta
### Sum of squares of residuals
   sig=t(res)%*%res/n
   npar=ncol(zc)
   if(include.mean)npar=npar-1
   score1=n*log(det(sig))+2*npar
   score2=det(sig*n)
   score=c(score1,score2)
###
   if(SD){
    sigls <- sig*(n/(n-p))
    s <- kronecker(sigls,ztzinv)
    se <- sqrt(diag(s)/n)
    coef <- cbind(c(beta),se)
#### MLE estimate of sigma 
    d1 = log(det(sig))
    npar=nrow(coef)
    if(include.mean)npar=npar-1
    aic <- n*d1+2*npar
    bic <- n*d1 + log(n)*npar
    hq <- n*d1 + 2*log(log(n))*npar
    infc=c(aic,bic,hq)
   }
 }
#
MlmNeSS <- list(beta=beta,residuals=res,sigma=sig,coef=coef,information=infc,score=score)
}

