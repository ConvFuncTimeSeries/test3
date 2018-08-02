#' Generate Univariate SETAR Models
#'
#' Generate univariate SETAR model for up to 3 regimes.
#' @param nob number of observations.
#' @param arorder AR-order for each regime. The length of arorder controls the number of regimes.
#' @param phi a 3-by-p matrix. Each row contains the AR coefficients for a regime.
#' @param d delay for threshold variable.
#' @param thr threshold values. 
#' @param sigma standard error for each regime.
#' @param cnst constant terms.
#' @param ini burn-in period.
#' @return uTAR.sim returns a list with components:
#' \item{series}{a time series following SETAR model.}
#' \item{at}{innovation of the time seres.}
#' \item{arorder}{AR-order for each regime.}
#' \item{thr}{threshold value.}
#' \item{phi}{a 3-by-p matrix. Each row contains the AR coefficients for a regime.}
#' \item{cnst}{constant terms}
#' \item{sigma}{standard error for each regime.}
#' @examples
#' arorder=rep(1,2)
#' ar.coef=matrix(c(0.7,-0.8),2,1)
#' y=uTAR.sim(100,arorder,ar.coef,1,0)
#' @export
"uTAR.sim" <- function(nob,arorder,phi,d=1,thr=c(0,0),sigma=c(1,1,1),cnst=rep(0,3),ini=500){
p <- max(arorder)
nT <- nob+ini
ist <- max(p,d)+1
et <- rnorm(nT)
nregime <- length(arorder)
if(nregime > 3)nregime <- 3
## For two-regime models, regime 3 is set to the same as regime 2.
if(nregime == 2)arorder=c(arorder[1],arorder[2],arorder[2])
k1 <- length(cnst)
if(k1 < 3)cnst <- c(cnst,rep(cnst[k1],(3-k1)))
k1 <- length(sigma)
if(k1 < 3) sigma <- c(sigma,rep(sigma[k1],(3-k1)))
if(!is.matrix(phi))phi <- as.matrix(phi)
k1 <- nrow(phi)
if(k1 < 3){
  for (j in (k1+1):3){
    phi <- rbind(phi,phi[k1,])
    }
 }
k1 <- length(thr)
if(k1 < 2) thr=c(thr,10^7)
if(nregime == 2)thr[2] <- 10^7
zt <- et[1:ist]*sigma[1]
at <- zt
##
for (i in ist:nT){
  if(zt[i-d] <= thr[1]){
    at[i] = et[i]*sigma[1]
    tmp = cnst[1]
    if(arorder[1] > 0){
      for (j in 1:arorder[1]){
       tmp <- tmp + zt[i-j]*phi[1,j]
       }
     }
    zt[i] <- tmp + at[i]
    }
    else{ if(zt[i-d] <= thr[2]){
           at[i] <- et[i]*sigma[2]
	   tmp <- cnst[2]
	   if(arorder[2] > 0){
	    for (j in 1:arorder[2]){
	     tmp <- tmp + phi[2,j]*zt[i-j]
	     }
	    }
	   zt[i] <- tmp + at[i]
	   }
	   else{ at[i] <- et[i]*sigma[3]
	         tmp <- cnst[3]
	         if(arorder[3] > 0){
		  for (j in 1:arorder[3]){
		    tmp <- tmp + phi[3,j]*zt[i-j]
		    }
		   }
		 zt[i] <- tmp + at[i]
		}
	 }

   }	  
uTAR.sim <- list(series = zt[(ini+1):nT], at=at[(ini+1):nT], arorder=arorder,thr=thr,phi=phi,cnst=cnst,sigma=sigma)
}



#' Estimation of a Univariate Two-Regime SETAR Model
#'
#' Estimation of a univariate two-regime SETAR model, including threshold value.
#' The procedure of Li and Tong (2016) is used to search for the threshold.
#' @param y a vector of time series.
#' @param p1,p2 AR-orders of regime 1 and regime 2.
#' @param d delay for threshold variable, default is 1.
#' @param thrV threshold variable. If thrV is not null, it must have the same length as that of y.
#' @param Trim lower and upper quantiles for possible threshold values.
#' @param k0 the maximum number of threshold values to be evaluated. If the sample size is large (> 3000), then k0 = floor(nT*0.5).The default is k0=300. But k0 = floor(nT*0.8) if nT < 300.
#' @param include.mean a logical value indicating whether constant terms are included.
#' @param thrQ lower and upper quantiles to search for threshold value.
#' @references
#' Li, D., and Tong. H. (2016) Nested sub-sample search algorithm for estimation of threshold models. \emph{Statisitca Sinica}, 1543-1554.
#' @return uTAR returns a list with components:
#' \item{data}{the data matrix, y.}
#' \item{arorder}{AR orders of regimes 1 and 2.}
#' \item{delay}{the delay for threshold variable.}
#' \item{residuals}{estimated innovations.}
#' \item{coef}{a 2-by-(p+1) matrices. The first row show the estimation results in regime 1, and the second row shows these in regime 2.}
#' \item{sigma}{estimated innovational covariance matrices of regimes 1 and 2.}
#' \item{nobs}{numbers of observations in regimes 1 and 2.}
#' \item{model1,model2}{estimated models of regimes 1 and 2.}
#' \item{thr}{threshold value.}
#' \item{D}{a set of threshold values.} 
#' \item{RSS}{RSS}
#' \item{AIC}{AIC value}
#' \item{cnst}{logical values indicating whether the constant terms are included in regimes 1 and 2.}
#' \item{sresi}{standardized residuals.}
#' @examples
#' arorder=rep(1,2)
#' ar.coef=matrix(c(0.7,-0.8),2,1)
#' y=uTAR.sim(100,arorder,ar.coef,1,0)
#' est=uTAR(y$series,1,1,1,y$series,c(0.2,0.8),100,TRUE,c(0.2,0.8))
#' @export
"uTAR" <- function(y,p1,p2,d=1,thrV=NULL,Trim=c(0.15,0.85),k0=300,include.mean=T,thrQ=c(0,1)){
if(is.matrix(y))y=y[,1]
p = max(p1,p2)
if(p < 1)p=1
ist=max(p,d)+1
nT <- length(y)
### regression framework
if(k0 > nT)k0 <- floor(nT*0.8)
if(nT > 3000)k0 <- floor(nT*0.5)
yt <- y[ist:nT]
nobe <- nT-ist+1
x1 <- NULL
for (i in 1:p){
 x1=cbind(x1,y[(ist-i):(nT-i)])
 }
if(length(thrV) < nT){
    tV <- y[(ist-d):(nT-d)]
    }else{
    tV <- thrV[ist:nT]
    }
D <- NeSS(yt,x1,x2=NULL,thrV=tV,Trim=Trim,k0=k0,include.mean=include.mean,thrQ=thrQ)
#cat("Threshold candidates: ",D,"\n")
ll = length(D)
### house keeping
RSS=NULL
resi <- NULL
sresi <- NULL
m1a <- NULL
m1b <- NULL
Size <- NULL
Phi <- matrix(0,2,p+1)
Sigma <- NULL
infc <- NULL
thr <- NULL
##
if(ll > 0){
for (i in 1:ll){
  thr = D[i]
  idx = c(1:nobe)[tV <= thr]
  X = x1[,c(1:p1)]
  if(include.mean){X=cbind(rep(1,nobe),X)}
  m1a <- lm(yt~-1+.,data=data.frame(X),subset=idx)
  X <- x1[,c(1:p2)]
  if(include.mean){X=cbind(rep(1,nobe),X)}
  m1b <- lm(yt~-1+.,data=data.frame(X),subset=-idx)
  RSS=c(RSS,sum(m1a$residuals^2,m1b$residuals^2))
  }
##cat("RSS: ",RSS,"\n")
j=which.min(RSS)
thr <- D[j]
cat("Estimated Threshold: ",thr,"\n")
### Final estimation results
  resi = rep(0,nobe)
  sresi <- rep(0,nobe)
  idx = c(1:nobe)[tV <= thr]
  n1 <- length(idx)
  n2 <- length(yt)-n1
  Size <- c(n1,n2)
  X = x1[,c(1:p1)]
  if(include.mean){X=cbind(rep(1,(nT-ist+1)),X)}
  m1a <- lm(yt~-1+.,data=data.frame(X),subset=idx)
  resi[idx] <- m1a$residuals
 X <- as.matrix(X)
 cat("Regime 1: ","\n")
 Phi[1,1:ncol(X)] <- m1a$coefficients
 coef1 <- summary(m1a)
 sigma1 <- coef1$sigma
 sresi[idx] <- m1a$residuals/sigma1
 print(coef1$coefficients)
 sig1 <- sigma1^2*(n1-ncol(X))/n1
 AIC1 <- n1*log(sig1) + 2*ncol(X) 
 cat("nob1, sigma1 and AIC: ",c(n1, sigma1, AIC1), "\n")
  X <- x1[,c(1:p2)]
  if(include.mean){X=cbind(rep(1,nobe),X)}
  m1b <- lm(yt~-1+.,data=data.frame(X),subset=-idx)
  resi[-idx] <- m1b$residuals
  X <- as.matrix(X)
  cat("Regime 2: ","\n")
  Phi[2,1:ncol(X)] <- m1b$coefficients
  coef2 <- summary(m1b)
  sigma2 <- coef2$sigma
  sresi[-idx] <- m1b$residuals/sigma2
  print(coef2$coefficients)
  sig2 <- sigma2^2*(n2-ncol(X))/n2
  AIC2 <- n2*log(sig2) + 2*ncol(X)
  cat("nob2, sigma2 and AIC: ",c(n2,sigma2, AIC2),"\n")
 fitted.values <- yt-resi
 Sigma <- c(sigma1,sigma2)
 ### pooled estimate
 sigmasq = (n1*sigma1^2+n2*sigma2^2)/(n1+n2)
 AIC <- AIC1+AIC2
 cat(" ","\n")
 cat("Overal MLE of sigma: ",sqrt(sigmasq),"\n")
 cat("Overall AIC: ",AIC,"\n")
 }else{cat("No threshold found: Try again with a larger k0.","\n")}
uTAR <- list(data=y,arorder=c(p1,p2),delay=d,residuals = resi, coefs = Phi, sigma=Sigma,
nobs = Size, model1 = m1a, model2 = m1b,thr=thr, D=D, RSS=RSS,AIC = AIC,
cnst = rep(include.mean,2), sresi=sresi)
}


#' Search for Threshold Value of A Two-Regime SETAR Model
#'
#' Search for the threshold of a SETAR model for a given range of candidates for threshold values,
#' and perform recursive LS estimation.
#' The program uses a grid to search for threshold value.
#' It is a conservative approach, but might be more reliable than the Li and Tong (2016) procedure.
#' @param y a vector of time series.
#' @param p1,p2 AR-orders of regime 1 and regime 2.
#' @param d delay for threshold variable, default is 1.
#' @param thrV threshold variable. if it is not null, thrV must have the same length as that of y.
#' @param thrQ lower and upper limits for the possible threshold values.
#' @param Trim lower and upper trimming to control the sample size in each regime.
#' @param include.mean a logical value for including constant term.
#' @references
#' Li, D., and Tong. H. (2016) Nested sub-sample search algorithm for estimation of threshold models. \emph{Statisitca Sinica}, 1543-1554.
#' @return uTAR.grid returns a list with components:
#' \item{data}{the data matrix, y.}
#' \item{arorder}{AR orders of regimes 1 and 2.}
#' \item{residuals}{estimated innovations.}
#' \item{coefs}{a 2-by-(p+1) matrices. The first row show the estimation results in regime 1, and the second row shows these in regime 2.}
#' \item{sigma}{estimated innovational covariance matrices of regimes 1 and 2.}
#' \item{nobs}{numbers of observations in regimes 1 and 2.}
#' \item{delay}{the delay for threshold variable.}
#' \item{model1,model2}{estimated models of regimes 1 and 2.}
#' \item{cnst}{logical values indicating whether the constant terms are included in regimes 1 and 2.}
#' \item{thr}{threshold value.}
#' \item{D}{a set of possible threshold values.}
#' \item{RSS}{residual sum of squares.}
#' \item{information}{information criterion.}
#' \item{sresi}{standardized residuals.} 
#' @examples
#' arorder=rep(1,2)
#' ar.coef=matrix(c(0.7,-0.8),2,1)
#' y=uTAR.sim(100,arorder,ar.coef,1,0)
#' uTAR.grid(y$series,1,1,1,y$series,c(0,1),c(0.2,0.8),TRUE)
#' @export
"uTAR.grid" <- function(y,p1,p2,d=1,thrV=NULL,thrQ=c(0,1),Trim=c(0.1,0.9),include.mean=T){
if(is.matrix(y))y=y[,1]
p = max(p1,p2)
if(p < 1)p=1
ist=max(p,d)+1
nT <- length(y)
yt <- y[ist:nT]
nobe <- nT-ist+1
Phi <- matrix(0,2,p+1)
x1 <- NULL
for (i in 1:p1){
 x1=cbind(x1,y[(ist-i):(nT-i)])
 }
x1 <- as.matrix(x1)
if(include.mean){x1 <- cbind(rep(1,nobe),x1)}
k1 <- ncol(x1)
x2 <- NULL
for (i in 1:p2){
 x2 <- cbind(x2,y[(ist-i):(nT-i)])
 }
x2 <- as.matrix(x2)
if(include.mean){x2 <- cbind(rep(1,nobe),x2)}
k2 <- ncol(x2)
#
if(length(thrV) < nT){
    tV <- y[(ist-d):(nT-d)]
    }else{
    tV <- thrV[ist:nT]
    }
### narrow down possible threshold range
q1=c(min(tV),max(tV))
if(thrQ[1] > 0)q1[1] <- quantile(tV,thrQ[1])
if(thrQ[2] < 1)q1[2] <- quantile(tV,thrQ[2])
idx <- c(1:length(tV))[tV <= q1[1]]
jdx <- c(1:length(tV))[tV >= q1[2]]
trm <- c(idx,jdx)
yt <- yt[-trm]
if(k1 == 1){x1 <- x1[-trm]
         }else{x1 <- x1[-trm,]}
if(k2 == 1){x2 <- x2[-trm]
         }else{x2 <- x2[-trm,]}
tV <- tV[-trm]
### sorting
DD <- sort(tV,index.return=TRUE)
D <- DD$x
Dm <- DD$ix
yt <- yt[Dm]
if(k1 == 1){x1 <- x1[Dm]
          }else{x1 <- x1[Dm,]}
if(k2 == 1){x2 <- x2[Dm]
          }else{x2 <- x2[Dm,]}
nobe <- length(yt)
q2 = quantile(tV,Trim)
jst <- max(c(1:nobe)[D < q2[1]])
jend <- min(c(1:nobe)[D > q2[2]])
cat("jst, jend: ",c(jst,jend),"\n")
nthr <- jend-jst+1
RSS=NULL
### initial estimate with thrshold D[jst]
if(k1 == 1){
          X <- as.matrix(x1[1:jst])
          }else{ X <- x1[1:jst,]}
Y <- as.matrix(yt[1:jst])
XpX = t(X)%*%X
XpY = t(X)%*%Y
Pk1 <- solve(XpX)
beta1 <- Pk1%*%XpY
resi1 <- Y - X%*%beta1
if(k2 == 1){
        X <- as.matrix(x2[(jst+1):nobe])
        }else{ X <- x2[(jst+1):nobe,]}
Y <- as.matrix(yt[(jst+1):nobe])
XpX <- t(X)%*%X
XpY <- t(X)%*%Y
Pk2 <- solve(XpX)
beta2 <- Pk2%*%XpY
resi2 <- Y - X%*%beta2
RSS <- sum(c(resi1^2,resi2^2))
#### recursive
#### Moving one observations from regime 2 to regime 1
for (i in (jst+1):jend){
 if(k1 == 1){xk1 <- matrix(x1[i],1,1)
          }else{
            xk1 <- matrix(x1[i,],k1,1)
          }
if(k2 == 1){xk2 <- matrix(x2[i],1,1)
          }else{
         xk2 <- matrix(x2[i,],k2,1)
         }
 tmp <- Pk1%*%xk1
 deno <- 1 + t(tmp)%*%xk1
 er <- t(xk1)%*%beta1 - yt[i]
 beta1 <- beta1 - tmp*c(er)/c(deno[1])
 Pk1 <- Pk1 - tmp%*%t(tmp)/c(deno[1])
 if(k1 == 1){X <- as.matrix(x1[1:i])
            }else{X <- x1[1:i,]}
 resi1 <- yt[1:i] - X%*%beta1
 tmp <- Pk2 %*% xk2
 deno <- t(tmp)%*%xk2 - 1
 er2 <- t(xk2)%*%beta2 - yt[i]
 beta2 <- beta2 - tmp*c(er2)/c(deno[1])
 Pk2 <- Pk2 - tmp%*%t(tmp)/c(deno[1])
 if(k2 == 1){X <- as.matrix(x2[(i+1):nobe])
            }else{X <- x2[(i+1):nobe,]}
 resi2 <- yt[(i+1):nobe] - X%*%beta2
 RSS <- c(RSS, sum(c(resi1^2,resi2^2)))
  }
##cat("RSS: ",RSS,"\n")
j=which.min(RSS)+(jst-1)
thr <- D[j]
cat("Estimated Threshold: ",thr,"\n")
### Final estimation results
  resi = rep(0,nobe)
  sresi <- rep(0,nobe)
  idx <- c(1:j)
  n1 <- j
  n2 <- nobe-j
  X <- x1
  m1a <- lm(yt~-1+.,data=data.frame(x1),subset=idx)
  resi[idx] <- m1a$residuals
 cat("Regime 1: ","\n")
 Phi[1,1:k1] <- m1a$coefficients
 coef1 <- summary(m1a)
 sresi[idx] <- m1a$residuals/coef1$sigma
 print(coef1$coefficients)
 cat("nob1 & sigma1: ",c(n1, coef1$sigma), "\n")
  m1b <- lm(yt~-1+.,data=data.frame(x2),subset=-idx)
  resi[-idx] <- m1b$residuals
  Phi[2,1:k2] <- m1b$coefficients
  cat("Regime 2: ","\n")
  coef2 <- summary(m1b)
  sresi[-idx] <- m1b$residuals/coef2$sigma
  print(coef2$coefficients)
  cat("nob2 & sigma2: ",c(n2,coef2$sigma),"\n")
 fitted.values <- yt-resi
 ### pooled estimate
 sigma2 = (n1*coef1$sigma+n2*coef2$sigma)/(n1+n2)
 d1 = log(sigma2)
 aic = d1 + 2*(p1+p2)/nT
 bic = d1 + log(nT)*(p1+p2)/nT
 hq = d1 + 2*log(log(nT))*(p1+p2)/nT
 infc = c(aic,bic,hq)
 cat(" ","\n")
 cat("The fitted TAR model ONLY uses data in the specified threshold range!!!","\n")
 cat("Overal MLE of sigma: ",sqrt(sigma2),"\n")
 cat("Overall information criteria(aic,bic,hq): ",infc,"\n")
uTAR.grid <- list(data=y,arorder=c(p1,p2), residuals = resi, coefs = Phi,
sigma=c(coef1$sigma,coef2$sigma), nobs=c(n1,n2),delay=d,model1 = m1a,cnst = rep(include.mean,2), 
model2 = m1b,thr=thr, D=D, RSS=RSS,information=infc,sresi=sresi)
}

#' General Estimation of TAR Models
#'
#' General estimation of TAR models with known threshold values. 
#' It perform LS estimation of a univariate TAR model, and can handle multiple regimes.
#' @param y time series.
#' @param arorder AR order of each regime. The number of regime is the length of arorder.
#' @param thr given threshould(s). There are k-1 threshold for a k-regime model.
#' @param d delay for threshold variable, default is 1.
#' @param thrV external threhold variable if any. If it is not NULL, thrV must have the same length as that of y.
#' @param include.mean a logical value indicating whether constant terms are included. Default is TRUE.
#' @param output a logical value for output. Default is TRUE.
#' @return uTAR.est returns a list with components:
#' \item{data}{the data matrix, y.}
#' \item{k}{the dimension of y.}
#' \item{arorder}{AR orders of regimes 1 and 2.}
#' \item{coefs}{a (p*k+1)-by-(2k) matrices. The first row show the estimation results in regime 1, and the second row shows these in regime 2.}
#' \item{sigma}{estimated innovational covariance matrices of regimes 1 and 2.}
#' \item{thr}{threshold value.}
#' \item{residuals}{estimated innovations.}
#' \item{sresi}{standardized residuals.}
#' \item{nobs}{numbers of observations in different regimes.}
#' \item{delay}{delay for threshold variable.}
#' \item{cnst}{logical values indicating whether the constant terms are included in different regimes.}
#' \item{AIC}{AIC value.}
#' @examples
#' arorder=rep(1,2)
#' ar.coef=matrix(c(0.7,-0.8),2,1)
#' y=uTAR.sim(100,arorder,ar.coef,1,0)
#' est=uTAR.est(y$series,arorder,0,1)
#' @export
"uTAR.est" <- function(y,arorder=c(1,1),thr=c(0),d=1,thrV=NULL,include.mean=c(TRUE,TRUE),output=TRUE){
if(is.matrix(y))y <- y[,1]
p <- max(arorder)
if(p < 1)p <-1
k <- length(arorder)
ist <- max(p,d)+1
nT <- length(y)
Y <- y[ist:nT]
X <- NULL
for (i in 1:p){
    X <- cbind(X,y[(ist-i):(nT-i)])
    }
X <- as.matrix(X)
nobe <- nT-ist+1
if(length(thrV) < nT){
 thrV <- y[(ist-d):(nT-d)]
 if(output) cat("Estimation of a SETAR model: ","\n")
 }else{thrV <- thrV[ist:nT]}
k1 <- length(thr)
### house keeping
 coefs <- NULL
 sigma <- NULL
 nobs <- NULL
 resi <- rep(0,nobe)
 sresi <- rep(0,nobe)
 npar <- NULL
 AIC <- 99999
#
if(k1 > (k-1)){
 cat("Only the first: ",k1," thresholds are used.","\n")
 THr <- min(thrV[1:nobe])-1
 THr <- c(THr,thr[1:(k-1)],(max(thrV[1:nobe])+1))
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
  x1 <- as.matrix(X)
  if(p > 1)x1 <- x1[,c(1:arorder[j])]
  cnst <- rep(1,nobe)
  if(include.mean[j]){x1 <- cbind(cnst,x1)}
  np <- ncol(x1)
  npar <- c(npar,np)
  beta <- rep(0,p+1)
  sig <- NULL
  if(n1 <= ncol(x1)){
   cat("Regime ",j," does not have sufficient observations!","\n")
   }
   else{
   mm <- lm(Y~-1+.,data=data.frame(x1),subset=kdx)
   c1 <- summary(mm)
   if(output){
    cat("Estimation of regime: ",j," with sample size: ",n1,"\n")
    print(c1$coefficients)
    cat("Sigma estimate: ",c1$sigma,"\n")
   }
   resi[kdx] <- c1$residuals
   sresi[kdx] <- c1$residuals/c1$sigma
   sig <- c1$sigma
   beta[1:np] <- c1$coefficients[,1]
   coefs <- rbind(coefs,beta)
   sigma <- c(sigma,sig)
   }
  }
  AIC <- 0
  for (j in 1:k){
   sig <- sigma[j]
   n1 <- nobs[j]
   np <- npar[j]
   if(!is.null(sig)){s2 <- sig^2*(n1-np)/n1
               AIC <- AIC + log(s2)*n1 + 2*np
	       }
	}
  if(output) cat("Overal AIC: ",AIC,"\n")
}
uTAR.est <- list(data=y,k = k, arorder=arorder,coefs=coefs,sigma=sigma,thr=thr,residuals=resi,
sresi=sresi,nobs=nobs,delay=d,cnst=include.mean, AIC=AIC)
}



#' Prediction of a Fitted Univariate TAR Model
#'
#' Prediction of a fitted univariate TAR model.
#' @param model univariate TAR model.
#' @param orig forecast origin.
#' @param h forecast horizon.
#' @param iternations number of iterations.
#' @param ci confidence level.
#' @param output a logical value for output, default is TRUE.
#' @return uTAR.pred returns a list with components:
#' \item{model}{univariate TAR model.}
#' \item{pred}{prediction.}
#' \item{Ysim}{fitted y.}
#' @examples
#' arorder=rep(1,2)
#' ar.coef=matrix(c(0.7,-0.8),2,1)
#' y=uTAR.sim(100,arorder,ar.coef,1,0)
#' est=uTAR.est(y$series,arorder,0,1)
#' pred=uTAR.pred(est,100,1,100,0.95,TRUE)
#' @export
"uTAR.pred" <- function(model,orig,h=1,iterations=3000,ci=0.95,output=TRUE){
## ci: probability of pointwise confidence interval coefficient
### only works for self-exciting uTAR models.
###
y <- model$data
arorder <- model$arorder
coefs <- model$coefs
sigma <- model$sigma
thr <- c(model$thr)
include.mean <- model$cnst
d <- model$delay
p <- max(arorder)
k <- length(arorder)
nT <- length(y)
if(orig < 1)orig <- nT
if(orig > nT) orig <- nT
if(h < 1) h <- 1
### compute the predictions. Use simulation if the threshold is predicted
Ysim <- matrix(0,h,iterations)
et <- rnorm(h*iterations)
for (it in 1:iterations){
  etcnt <- (it-1)*h
  yp <- y[1:orig]
  for (ii in 1:h){
   t <- orig+ii
   thd <- yp[(t-d)]
   JJ <- 1
   for (j in 1:(k-1)){
     if(thd > thr[j]){JJ <- j+1}
     }
   ino <- etcnt+ii
   at <- et[ino]*sigma[JJ]
   x <- NULL
   pJ <- arorder[JJ]
   npar <- pJ
   if(include.mean[JJ]){ x <- 1
                         npar <- npar+1
		        }
   for (i in 1:pJ){
    x <- c(x,yp[(t-i)])
    }
   yhat <- sum(x*coefs[JJ,1:npar])+at
   yp <- rbind(yp,yhat)
   Ysim[ii,it] <- yhat
   }
 }
### summary
 pred <- NULL
 CI <- NULL
 pr <- (1-ci)/2
 prob <- c(pr, 1-pr)
 for (ii in 1:h){
   pred <- rbind(pred,mean(Ysim[ii,]))
   int <- quantile(Ysim[ii,],prob=prob)
   CI <- rbind(CI,int)
 }
if(output){
cat("Forecast origin: ",orig,"\n")
cat("Predictions: 1-step to ",h,"-step","\n")
Pred <- cbind(1:h,pred)
colnames(Pred) <- c("step","forecast")
print(Pred)
bd <- cbind(c(1:h),CI)
colnames(bd) <- c("step", "Lowb","Uppb")
cat("Pointwise ",ci*100," % confident intervals","\n")
print(bd)
}
uTAR.pred <- list(data=y,pred = pred,Ysim=Ysim)
}

#' @export
"NeSS" <- function(y,x1,x2=NULL,thrV=NULL,Trim=c(0.15,0.85),k0=300,include.mean=TRUE,thrQ=c(0,1),score="AIC"){
### Based on Li and Tong (2016)method to narrow down the candidates for 
### threshold value.
### y: dependent variable (can be multivariate)
### x1: the regressors that allow for coefficients to change
### x2: the regressors for explanatory variables
### thrV: threshold variable,default is simply the time index
### Trim: quantiles for lower and upper bound of threshold
### k0: the maximum number of candidates for threhold at the output
### include.mean: switch for including the constant term.
### thrQ: lower and upper quantiles to search for threshold
####
#### return: a set of possible threshold values
####
if(!is.matrix(y))y <- as.matrix(y)
ky <- ncol(y)
if(!is.matrix(x1))x1 <- as.matrix(x1)
if(length(x2) > 0){
 if(!is.matrix(x2))x2 <- as.matrix(x2)
  }
n <- nrow(y)
n1 <- nrow(x1)
n2 <- n1
if(!is.null(x2))n2 <- nrow(x2)
n <- min(n,n1,n2)
##
k1 <- ncol(x1)
k2 <- 0
if(!is.null(x2)) k2 <- ncol(x2)
k <- k1+k2
##
if(length(thrV) < 1){thrV <- c(1:n)
                   }else{thrV <- thrV[1:n]}
### set threshold range
idx <- NULL
jdx <- NULL
if(thrQ[1] > 0){low <- quantile(thrV,thrQ[1])
                idx <- c(1:length(thrV))[thrV < low]
	}
if(thrQ[2] < 1){upp <- quantile(thrV,thrQ[2])
                jdx <- c(1:length(thrV))[thrV > upp]
	}
jrm <- c(idx,jdx)
if(!is.null(jrm)){
		if(ky == 1){y <- as.matrix(y[-jrm])
		          }else{y <- y[-jrm,]}
                if(k1 == 1){x1 <- as.matrix(x1[-jrm])
		          }else{x1 <- x1[-jrm,]}
	        if(k2 > 0){
		     if(k2 == 1){x2 <- as.matrix(x2[-jrm])
		                }else{x2 <- x2[-jrm,]}
		     }
		thrV <- thrV[-jrm]
                }
thr <- sort(thrV)
n <- length(thrV)
### use k+2 because of including the constant term.
n1 <- floor(n*Trim[1])
if(n1 < (k+2)) n1 = k+2
n2 <- floor(n*Trim[2])
if((n-n2) > (k+2)){n2 = n -(k+2)}
D <- thr[n1:n2]
####
X=cbind(x2,x1)
if(include.mean){k = k+1
    X=cbind(rep(1,n),X)}
qprob=c(0.25,0.5,0.75)
#### Scalar case
if(ky == 1){
X=data.frame(X)
Y <- y[,1]
m1 <- lm(Y~-1+.,data=X)
R1 <- sum(m1$residuals^2)
Jnr <- rep(R1,length(qprob))
#cat("Jnr: ",Jnr,"\n")
while(length(D) > k0){
  qr <- quantile(D,qprob)
##  cat("qr: ",qr,"\n")
  for (i in 1:length(qprob)){
    idx <- c(1:n)[thrV <= qr[i]]
     m1a <- lm(Y~-1+.,data=X,subset=idx)
     m1b <- lm(Y~-1+.,data=X,subset=-idx)
     Jnr[i] = Jnr[i] - sum(c(m1a$residuals^2,m1b$residuals^2))
     }
##     cat("in Jnr: ",Jnr,"\n")
  if(Jnr[1] >= max(Jnr[-1])){
    D <- D[D <= qr[2]]
    }
    else{ if(Jnr[2] >= max(Jnr[-2])){
           D <- D[((qr[1] <= D) && (D <= qr[3]))]
            }
            else{ D <- D[D >= qr[2]]}
     }
## cat("n(D): ",length(D),"\n")
 }
}
#### Multivariate case
if(ky > 1){
ic = 2
if(score=="AIC")ic=1
m1 <- MlmNeSS(y,X,SD=FALSE,include.mean=include.mean)
R1 <- m1$score[ic]
Jnr <- rep(R1,length(qprob))
while(length(D) > k0){
  qr <- quantile(D,qprob)
  for (i in 1:length(qprob)){
    idx <- c(1:n)[thrV <= qr[i]]
     m1a <- MlmNeSS(y,X,subset=idx,SD=FALSE,include.mean=include.mean)
     m1b <- MlmNeSS(y,X,subset=-idx,SD=FALSE,include.mean=include.mean)
     Jnr[i] <- Jnr[i] - sum(m1a$score[ic],m1b$score[ic])
     }
  if(Jnr[1] >= max(Jnr[-1])){
    D <- D[D <= qr[2]]
    }
    else{ if(Jnr[2] >= max(Jnr[-2])){
           D <- D[D <= qr[3]]
	   D <- D[D >= qr[1]]
            }
            else{ D <- D[D >= qr[2]]}
     }
 }
### end of Multivariate case
}
D 
}


