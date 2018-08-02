#' Generate Univeraite 2-regime Markov Switching Models
#'
#' Generate univeraite 2-regime Markov switching models.
#' @param nob number of observations.
#' @param order AR order for each regime.
#' @param phi1,phi2 AR coefficients.
#' @param epsilon transition probabilities (switching out of regime 1 and 2).
#' @param sigma standard errors for each regime.
#' @param cnst constant term for each regime.
#' @param ini burn-in period.
#' @return MSM.sim returns a list with components:
#' \item{series}{a time series following SETAR model.}
#' \item{at}{innovation of the time seres.}
#' \item{state}{states for the time series.}
#' \item{epsilon}{transition probabilities (switching out of regime 1 and 2).}
#' \item{sigma}{standard error for each regime.}
#' \item{cnst}{constant terms.}
#' \item{order}{AR-order for each regime.}
#' \item{phi1, phi2}{the AR coefficients for two regimes.}
#' @examples
#' y=MSM.sim(100,c(1,1),0.7,-0.5,c(0.5,0.6),c(1,1),c(0,0),500)
#' @export
"MSM.sim" <- function(nob,order=c(1,1),phi1=NULL,phi2=NULL,epsilon=c(0.1,0.1),sigma=c(1,1),cnst=c(0,0),ini=500){
nT <- nob+ini
et <- rnorm(nT)
p <- max(order)
if(p < 1) p <- 1
if(order[1] < 0)order[1] <- 1
if(order[2] < 0)order[2] <- 1
if(sigma[1] < 0)sigma[1] <- 1
if(sigma[2] < 0)sigma[2] <- 1
ist <- p+1
y <- et[1:p]
at <- et
state <- rbinom(p,1,0.5)+1
for (i in ist:nT){
  p1 <- epsilon[1]
  sti <- state[i-1]
  if(sti == 2)p1 <- epsilon[2]
  sw <- rbinom(1,1,p1)
  st <- sti
  if(sti == 1){
   if(sw == 1)st <- 2
   }
   else{
    if(sw == 1)st <- 1
    }
 if(st == 1){
  tmp <- cnst[1]
  at[i] <- et[i]*sigma[1]
  if(order[1] > 0){
   for (j in 1:order[1]){
    tmp = tmp + phi1[j]*y[i-j]
    }
   }
   y[i] <- tmp+at[i]
   }
   else{ tmp <- cnst[2]
   at[i] <- et[i]*sigma[2]
   if(order[2] > 0){
    for (j in 1:order[2]){
     tmp <- tmp + phi2[j]*y[i-j]
     }
    }
    y[i] <- tmp + at[i]
   }
   state <- c(state,st)
 }
MSM.sim <- list(series=y[(ini+1):nT],at = at[(ini+1):nT], state=state[(ini+1):nT], epsilon=epsilon,
sigma=sigma,cnst=cnst,order=order,phi1=phi1,phi2=phi2)
}


#' Fitting Univariate Autoregressive Markov Switching Models
#'
#' Fit autoregressive Markov switching models to a univariate time series using the package MSwM.
#' @param y a time series.
#' @param p AR order.
#' @param nregime the number of regimes.
#' @param include.mean a logical value for including constant terms.
#' @param sw logical values for whether coefficients are switching. The length of \code{sw} has to be equal of the number of coefficients in the model plus include.mean.
#' @return \code{MSM.fit} returns an object of class code{MSM.lm} or \code{MSM.glm}, depending on the input model.
#' @examples
#' y=MSM.sim(100,c(1,1),0.7,-0.5,c(0.5,0.6),c(1,1),c(0,0),500)
#' library(parallel)
#' library(MSwM)
#' MSM.fit(y$series,1,2,TRUE,c(TRUE,TRUE,TRUE))
#' @export
"MSM.fit" <- function(y,p,nregime=2,include.mean=T,sw=NULL){
#require(MSwM)
if(is.matrix(y))y <- y[,1]
if(p < 0)p <- 1
ist <- p+1
nT <- length(y)
X <- y[ist:nT]
if(include.mean) X <- cbind(X,rep(1,(nT-p)))
for (i in 1:p){
 X <- cbind(X,y[(ist-i):(nT-i)])
 }
if(include.mean){
  colnames(X) <- c("y","cnst",paste("lag",1:p))
  }else{
  colnames(X) <- c("y",paste("lag",1:p))
  }
X <- data.frame(X)
mo <- lm(y~-1+.,data=X)
### recall that dependent variable "y" is a column in X.
npar <- ncol(X)
if(is.null(sw))sw = rep(TRUE,npar)
mm <- msmFit(mo,k=nregime,sw=sw)
mm
}
