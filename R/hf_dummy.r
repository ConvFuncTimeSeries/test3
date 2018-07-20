#' Create Dummy Variables for High-Frequency iIntraday Seasonality
#'
#' @param int length of time interval in minutes.
#' @param Fopen number of dummies from the market open.
#' @param Tend number of dummies to the market close.
#' @param days number of trading days in the data.
#' @param pooled
#' @param skipmin
#' @export
"hfDummy" <- function(int=1,Fopen=10,Tend=10,days=1,pooled=1,skipmin=0){
nintval=(6.5*60-skipmin)/int
Days=matrix(1,days,1)
X=NULL
ini=rep(0,nintval)
for (i in 1:Fopen){
x1=ini
x1[i]=1
Xi=kronecker(Days,matrix(x1,nintval,1))
X=cbind(X,Xi)
}
for (j in 1:Tend){
x2=ini
x2[nintval-j+1]=1
Xi=kronecker(Days,matrix(x2,nintval,1))
X=cbind(X,Xi)
}
X1=NULL
if(pooled > 1){
nh=floor(Fopen/pooled)
rem=Fopen-nh*pooled
if(nh > 0){
for (ii in 1:nh){
ist=(ii-1)*pooled
y=apply(X[,(ist+1):(ist+pooled)],1,sum)
X1=cbind(X1,y)
}
}
if(rem > 0){
X1=cbind(X1,X[,(Fopen-rem):Fopen])
}
nh=floor(Tend/pooled)
rem=Tend-nh*pooled
if(nh > 0){
for (ii in 1:nh){
ist=(ii-1)*pooled
y=apply(X[,(Fopen+ist+1):(Fopen+ist+pooled)],1,sum)
X1=cbind(X1,y)
}
}
if(rem > 0){
X1=cbind(X1,X[,(Fopen+Tend-rem):(Fopen+Tend)])
}

X=X1
}

hfDummy <- X
}


######### factorial product #####
#' @export
"factorial" <- function(n,log=T){
x=c(1:n)
if(log){
x=log(x)
y=cumsum(x)
}
else{
y=cumprod(x)
}
y[n]
}



