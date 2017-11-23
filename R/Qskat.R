#' Fit a null linear regression model
#'
#' Fit a null linear model to be used for variant set association test
#' @param  Y continuous outcome
#' @param  X covariates to be adjusted, setting X=NULL with no covariate
#' @keywords KAT.cnull
#' @export
KAT.cnull <- function(Y,X){
  if(is.null(X)){
    U0 = Y - mean(Y)
    DF0 = length(Y)-1
    s2 = sum(U0^2)/DF0
    Ux = matrix(1/sqrt(length(Y)), length(Y),1)
  } else{
    X0 = cbind(1,X)
    DF0 = length(Y)-dim(X0)[2]
    a0 = svd(X0,nv=0)
    U0 = lm(Y~X)$res
    s2 = sum(U0^2)/DF0
    Ux = a0$u
  }
  return(list(U0=U0,s2=s2,Ux=Ux, DF0=DF0, mode='C'))
}


#' Sequence kernel association test (SKAT) for quantitative trait based on marginal Wald-statistics
#'
#' Compute the significance p-value for SKAT based on marginal Wald-statistics
#' @param  obj a fitted null binomial model using KAT.cnull()
#' @param  G genotype matrix, sample in rows, variant in columns
#' @param  W.beta Beta parameters for variant weights
#' @return SKAT p-value
#' @keywords SKAT
#' @export
#' @references
#' Wu, M. C., Lee, S., Cai, T., Li, Y., Boehnke, M., and Lin, X. (2011) Rare Variant Association Testing for Sequencing Data Using the Sequence Kernel Association Test (SKAT). American Journal of Human Genetics, 89, 82-93.
#'
#' Wu, M. C., Kraft, P., Epstein, M. P.,Taylor, D., M., Chanock, S. J., Hunter, D., J., and Lin, X. (2010) Powerful SNP Set Analysis for Case-Control Genome-wide Association Studies. American Journal of Human Genetics, 86, 929-942.
#'
#' Duchesne, P. and Lafaye De Micheaux, P. (2010) Computing the distribution of quadratic forms: Further comparisons between the Liu-Tang-Zhang approximation and exact methods, Computational Statistics and Data Analysis, 54, 858-862.
#'
#' Wu,B., Pankow,J.S., Guan,W. (2015) Sequence kernel association analysis of rare variant set based on the marginal regression model for binary traits. Genetic Epidemiology, 39(6), 399-405.
#'
#' Wu,B., Guan,W., Pankow,J.S. (2016) On efficient and accurate calculation of significance p-values for sequence kernel association test of variant set. Annals of human genetics, 80(2), 123-135.
SKATT <- function(obj,G, W.beta=c(1.25,25.5)){
  N = dim(G)[2]; maf = colMeans(G)/2
  W = maf^(W.beta[1]-1)*(1-maf)^(W.beta[2]-1);  W = W/sum(W)*N
  SSE = obj$DF0*obj$s2
  tmp = t(obj$Ux)%*%G
  Ge = G - obj$Ux%*%tmp
  Ge1 = sqrt(colSums(Ge^2))
  bta = colSums(obj$U0*G)/Ge1^2
  Ts = bta/sqrt((SSE-bta^2*Ge1^2)/(obj$DF0-1))*Ge1
  ## Zs = -sign(Ts)*qnorm(pt(-abs(Ts),obj$DF0-1))
  Q = sum(Ts^2*W^2)
  Gs = t(Ge)%*%Ge
  Gs = t(Gs/Ge1)/Ge1
  R = t(Gs*W)*W
  lam = svd(R, nu=0,nv=0)$d
  KAT.pval(Q, lam)
}

## Sequence kernel association test (SKAT) for quantitative trait based on several marginal statistics
##
## Compute the significance p-value for SKAT based on several marginal statistics
## @param  obj a fitted null binomial model using KAT.cnull()
## @param  G genotype matrix, sample in rows, variant in columns
## @param  W.beta Beta parameters for variant weights
## @return SKAT p-values based on three marginal statistics
## @export
## several with moment adjustment
SKATM <- function(obj,G, W.beta=c(1.25,25)){
  N = dim(G)[2]; maf = colMeans(G)/2
  W = maf^(W.beta[1]-1)*(1-maf)^(W.beta[2]-1);  W = W/sum(W)*N
  SSE = obj$DF0*obj$s2
  tmp = t(obj$Ux)%*%G
  Ge = G - obj$Ux%*%tmp
  Ge1 = sqrt(colSums(Ge^2))
  bta = colSums(obj$U0*G)/Ge1^2
  Ts = bta/sqrt((SSE-bta^2*Ge1^2)/(obj$DF0-1))*Ge1
  Gs = t(Ge)%*%Ge
  Gs = t(Gs/Ge1)/Ge1
  R = t(Gs*W)*W
  lam = svd(R, nu=0,nv=0)$d
  Q0 = sum(W^2*Ts^2)
  Q1 = sum(W^2*Ts^2*(obj$DF0-3)/(obj$DF0-1))
  Zs = -sign(Ts)*qnorm(pt(-abs(Ts),obj$DF0-1))
  Q2 = sum(W^2*Zs^2)
  KAT.pval(c(Q0,Q1,Q2), lam)
}

#' Optimal sequence kernel association test (SKAT-O) for quantitative trait based on marginal t-statistics
#'
#' Compute the significance p-value for SKAT-O based on marginal t-statistics. The computational algorithm
#' is described in detail at Wu et. al (2015).
#' @param  obj a fitted null linear model using KAT.cnull()
#' @param  G genotype matrix, sample in rows, variant in columns
#' @param  W.beta Beta parameters for variant weights
#' @param  rho weights for burden test
#' @return SKATOT p-value
#' @keywords SKATO-T
#' @export
#' @references
#' Lee, S., Wu, M. C., and Lin, X. (2012) Optimal tests for rare variant effects in sequencing association studies. Biostatistics, 13, 762-775.
#' 
#' Wu,B., Pankow,J.S., Guan,W. (2015) Sequence kernel association analysis of rare variant set based on the marginal regression model for binary traits. Genetic Epidemiology, 39(6), 399-405.
#'
#' Wu,B., Guan,W., Pankow,J.S. (2016) On efficient and accurate calculation of significance p-values for sequence kernel association test of variant set. Annals of human genetics, in press.
#' @examples
#' library(CompQuadForm)
#' Y = rnorm(5000); X = matrix(rnorm(10000),5000,2)
#' G = matrix(rbinom(100000,2,0.01), 5000,10)
#' Y = Y + G[,1]*0.4 + G[,10]*0.3 + G[,2]*0.2
#' SKATT(KAT.cnull(Y,X), G, c(1.5,25.5))
#' SKATOT(KAT.cnull(Y,X), G, c(1.5,25.5))
#' ## library(SKAT)
#' ## SKAT(G, SKAT_Null_Model(Y~X, out_type='C'), method='davies')$p.value
#' ## SKAT(G, SKAT_Null_Model(Y~X, out_type='C'), method='optimal.adj')$p.value
SKATOT = function(obj,G, W.beta=c(1.5,25.5), rho=c(0,0.1^2,0.2^2,0.3^2,0.4^2,0.5^2,0.5,1)){
  N = dim(G)[2]; maf = colMeans(G)/2
  W = maf^(W.beta[1]-1)*(1-maf)^(W.beta[2]-1);  W = W/sum(W)*N
  SSE = obj$DF0*obj$s2
  tmp = t(obj$Ux)%*%G
  Ge = G - obj$Ux%*%tmp
  Ge1 = sqrt(colSums(Ge^2))
  bta = colSums(obj$U0*G)/Ge1^2
  Ts = bta/sqrt((SSE-bta^2*Ge1^2)/(obj$DF0-1))*Ge1
  ## Zs = -sign(Ts)*qnorm(pt(-abs(Ts),obj$DF0-1))
  ## Z = Zs*W
  Z = Ts*W
  Gs = t(Ge)%*%Ge
  Gs = t(Gs/Ge1)/Ge1
  R = t(Gs*W)*W
  ##
  K = length(rho); K1 = K
  Qs = sum(Z^2); Qb = sum(Z)^2; Qw = (1-rho)*Qs + rho*Qb
  pval = rep(0,K)
  Rs = rowSums(R); R1 = sum(Rs); R2 = sum(Rs^2); R3 = sum(Rs*colSums(R*Rs))
  RJ2 = outer(Rs,Rs,'+')/N
  ## min-pval
  if(rho[K]>=1){
    K1 = K-1
    pval[K] = pchisq(Qb/R1, 1, lower.tail=FALSE)
  }
  Lamk = vector('list', K1);  rho1 = rho[1:K1]
  tmp = sqrt(1-rho1+N*rho1) - sqrt(1-rho1)
  c1 = sqrt(1-rho1)*tmp;  c2 = tmp^2*R1/N^2
  for(k in 1:K1){
    mk = (1-rho[k])*R + c1[k]*RJ2 + c2[k]
    Lamk[[k]] = pmax(svd(mk, nu=0,nv=0)$d, 0)
    pval[k] = KAT.pval(Qw[k],Lamk[[k]])
  }
  Pmin = min(pval)
  qval = rep(0,K1)
  for(k in 1:K1) qval[k] = Liu.qval.mod(Pmin, Lamk[[k]])
  lam = pmax(svd(R-outer(Rs,Rs)/R1, nu=0,nv=0)$d[-N], 0)
  tauk = (1-rho1)*R2/R1 + rho1*R1;  vp2 = 4*(R3/R1-R2^2/R1^2)
  MuQ = sum(lam);  VarQ = sum(lam^2)*2
  sd1 = sqrt(VarQ)/sqrt(VarQ+vp2)
  if(K1<K){
    q1 = qchisq(Pmin,1,lower=FALSE)
    T0 = Pmin
  } else{
    tmp = ( qval-(1-rho)*MuQ*(1-sd1)/sd1 )/tauk
    q1 = min(tmp)
    T0 = pchisq(q1,1,lower=FALSE)
  }
  katint = function(xpar){
    eta1 = sapply(xpar, function(eta0) min((qval-tauk*eta0)/(1-rho1)))
    x = (eta1-MuQ)*sd1 + MuQ
    KAT.pval(x,lam)*dchisq(xpar,1)
  }
  p.value = try({ T0 + integrate(katint, 0,q1,  subdivisions=1e3,abs.tol=1e-25)$val }, silent=TRUE)
  prec = 1e-4
  while(class(p.value)=='try-error'){
    p.value = try({ T0 + integrate(katint, 0,q1, abs.tol=Pmin*prec)$val }, silent=TRUE)
    prec = prec*2
  }
  return( min(p.value, Pmin*K) )
}

