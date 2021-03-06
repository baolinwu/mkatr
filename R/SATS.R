#' SNP-set association tests using GWAS Summary data
#'
#' Compute p-values for the SNP-set tests using GWAS Z-statistics: squared sum test (S2T), sum test (ST), and adaptive test (AT).
#' With less than two values of \eqn{\rho}'s, it outputs only p-values for S2T and ST.
#' @param  Z GWAS summary Z-statistics for a set of SNPs
#' @param  R SNP pairwise LD matrix
#' @param  W SNP weights. Default to equal weights
#' @param  rho weights for the ST
#' @return 
#' \describe{
#'   \item{p.value}{ p-values for AT, S2T, and ST }
#'   \item{pval}{ the list of all p-values }
#'   \item{rho.est}{ estimated optimal \eqn{\rho} value }
#' }
#' @keywords sats
#' @export
#' @references
#' Guo,B. and Wu,B.(2018). Statistical methods to detect novel genetic variants using publicly available GWAS summary data. \emph{CBC}, to appear.
#' @examples
#' R = cor(matrix(rnorm(500),100,5)*sqrt(0.8)+rnorm(100)*sqrt(0.2))
#' Z = rnorm(5) + 0:4
#' sats(Z,R)
sats <- function(Z,R,W=NULL, rho=c((0:5/10)^2,0.5,1)){
  M = length(Z)
  if(is.null(W)) W = rep(1,M)
  ##
  Zw = Z*W; Rw = t(R*W)*W
  eR = eigen(Rw,sym=TRUE);  lamR = abs(eR$val); eta = colSums(eR$vec)*sqrt(lamR)
  R1 = sum(eta^2); R2 = sum(eta^2*lamR)
  c2 = outer(eta,eta)
  Lamq = eigen(diag(lamR) - R2/R1^2*c2, symmetric=TRUE,only.values=TRUE)$val
  ## ST
  Qb = sum(Zw)^2
  pvalb = pchisq(Qb/R1, 1,lower=FALSE)
  ## S2T
  Qv = sum(Zw^2)
  pvalv = KATpval(Qv,lamR)
  ## AT
  L = length(rho)
  if(L<=2){
    return(list(p.value=c(A=NULL, S2=pvalv, S=pvalb), pval=pval) )
  } 
  L1 = L-1; rho1 = rho[-L]
  Qw = (1-rho)*Qv + rho*Qb
  pval = rep(1, L)
  pval[1] = pvalv; pval[L] = pvalb
  Lamk = vector('list', L)
  Lamk[[L]] = R1;  Lamk[[1]] = lamR
  for(k in 2:L1){
    mk = rho[k]*c2;  diag(mk) = diag(mk) + (1-rho[k])*lamR
    aak = zapsmall( abs(eigen(mk, sym=TRUE, only.val=TRUE)$val) )
    Lamk[[k]] = aak[aak>0]
    pval[k] = KATpval(Qw[k],Lamk[[k]])
  }
  minP = min(pval)
  ##
  L = length(rho)
  qval = rep(0,L1)
  for (k in 1:L1) qval[k] = Liu0.qval(minP, Lamk[[k]])
  q1 = qchisq(minP,1,lower=FALSE)
  tauk = (1-rho1)*R2/R1 + rho1*R1
  katint = function(xpar){
    eta1 = sapply(xpar, function(eta0) min((qval-tauk*eta0)/(1-rho1)))
    KATpval(eta1,Lamq)*dchisq(xpar,1)
  }
  prec = 1e-4
  p.value = try({ minP + integrate(katint, 0,q1,  subdivisions=1e3,abs.tol=minP*prec)$val }, silent=TRUE)
  while(class(p.value)=='try-error'){
    prec = prec*2
    p.value = try({ minP + integrate(katint, 0,q1, abs.tol=minP*prec)$val }, silent=TRUE)
  }
  p.value = min(p.value,minP*L)
  return(list(p.value=c(A=p.value, S2=pvalv, S=pvalb), pval=pval, rho.est=rho[which.min(pval)]) )
}



