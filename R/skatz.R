#' SNP-set association tests using GWAS summary Z-statistics
#'
#' Compute p-values for the SNP-set tests (VC, burden test, adaptive test) using GWAS Z-statistics.
#' With less than two values of rho's, it outputs only p-values for VC and BT.
#' @param  Z summary Z-statistics for a set of SNPs from GWAS
#' @param  R SNP pairwise LD matrix
#' @param  W SNP weights. Default to equal weights
#' @param  rho weights for burden test
#' @return 
#' \describe{
#'   \item{p.value}{ p-values for Adaptive test, VC, and BT }
#'   \item{pval}{ the list of all p-values }
#'   \item{rho.est}{ estimated optimal \eqn{\rho} value }
#' }
#' @keywords ASATZ
#' @export
#' @references
#' Guo,B., Wu,B.(2017). Statistical methods to detect novel genetic variants using publicly available GWAS summary data. tech report.
#' @examples
#' R = cor(matrix(rnorm(500),100,5)*sqrt(0.8)+rnorm(100)*sqrt(0.2))
#' Z = rnorm(5) + 0:4
#' ASATZ(Z,R)
ASATZ <- function(Z,R,W=NULL, rho=c((0:5/10)^2,0.5,1)){
  M = length(Z)
  if(is.null(W)) W = rep(1,M)
  ##
  Zw = Z*W; Rw = t(R*W)*W
  eR = eigen(Rw,sym=TRUE);  lamR = abs(eR$val); eta = colSums(eR$vec)*sqrt(lamR)
  R1 = sum(eta^2); R2 = sum(eta^2*lamR)
  c2 = outer(eta,eta)
  Lamq = eigen(diag(lamR) - R2/R1^2*c2, symmetric=TRUE,only.values=TRUE)$val
  ## BT
  Qb = sum(Zw)^2
  pvalb = pchisq(Qb/R1, 1,lower=FALSE)
  ## VC 
  Qv = sum(Zw^2)
  pvalv = KATpval(Qv,lamR)
  ## A
  L = length(rho)
  if(L<=2){
    return(list(p.value=c(A=NULL, V=pvalv, B=pvalb), pval=pval) )
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
  qval = rep(0,L1)
  for(k in 1:L1) qval[k] = Liu0.qval(minP, Lamk[[k]])
  Rs = rowSums(Rw); R3 = sum(Rs*colSums(Rw*Rs))
  lam = eigen(Rw-outer(Rs,Rs)/R1,sym=TRUE,only.val=TRUE)$val
  tauk = (1-rho1)*R2/R1 + rho1*R1;  vp2 = 4*(R3/R1-R2^2/R1^2)
  MuQ = sum(lam);  VarQ = sum(lam^2)*2
  sd1 = sqrt(VarQ)/sqrt(VarQ+vp2)
  q1 = qchisq(minP,1,lower=FALSE)
  katint = function(xpar){
    eta1 = sapply(xpar, function(eta0) min((qval-tauk*eta0)/(1-rho1)))
    x = (eta1-MuQ)*sd1 + MuQ
    KATpval(x,lam)*dchisq(xpar,1)
  }
  prec = 1e-4
  p.value = try({ minP + integrate(katint, 0,q1,  subdivisions=1e3,abs.tol=minP*prec)$val }, silent=TRUE)
  while(class(p.value)=='try-error'){
    p.value = try({ minP + integrate(katint, 0,q1, abs.tol=minP*prec)$val }, silent=TRUE)
    prec = prec*2
  }
  return(list(p.value=c(A=p.value, V=pvalv, B=pvalb), pval=pval, rho.est=rho[which.min(pval)]) )
}

