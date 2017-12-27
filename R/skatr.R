#' Reproduce the p-value from the SKAT R package
#' @export
SKAT.pval <- function(Q.all, lambda){
  pval = rep(0, length(Q.all))
  i1 = which(is.finite(Q.all))
  for(i in i1){
    tmp = davies(Q.all[i],lambda,acc=1e-6,lim=1e4); pval[i] = tmp$Qq
    if((tmp$ifault>0)|(pval[i]<=0)|(pval[i]>=1)) pval[i] = Liu.pval(Q.all[i],lambda)
  }
  return(pval)
}
####
#' Sequence kernel association test (SKAT) using variant test statistics
#'
#' Reproduce the SKAT p-value (from the SKAT R package) using marginal variant score statistics
#' 
#' We compute the SKAT based on the variant test statistics (typically of much smaller dimension), which
#' leads to much more efficient computations. The same algorithm as the SKAT R package
#' is used to compute the tail probability of 1-DF chi-square mixtures.
#'
#' @param  obj a fitted null model using KAT.null() or KAT.cnull()
#' @param  G genotype matrix, sample in rows, variant in columns
#' @param  W.beta Beta parameters for variant weights
#' @return SKAT p-value
#' @keywords SKAT
#' @export
#' @references
#' Wu, M.C., Lee, S., Cai, T., Li, Y., Boehnke, M., and Lin, X. (2011) Rare Variant Association Testing for Sequencing Data Using the Sequence Kernel Association Test (SKAT). American Journal of Human Genetics, 89, 82-93.
#'
#' Wu,B., Guan,W., and Pankow,J.S. (2016) On efficient and accurate calculation of significance p-values for sequence kernel association test of variant set. Annals of Human Genetics, 80(2), 123-135.
SKATv <- function(obj,G, W.beta=c(1,25)){
  N = dim(G)[2]; maf = colMeans(G)/2
  W = maf^(W.beta[1]-1)*(1-maf)^(W.beta[2]-1);  W = W/sum(W)*N
  tmp = t(obj$Ux)%*%G
  if(obj$mode=='C'){
    Gs = t(G)%*%G - t(tmp)%*%tmp
    Zs = colSums(obj$U0*G)/sqrt(obj$s2)
  } else{
    Gs = t(G*obj$Yv)%*%G - t(tmp)%*%tmp
    Zs = colSums(obj$U0*G)
  }
  R = t(Gs*W)*W
  Z = Zs*W
  lam = eigen(R, sym=TRUE,only.val=TRUE)$val
  SKAT.pval(sum(Z^2), lam)
}
####
#' Sequence kernel association test (SKAT) with linear kernel using variant test statistics
#'
#' Compute accurate SKAT (linear kernel) p-value based on marginal variant score statistics
#'
#' We compute the SKAT based on the variant test statistics (typically of much smaller dimension), which
#' leads to much more efficient computations. Davies' method is used to compute the tail probability
#' of 1-DF chi-square mixtures with more stringent convergence criteria (acc=1e-9,lim=1e6). When it
#' fails, we then switch to the saddlepoint approximation.
#'
#' @param  obj a fitted null model using KAT.null() or KAT.cnull()
#' @param  G genotype matrix, sample in rows, variant in columns
#' @param  W.beta Beta parameters for variant weights
#' @return more accurate SKAT p-value
#' @keywords SKAT
#' @export
#' @references
#' Wu, M.C., Lee, S., Cai, T., Li, Y., Boehnke, M., and Lin, X. (2011) Rare Variant Association Testing for Sequencing Data Using the Sequence Kernel Association Test (SKAT). American Journal of Human Genetics, 89, 82-93.
#'
#' Wu,B., Guan,W., and Pankow,J.S. (2016) On efficient and accurate calculation of significance p-values for sequence kernel association test of variant set. Annals of Human Genetics, 80(2), 123-135.
SKATh <- function(obj,G, W.beta=c(1,25)){
  N = dim(G)[2]; maf = colMeans(G)/2
  W = maf^(W.beta[1]-1)*(1-maf)^(W.beta[2]-1);  W = W/sum(W)*N
  tmp = t(obj$Ux)%*%G
  if(obj$mode=='C'){
    Gs = t(G)%*%G - t(tmp)%*%tmp
    Zs = colSums(obj$U0*G)/sqrt(obj$s2)
  } else{
    Gs = t(G*obj$Yv)%*%G - t(tmp)%*%tmp
    Zs = colSums(obj$U0*G)
  }
  R = t(Gs*W)*W
  Z = Zs*W
  lam = eigen(R, sym=TRUE,only.val=TRUE)$val
  KATpval(sum(Z^2), lam)
}
