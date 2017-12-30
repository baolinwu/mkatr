# mkatr
 - An R package implementing various statistical methods for testing variant-set association

------
## SNP-set association tests using GWAS summary data
 - Reference
    - Guo,B. and Wu,B. (2017) Statistical methods to detect novel genetic variants using publicly available GWAS summary data. *tech report*
 - Implemented in the "SATS()" R function returning three test p-values: adaptive test (AT), variance components test (VC), sum test (ST)
 - Sample R codes
 ```r
 library(mkatr)
 R = cor(matrix(rnorm(500),100,5)*sqrt(0.8)+rnorm(100)*sqrt(0.2))
 Z0 = rnorm(5)
 SATS(Z0,R)
 Z = Z0 + c(2,3,-1,2,-3)
 SATS(Z,R)
 Z = Z0 + c(2,3,4,2,3)
 SATS(Z,R)
   ```
