# mkatr
 - An R package implementing various statistical methods for testing variant-set association

------
## SNP-set association tests using GWAS summary data
 - Reference
    - Guo,B. and Wu,B. (2017) Statistical methods to detect novel genetic variants using publicly available GWAS summary data. *tech report*
 - Sample R codes
 ```r
 library(mkatr)
 R = cor(matrix(rnorm(500),100,5)*sqrt(0.8)+rnorm(100)*sqrt(0.2))
 Z = rnorm(5) + 0:4
 SATS(Z,R)
 ```
