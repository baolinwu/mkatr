# mkatr
 - An R package implementing various statistical methods for testing variant-set association

------
## On sequence-kernel association test of rare variant set
 - References
    - Wu,B., Pankow,J.S., Guan,W. (2015) Sequence kernel association analysis of rare variant set based on the marginal regression model for binary traits. *Genetic Epidemiology*, 39(6), 399-405.
    - Wu,B., Guan,W., Pankow,J.S. (2016) On efficient and accurate calculation of significance p-values for sequence kernel association test of variant set. *Annals of human genetics*, 80(2), 123-135.
 - Sample R codes
```r
library(mkatr)
### simulate outcomes
X = matrix(rnorm(10000), 5000,2)
D = rbinom(5000,1,0.5)
G = matrix(rbinom(100000,2,0.01), 5000,10)
### SKATL and SKATOL tests
SKATL(KAT.null(D,X), G, c(1.5,25.5))
SKATOL(KAT.null(D,X), G, c(1.5,25.5))
### compared to SKAT package
## library(SKAT)
## SKAT(G, SKAT_Null_Model(D ~ X, out_type="D"))$p.value
## SKAT(G, SKAT_Null_Model(D ~ X, out_type="D"), method="optimal.adj")$p.value
```

------
## SNP-set association tests using GWAS summary data
 - Reference
    - Guo,B. and Wu,B. (2018) Statistical methods to detect novel genetic variants using publicly available GWAS summary data. *CBC*, to appear.
 - Implemented in the "sats()" R function returning three test p-values: adaptive test (AT), squared sum test (S2T), sum test (ST)
 - Sample R codes
 ```r
 library(mkatr)
 R = cor(matrix(rnorm(500),100,5)*sqrt(0.8)+rnorm(100)*sqrt(0.2))
 Z0 = rnorm(5)
 sats(Z0,R)
 Z = Z0 + c(2,3,-1,2,-3)
 sats(Z,R)
 Z = Z0 + c(2,3,4,2,3)
 sats(Z,R)
   ```
