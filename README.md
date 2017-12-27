# mkatr
 - An R pacakge implementing various statistical methods for testing variant-set association

------
## On sequence-kernel association test of rare variant set
 - Reference
    - Wu,B., Pankow,J.S., Guan,W. (2015) Sequence kernel association analysis of rare variant set based on the marginal regression model for binary traits. *Genetic Epidemiology*, 39(6), 399-405.
    - Wu,B., Guan,W., Pankow,J.S. (2016) On efficient and accurate calculation of significance p-values for sequence kernel association test of variant set. *Annals of human genetics*, 80(2), 123-135.
 
 - Sample R codes
```r
library(CompQuadForm)
library(SKAT)
library(mkatr)
## simulate outcomes
X = matrix(rnorm(10000), 5000,2)
D = rbinom(5000,1,0.5)
G = matrix(rbinom(100000,2,0.01), 5000,10)
## SKATL and SKATOL tests
SKATL(KAT.null(D,X), G, c(1.5,25.5))
SKATOL(KAT.null(D,X), G, c(1.5,25.5))
## compared to SKAT package
SKAT(G, SKAT_Null_Model(D ~ X, out_type="D"))$p.value
SKAT(G, SKAT_Null_Model(D ~ X, out_type="D"), method="optimal.adj")$p.value
```

------
## SNP-set association tests using GWAS summary data
 - Reference
    - Guo,B. and Wu,B. (2017) Statistical methods to detect novel genetic variants using publicly available GWAS summary data. *tech report*
 - Sample R codes
 ```r
 library(mkatr)
 R = cor(matrix(rnorm(500),100,5)*sqrt(0.8)+rnorm(100)*sqrt(0.2))
 Z = rnorm(5) + 0:4
 ASATZ(Z,R)
 ```
