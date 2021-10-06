# Exploration of saddlepoint approximations in binary GWAS

Detailed analysis of saddlepoint approximation of the score test statistic in a GWAS with a binary reponse
such as presence or absence of a disease.

## Usage

Perform the GWAS using normal approximation, single saddlepoint approximation of the efficient score statistic
with second continuity correction, or double saddlepoint approximation with second continuity correction. For the saddlepoint approximations, choose
either using the CDF approximation from Lugananni-Rice of from Barndorff-Nielsen.

Given genotype matrix genos with each variant along the columns, covariate matrix cov with covariates along the columns, and a significance level 
alpha where the null hypothesis is rejected.
Choose methods which is either normal approximation "noraprx", single saddlepoint of efficient score "SPASCC", or
double saddlepoint with second contiuity correction "DSPASCC". For a fast version of double saddlepoint with second continuity correction,
use "DSPASCC_FAST". Compute p-value using Lugananni-Rice "LR" or Barndorff-Nielsen "BN".
For instance, with double saddlepoint approximation with second continuity correction using Barndorff-Nielsen p-value computation: 
```R
p.values = ScoreTest(genos,pheno,cov,alpha=5*10^-8,pval.comp = "BN",methods = "DSPASCC")

```
