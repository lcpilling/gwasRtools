# gwasRtools
Some useful R functions for processing GWAS output

[![](https://img.shields.io/badge/version-0.1.0-informational.svg)](https://github.com/lukepilling/gwasRtools)
[![](https://img.shields.io/github/last-commit/lukepilling/gwasRtools.svg)](https://github.com/lukepilling/gwasRtools/commits/master)
[![](https://img.shields.io/badge/lifecycle-experimental-orange)](https://www.tidyverse.org/lifecycle/#experimental)

## List of functions
  - [lambda_gc()](#lambda_gc)
  - [get_loci()](#get_loci)
  - [get_nearest_gene()](#get_nearest_gene)



## lambda_gc()
Estimate inflation of test statistics. Lambda GC compares the median test statistic against the expected median test statistic under the null hypothesis of no association. For well-powered quantitative traits with a known polygenic inheritance, we expect inflation of lambda GC, but for traits with no expected association, we expect lambda GC to be around 1.

```r
lambda_gc(p_values)
```



## get_loci()
Determine loci from a GWAS summary statistics file. Use distance from lead significant SNP to estimate independet loci in GWAS summary stats. Uses -log10(p) derived from BETA/SE so does not need P as input.

```r
get_loci(
  gwas,
  snp_col = "SNP",
  chr_col = "CHR",
  pos_col = "BP",
  maf_col = "MAF",
  beta_col = "BETA",
  se_col = "SE",
  n_bases = 1e+06
)
```



## get_nearest_gene()


