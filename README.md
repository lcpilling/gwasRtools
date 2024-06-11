
# gwasRtools <a href="https://lcpilling.github.io/gwasRtools/"><img src="man/figures/gwasRtools.png" align="right" width="150" /></a>

Some useful R functions for processing GWAS output

<!-- badges: start -->
[![](https://img.shields.io/badge/version-0.1.6.9000-informational.svg)](https://github.com/lcpilling/gwasRtools)
[![](https://img.shields.io/github/last-commit/lcpilling/gwasRtools.svg)](https://github.com/lcpilling/gwasRtools/commits/master)
[![](https://img.shields.io/badge/lifecycle-experimental-orange)](https://www.tidyverse.org/lifecycle/#experimental)
[![DOI](https://zenodo.org/badge/655790727.svg)](https://zenodo.org/badge/latestdoi/655790727)
<!-- badges: end -->

## List of functions
  - [get_loci()](#get_loci)
  - [get_nearest_gene()](#get_nearest_gene)
  - [lambda_gc()](#lambda_gc)


## Installation
To install the development version from GitHub use the `remotes` package:

```r
remotes::install_github("lcpilling/gwasRtools")

remotes::install_github("lcpilling/gwasRtools@*release") # To install the latest release
remotes::install_github("lcpilling/gwasRtools@v0.1.3")   # To install a specific version (see tags)
```


## Example dataset
The package includes a subset of variants from the Graham et al. 2021 GWAS of LDL in 1,320,016 Europeans (GWAS catalog GCST90239658). I will use this throughout.


```r
library(gwasRtools)
head(gwas_example)
#>           SNP CHR        BP    A1 A2      MAF        BETA        SE     P
#> 1 rs370704455   1 104275751 CATCT  C 7.21e-02 -0.00350057 0.0040782 0.391
#> 2 rs545633289   1 104275758     C  G 5.25e-04 -0.36329200 0.2903640 0.211
#> 3 rs367651119   1 104275797     G  A 5.99e-03 -0.01142350 0.0164861 0.488
#> 4 rs560138913   1 104296901     T  G 3.67e-03 -0.05503450 0.0318479 0.084
#> 5 rs533558585   1 104296979     T  C 7.57e-05  0.00173860 0.1826840 0.992
#> 6 rs551694811   1 104297001     T  A 1.97e-03 -0.01075230 0.0316984 0.734
```


## get_loci()
Determine loci from GWAS summary statistics. Use distance from lead significant SNP to estimate independent loci [default distance = 500kb]. Uses -log10(p) derived from BETA/SE so does not need P as input. Example below with default input.

```r
gwas_loci = get_loci(gwas_example)
#> 
#> Locus size (bases) = 5e+05
#> P-value threshold = 5e-08
#> 
#> N variants = 319732
#> N variants p<threshold = 4132
#> N loci = 15

head(gwas_loci)
#>               SNP CHR        BP   A1 A2   MAF   BETA     SE        P locus  lead
#> 57882  rs12046439   1 107536799    T  C 0.248 0.0099 0.0017 5.01e-09     1 FALSE
#> 57900 rs143849791   1 107537916 CATG  C 0.325 0.0128 0.0016 5.85e-15     1 FALSE
#> 57922 rs113329442   1 107539252    A  G 0.330 0.0110 0.0014 1.27e-13     1 FALSE
#> 57987   rs3861909   1 107544176    G  A 0.327 0.0118 0.0015 3.52e-15     1 FALSE
#> 58025  rs17496332   1 107546375    A  G 0.331 0.0111 0.0014 8.70e-14     1 FALSE
#> 58091   rs2878349   1 107549245    G  A 0.327 0.0118 0.0014 2.33e-15     1 FALSE

gwas_loci |> dplyr::filter(lead==TRUE) |> head()
#>           SNP CHR        BP A1 A2   MAF    BETA     SE         P locus lead
#> 1 rs111232683   1 107566149  G  C 0.343  0.0135 0.0016  5.43e-17     1 TRUE
#> 2 rs114254196   1 108635400  C  T 0.008 -0.0448 0.0081  4.38e-08     2 TRUE
#> 3 rs115292790   1 109310728  G  A 0.013 -0.0563 0.0060  2.01e-20     3 TRUE
#> 4  rs12740374   1 109817590  G  T 0.219 -0.1482 0.0016 4.73e-305     4 TRUE
#> 5 rs140266316   1 110326545  G  A 0.016 -0.0577 0.0059  4.73e-22     5 TRUE
#> 6    rs657801   1 111736389  T  C 0.315  0.0090 0.0015  1.88e-09     6 TRUE
```

 - Loci are numbered. Variants within a locus (i.e., significant below the `p_threshold` and less than `n_bases` from last significant variant).
 - Lead variant for each locus is highlighted where `lead==TRUE` (i.e., smallest p-value for any variant within a locus)

The HLA region can be treated as one continuous locus by setting `single_hla_locus` to TRUE. 

### Use LD clumping to identify independent SNPs at the same locus 

Setting option `ld_clump=TRUE` will use {[ieugwasr](https://github.com/MRCIEU/ieugwasr)} package `ld_clump()` function to run Plink LD clumping. 

Default is to use a local Plink installation (this is faster) with EUR reference panel. But setting option `ld_clump_local` to FALSE will use the online IEU API. See the {ieugwasr} docs for details. Default R2 threshold for LD pruning is 0.01 (modify with `ld_pruning_r2` option). 


```r
gwas_loci = get_loci(gwas_example, ld_clump=TRUE)
#> ** Performing LD clumping. Can take a few minutes
#> ** Local Plink installation will be called -- output appears in your R terminal
#> 
#> Locus size (bases) = 5e+05
#> P-value threshold = 5e-08
#> 
#> N variants = 319732
#> N variants p<threshold = 4132
#> N loci = 15
#> N independent variants (LD R2 threshold 0.01) = 153

head(gwas_loci)
#>               SNP CHR        BP   A1 A2   MAF   BETA     SE        P locus  lead lead_dist lead_ld
#> 57882  rs12046439   1 107536799    T  C 0.248 0.0099 0.0017 5.01e-09     1 FALSE     FALSE   FALSE
#> 57900 rs143849791   1 107537916 CATG  C 0.325 0.0128 0.0016 5.85e-15     1 FALSE     FALSE   FALSE
#> 57922 rs113329442   1 107539252    A  G 0.330 0.0110 0.0014 1.27e-13     1 FALSE     FALSE   FALSE
#> 57987   rs3861909   1 107544176    G  A 0.327 0.0118 0.0015 3.52e-15     1 FALSE     FALSE   FALSE
#> 58025  rs17496332   1 107546375    A  G 0.331 0.0111 0.0014 8.70e-14     1 FALSE     FALSE   FALSE
#> 58091   rs2878349   1 107549245    G  A 0.327 0.0118 0.0014 2.33e-15     1 FALSE     FALSE   FALSE

gwas_loci |> dplyr::filter(lead==TRUE) |> head()
#>           SNP CHR        BP A1 A2   MAF    BETA     SE        P locus lead lead_dist lead_ld
#> 1 rs111232683   1 107566149  G  C 0.343  0.0135 0.0016 5.43e-17     1 TRUE      TRUE    TRUE
#> 2 rs114254196   1 108635400  C  T 0.008 -0.0448 0.0081 4.38e-08     2 TRUE      TRUE    TRUE
#> 3 rs140300970   1 109020060  A  T 0.022 -0.0278 0.0049 2.13e-08     3 TRUE     FALSE    TRUE
#> 4 rs148503795   1 109166178  C  G 0.010 -0.0423 0.0070 2.47e-09     3 TRUE     FALSE    TRUE
#> 5  rs74896173   1 109167705  T  C 0.009 -0.0429 0.0076 2.08e-08     3 TRUE     FALSE    TRUE
#> 6 rs111751551   1 109242056  G  A 0.010 -0.0505 0.0069 4.32e-13     3 TRUE     FALSE    TRUE
```


Where before, locus 3 would only have had one lead SNP (based on distance/lowest p-value) `ld_clump()` has identified multiple independent variants in the region.

Note that there are now three `lead` columns:
 - `lead_dist` is the original `lead` column, simply based on distance and p-values
 - `lead_ld` is the direct results from `ld_clump()`
 - The `lead` column combines the two (some SNPs are missing from LD panel so it does not always choose the lowest p-value if only 1 variant identified at a locus).

** Note that `ld_clump()` only considers R^2 when defining independent variants, not D' -- you should perform additional checking/conditional analysis where relevant for `lead` variants in close proximity.


## get_nearest_gene()
Get nearest gene from a set of variants using GENCODE data. Need to provide a data.frame of variant IDs (e.g., rsids), CHR and POS. Default column names are the same as for `get_loci()`. Default max distance from variant to gene is 100kb.


```r
gwas_loci = get_nearest_gene(gwas_loci, build=37)
#> Using human genome build 37
#> Getting nearest gene for 4131 unique variants
#> (Removed 1 duplicated or missing variant IDs/positions)

head(gwas_loci)
#>           SNP CHR        BP   A1 A2   MAF   BETA     SE        P locus  lead  gene  dist
#> 1  rs12046439   1 107536799    T  C 0.248 0.0099 0.0017 5.01e-09     1 FALSE PRMT6 62468
#> 2 rs143849791   1 107537916 CATG  C 0.325 0.0128 0.0016 5.85e-15     1 FALSE PRMT6 61351
#> 3 rs113329442   1 107539252    A  G 0.330 0.0110 0.0014 1.27e-13     1 FALSE PRMT6 60015
#> 4   rs3861909   1 107544176    G  A 0.327 0.0118 0.0015 3.52e-15     1 FALSE PRMT6 55091
#> 5  rs17496332   1 107546375    A  G 0.331 0.0111 0.0014 8.70e-14     1 FALSE PRMT6 52892
#> 6   rs2878349   1 107549245    G  A 0.327 0.0118 0.0014 2.33e-15     1 FALSE PRMT6 50022

gwas_loci |> dplyr::filter(lead==TRUE) |> head()
#>           SNP CHR        BP A1 A2   MAF    BETA     SE        P locus lead     gene   dist
#> 1 rs111232683   1 107566149  G  C 0.343  0.0135 0.0016 5.43e-17     1 TRUE    PRMT6  33118
#> 2 rs114254196   1 108635400  C  T 0.008 -0.0448 0.0081 4.38e-08     2 TRUE SLC25A24  41258
#> 3 rs140300970   1 109020060  A  T 0.022 -0.0278 0.0049 2.13e-08     3 TRUE    NBPF6   6436
#> 4 rs148503795   1 109166178  C  G 0.010 -0.0423 0.0070 2.47e-09     3 TRUE  FAM102B -63467
#> 5  rs74896173   1 109167705  T  C 0.009 -0.0429 0.0076 2.08e-08     3 TRUE  FAM102B -64994
#> 6 rs111751551   1 109242056  G  A 0.010 -0.0505 0.0069 4.32e-13     3 TRUE  PRPF38B  -7111
```
 - If `dist` is positive, the variant is intergenic, and this is the distance to the closest gene.
 - If `dist` is negative, the variant is within a gene, and this is the distance to the start of the gene.
 - If `dist` is NA, the variant is not within `n_bases` of a gene in GENCODE.



## lambda_gc()
Estimate inflation of test statistics. Lambda GC compares the median test statistic against the expected median test statistic under the null hypothesis of no association. For well-powered GWAS of traits with a known polygenic inheritance, we expect inflation of lambda GC. For traits with no expected association, we expect lambda GC to be around 1.


```r
lambda_gc(gwas_example$P)
#> [1] 1.41112
```
