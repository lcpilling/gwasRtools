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
Determine loci from a GWAS summary statistics file. Use distance from lead significant SNP to estimate independet loci in GWAS summary stats. Uses -log10(p) derived from BETA/SE so does not need P as input. Example below with default input (output truncated):

```r
get_loci(
  gwas,
  snp_col = "SNP",
  chr_col = "CHR",
  pos_col = "BP",
  maf_col = "MAF",
  beta_col = "BETA",
  se_col = "SE",
  n_bases = 5e5,
  p_threshold = 5e-8
)

# example using BOLT-LMM output:
gwas_loci = get_loci(gwas, maf_col="A1FREQ")

head(gwas_loci, 5)
#> # A tibble: 5 × 11
#>    SNP               CHR        BP ALLELE1 ALLELE0 A1FREQ    BETA      SE  P_BOLT_LMM locus lead 
#>    <chr>           <dbl>     <dbl> <chr>   <chr>    <dbl>   <dbl>   <dbl>       <dbl> <dbl> <lgl>
#>  1 rs3041466           2 142820625 C       CAA      0.492 -0.0177 0.00317 0.000000021     1 TRUE 
#>  2 rs5858140           4  49039586 C       CA       0.429 -0.0182 0.00332 0.000000036     2 FALSE
#>  3 rs2605231           4  49057424 A       T        0.293 -0.0179 0.00329 0.000000044     2 FALSE
#>  4 4:49057608_AT_A     4  49057608 AT      A        0.301 -0.0188 0.00338 0.000000024     2 FALSE
#>  5 4:49057930_CT_C     4  49057930 CT      C        0.286 -0.0186 0.00336 0.000000029     2 FALSE

head(gwas_loci |> filter(lead==TRUE), 5)
#> # A tibble: 5 × 11
#>    SNP          CHR        BP ALLELE1 ALLELE0 A1FREQ    BETA      SE P_BOLT_LMM locus lead 
#>    <chr>      <dbl>     <dbl> <chr>   <chr>    <dbl>   <dbl>   <dbl>      <dbl> <dbl> <lgl>
#>  1 rs3041466      2 142820625 C       CAA      0.492 -0.0177 0.00317   2.10e- 8     1 TRUE 
#>  2 rs11722559     4  49227587 T       C        0.225 -0.0240 0.00391   7.70e-10     2 TRUE 
#>  3 rs4865414      4  52758294 C       T        0.684  0.0200 0.00322   5.10e-10     3 TRUE 
#>  4 rs74405522     6 157205286 T       C        0.962  0.0432 0.00788   4.20e- 8     4 TRUE 
#>  5 rs55730499     6 161005610 C       T        0.918  0.0378 0.00548   5.80e-12     5 TRUE
```

 - Loci are numbered. Variants within a locus (i.e., significant below the `p_threshold` and less than `n_bases` from last significant variant).
 - Lead variant for each locus is highlighted where `lead==TRUE` (i.e., smallest p-value for any variant within a locus)

## get_nearest_gene()
Get nearest gene from a set of variants using GENCODE data. Need to provide a data.frame of variant IDs (e.g., rsids), CHR and POS. Defaults below, with example output:

```r
get_nearest_gene(
  snps,
  snp_col = "SNP",
  chr_col = "CHR",
  pos_col = "BP",
  build   = 37,
  n_bases = 1e5
)
#> # A tibble: 3 × 5
#>   SNP          CHR        BP gene        dist
#>   <chr>      <dbl>     <dbl> <chr>      <dbl>
#> 1 rs55730499     6 161005610 LPA       -53095
#> 3 rs814573      19  45424351 APOC1       1745
#> 3 rs123456      20  98765432 NA            NA
```
 - If `dist` is positive, the variant is intergenic, and this is the distance to the closest gene.
 - If `dist` is negative, the variant is within a gene, and this is the distance to the start of the gene.
 - If `dist` is NA, the variant is not within `n_bases` of a gene in GENCODE.

