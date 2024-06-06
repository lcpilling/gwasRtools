# gwasRtools 0.1.6.9000 (6 June 2024)

### Changes
* Change URL to reflect my GitHub username change from `lukepilling` to `lcpilling` to be more consistent between different logins, websites, and social media
 - https://lcpilling.github.io/gwasRtools
 - https://github.com/lcpilling/gwasRtools


# gwasRtools 0.1.6 (20 April 2024)

* `get_loci()` -- update errors and use p-value if provided


# gwasRtools 0.1.5 (6 March 2024)

* Change how functions detect SAIGE input
* Added warnings for depreciated options in `get_loci()`


# gwasRtools 0.1.4 (17 Jan 2024)

## Updates
* `get_loci()` and `get_nearest_gene()` detect BOLT-LMM, SAIGE and REGENIE input automatically. Can be disabled with `detect_headers=FALSE`

## Breaking changes 
* `get_loci()` option `get_ld_indep` is renamed `ld_clump` to be more consistent with the {ieugwasr} package (https://github.com/MRCIEU/ieugwasr)
* `get_loci()` option `ld_pruning_r2` is renamed `ld_clump_r2` to be more consistent with the {ieugwasr} package (https://github.com/MRCIEU/ieugwasr)
* `get_loci()` option `exclude_hla` is renamed `single_hla_locus` because we are not actually excluding it, just treating it as one continuous locus


# gwasRtools 0.1.3 (27 Nov 2023)

* `get_loci()` now has the option `exclude_hla` to treat the HLA region as one continuous locus (default = FALSE). Coordinates are specified with `hla_pos = c(25e6, 34e6)` [these are the defaults]
* Minor documentation updates


# gwasRtools 0.1.2 (15 Sept 2023)

* Added `gwas_example` data for running examples
* `get_loci()` Genetic loci are now defined around variants with p-value < 5 × 10−8. The locus borders were set 500kb (`n_bases`) to each side of the highest genome-wide significant variant in each region. 
* `get_loci()` "lead" column is the best lead SNP from distance (i.e., all input SNPs) or from LD clumping (subset in reference dataset panel)
* Fixes for `get_loci()` niche cases
* Fixes for `get_nearest_gene()` niche cases


# gwasRtools 0.1.1 (30 Aug 2023)

* Fix internal -log10p calculation
* Swap most "Depends" to "Imports"
* Add dependency for {ieugwasr}
* Add `ld_clump` option to `get_loci()` function, to get independent SNPs at the same locus
* Add `use_pvalue` option to `get_loci()` function


# gwasRtools 0.1.0 (30 June 2023)

* Added a `NEWS.md` file to track changes to the package.
