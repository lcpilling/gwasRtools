# gwasRtools 0.1.2

* `get_loci()` "lead" column is the best lead SNP from distance (i.e., all input SNPs) or from LD clumping (subset in EUR panel)
* Fixes for `get_loci()` niche cases
* Fixes for `get_nearest_gene()` niche cases

# gwasRtools 0.1.1

* Fix internal -log10p calculation
* Swap most "Depends" to "Imports"
* Add dependency for {ieugwasr}
* Add `ld_clump` option to `get_loci()` function, to get independent SNPs at the same locus
* Add `use_pvalue` option to `get_loci()` function

# gwasRtools 0.1.0

* Added a `NEWS.md` file to track changes to the package.
