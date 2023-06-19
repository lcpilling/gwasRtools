#' Get nearest gene from a set of variants
#'
#' @description Use distance variant position to get the nearest gene. Uses {snpsettest} function `map_snp_to_gene` to identify genes from GENCODE databases (https://github.com/HimesGroup/snpsettest).
#'
#' @return Returns a data frame of variant IDs mapped to genes (with distance)
#'
#' @author Luke Pilling
#'
#' @name get_nearest_gene
#'
#' @param gwas A data.frame. Contains the GWAS summary statistics.
#' @param snp_col A string. Default="SNP". The RSID/variantID column name.
#' @param chr_col A string. Default="CHR". The chromosome column name.
#' @param pos_col A string. Default="BP". The base pair/position column name.
#' @param build An integer. Default=37. Genome build to use (can only be 37 or 38).
#' @param n_bases An interger. Default=1e5. The max distance in base-pairs between a variant and a gene to annotate
#'
#' @examples
#' get_nearest_gene(gwas)
#'
#' @export
#'

# function to get nearest gene from mapped output
get_nearest_gene = function(gwas,
                            snp_col = "SNP",
                            chr_col = "CHR",
                            pos_col = "BP",
                            build   = 37,
                            n_bases = 1e5)  {
	
	# check input
	if (! build %in% c(37,38))  stop("Build has to be 37 or 38")
	if (! snp_col %in% colnames(gwas)) stop(paste0("`snp_col` \"", snp_col, "\" not in provided data frame"))
	if (! chr_col %in% colnames(gwas)) stop(paste0("`chr_col` \"", chr_col, "\" not in provided data frame"))
	if (! pos_col %in% colnames(gwas)) stop(paste0("`pos_col` \"", pos_col, "\" not in provided data frame"))
	
	# get gene list
	gene_curated = snpsettest::gene.curated.GRCh37
	if (build == 38)  gene_curated = snpsettest::gene.curated.GRCh38
	
	# get SNPs actually in the genes
	gwas_map = gwas |> dplyr::mutate(id=!! rlang::sym(snp_col), chr=!! rlang::sym(chr_col), pos=!! rlang::sym(pos_col)) |> as.data.frame()
	gwas_mapped = snpsettest::map_snp_to_gene(gwas_map, gene_curated, extend_start=n_bases/1000, extend_end=n_bases/1000)
	gwas_mapped = gwas_mapped$map
	
	# merge gene symbols to IDs
	map = dplyr::left_join(gwas_mapped, gene_curated |> dplyr::select(gene.id, gene.name), by="gene.id")
	
	# get unique pos & distances
	map = map |> dplyr::mutate(
		dist_start = gene.start - pos, # positive integer = before start
		dist_end   = gene.end - pos,   # positive integer = before end
		dist       = NA)               # if variant is outside gene (both positive) update with distance
	map = map |> dplyr::mutate(
		dist = dplyr::case_when(
			dist_start > 0  & dist_end >  0 & dist_start > dist_end ~ dist_end,
			dist_start > 0  & dist_end >  0 & dist_start < dist_end ~ dist_start,
			dist_start < 0  & dist_end <  0 & dist_start < dist_end ~ abs(dist_end),
			dist_start < 0  & dist_end <  0 & dist_start > dist_end ~ abs(dist_start),
			dist_start <= 0 & dist_end >= 0 ~ 0, # in gene
			TRUE ~ NA_real_))
	
	# function to get nearest gene for each unique variant in `map`
	get_nearest_gene1 = function(id, map) {
		
		gene = NA 
		dist = 0
		
		# any gene within distance of variant?
		if (any(!is.na(map[map[,"id"]==id,"gene.name"])))  {
		
			# if variant in a gene
			if (any( map[map[,"id"]==id,"dist"] == 0 ) ) {
			
				# get first gene name where variant is in the gene
				gene = map[map[,"id"]==id,"gene.name"][ map[map[,"id"]==id,"dist"] == 0 ][1]
				
				# get distance to start (negative)
				dist = map[map[,"id"]==id,"dist_start"][ map[map[,"id"]==id,"dist"] == 0 ][1]
			
			} else {
			
				# closest to start OR end
				jj = which(map[map[,"id"]==id,"dist"] == min(map[map[,"id"]==id,"dist"]))
				gene = map[map[,"id"]==id,"gene.name"][ jj ]
				dist = map[map[,"id"]==id,"dist"][ jj ]
			
			}
		}
		
		if (is.na(gene)) dist = NA 
		
		mapped = data.frame(id, gene, dist)
		return(mapped)
		
	}
	
	# use `map` from {purrr} to make this quick and easy
	map2 = purrr::map(unique(map$id), \(id) get_nearest_gene1(id, map)) |> purrr::list_rbind()
	
	# merge original GWAS file with genes 
	gwas = gwas |> dplyr::mutate(ID=!! rlang::sym(snp_col))
	gwas_genes = dplyr::left_join(gwas, map2 |> dplyr::select(id, gene, dist), by=c("ID"="id"))
	gwas_genes = gwas_genes |> dplyr::select(-ID)
	
	# return 
	return(gwas_genes)
}



