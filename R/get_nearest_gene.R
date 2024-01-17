#' Get nearest gene from a set of variants
#'
#' @description Use distance variant position to get the nearest gene. Uses {snpsettest} function `map_snp_to_gene` to identify genes from GENCODE databases (https://github.com/HimesGroup/snpsettest).
#'
#' @return Returns a data frame of variant IDs mapped to genes (with distance). 
#'
#' - If `dist` is positive, the variant is intergenic, and this is the distance to the closest gene.
#' - If `dist` is negative, the variant is within a gene, and this is the distance to the start of the gene.
#' - If `dist` is NA, the variant is not within `n_bases` of a gene in GENCODE.
#'
#' @author Luke Pilling
#'
#' @name get_nearest_gene
#'
#' @param variants A data.frame. Contains the variants (e.g., in summary statistics).
#' @param detect_headers Logical. Default=TRUE. Search input headers to see if BOLT-LMM, SAIGE, or REGENIE input (user therefore doesn't need to provide). 
#' @param snp_col A string. Default="SNP". The RSID/variantID column name.
#' @param chr_col A string. Default="CHR". The chromosome column name.
#' @param pos_col A string. Default="BP". The base pair/position column name.
#' @param build An integer. Default=37. Genome build to use (can only be 37 or 38).
#' @param n_bases An interger. Default=1e5. The max distance in base-pairs between a variant and a gene to annotate
#'
#' @examples
#' gwas_loci = get_loci(gwas_example)
#'
#' gwas_loci_genes = get_nearest_gene(gwas_loci)
#'
#' head(gwas_loci_genes)
#'
#' head(gwas_loci_genes[ gwas_loci_genes$lead==TRUE , ])
#'
#' @export
#'

# function to get nearest gene from mapped output
get_nearest_gene = function(variants,
                            detect_headers  = TRUE,
                            snp_col = "SNP",
                            chr_col = "CHR",
                            pos_col = "BP",
                            build   = 37,
                            n_bases = 1e5)  {
	
	# check build input
	if (! build %in% c(37,38))  stop("Build has to be 37 or 38")
	
	# check headers. Default is BOLT-LMM. Is this SAIGE, or REGENIE output? If not, user needs to specify
	col_names = colnames(variants)
	if (detect_headers)  {
		if ("SNPID" %in% col_names & "CHR" %in% col_names & "POS" %in% col_names)  {
			cat("Detected SAIGE input. Using default headers. Disable with `detect_headers=FALSE`\n\n")
			snp_col  = "SNPID"
			chr_col  = "CHR"
			pos_col  = "POS"
		}
		if ("ID" %in% col_names & "CHROM" %in% col_names & "GENPOS" %in% col_names)  {
			cat("Detected REGENIE input. Using default headers. Disable with `detect_headers=FALSE`\n\n")
			snp_col  = "ID"
			chr_col  = "CHROM"
			pos_col  = "GENPOS"
		}
	}
	
	# check headers
	if (! snp_col %in% col_names) stop(paste0("`snp_col` \"", snp_col, "\" not in provided data frame"))
	if (! chr_col %in% col_names) stop(paste0("`chr_col` \"", chr_col, "\" not in provided data frame"))
	if (! pos_col %in% col_names) stop(paste0("`pos_col` \"", pos_col, "\" not in provided data frame"))
	
	# get gene list
	gene_curated = snpsettest::gene.curated.GRCh37
	if (build == 38)  gene_curated = snpsettest::gene.curated.GRCh38
	
	# messages
	cat(paste0("Using human genome build ", build, "\n"))
	
	# get data frame of variant IDs and positions - remove duplicates
	variants_map = variants |> 
		dplyr::mutate(id=!! rlang::sym(snp_col), chr=!! rlang::sym(chr_col), pos=!! rlang::sym(pos_col)) |> 
		dplyr::select(id, chr, pos) |> 
		as.data.frame() |> 
		dplyr::distinct() |> 
		na.omit()
	cat(paste0("Getting nearest gene for ", nrow(variants_map), " unique variants\n"))
	if (nrow(variants)>nrow(variants_map))  cat(paste0("(Removed ", nrow(variants)-nrow(variants_map), " duplicated or missing variant IDs/positions)\n"))
	
	# get GENCODE genes from list of variants
	variants_map = snpsettest::map_snp_to_gene(variants_map, gene_curated, extend_start=n_bases/1000, extend_end=n_bases/1000)
	variants_map = variants_map$map
	
	# merge gene symbols to IDs
	map = dplyr::left_join(variants_map, gene_curated |> dplyr::select(gene.id, gene.name), by="gene.id")
	
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
	
	# use `map` from {purrr} to make this quick and easy
	map2 = purrr::map(unique(map$id), \(id) gwasRtools:::get_nearest_gene1(id, map)) |> 
		purrr::list_rbind() |> 
		dplyr::select(id, gene, dist)
	
	# merge original variants file with genes 
	variants = variants |> dplyr::mutate(ID=!! rlang::sym(snp_col))
	variants = dplyr::left_join(variants, map2, by=c("ID"="id"))
	variants = variants |> dplyr::select(-ID)
	
	# return 
	return(variants)
}



#' Internal function to get nearest gene for each unique variant in `map`
#' @param id Identifier for the variant
#' @param map Modified output from snpsettest::map_snp_to_gene()
#' @noRd
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
			jj = which(map[map[,"id"]==id,"dist"] == min(map[map[,"id"]==id,"dist"]))[1]
			gene = map[map[,"id"]==id,"gene.name"][ jj ]
			dist = map[map[,"id"]==id,"dist"][ jj ]
		
		}
	}
	
	if (is.na(gene)) dist = NA 
	
	mapped = data.frame(id, gene, dist)
	return(mapped)
	
}

