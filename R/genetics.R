#' Estimate Lambda GC from p-values
#'
#' @description Estimate inflation of test statistics. Lambda GC compares the median test statistic against the expected median test statistic under the null hypothesis of no association. For well-powered quantitative traits with a known polygenic inheritance, we expect inflation of lambda GC, but for traits with no expected association, we expect lambda GC to be around 1.
#'
#' Original function by Jing Hua Zhao (https://jinghuazhao.github.io/)
#'
#' @return Returns a vector
#'
#' @author Jing Hua Zhao
#'
#' @name lambda_gc
#'
#' @param p a numeric vector of p-values
#'
#' @examples
#' lambda_gc(p)
#'
#' @export
#'

lambda_gc = function(p) {
	if (!is.numeric(p))  stop("p needs to be numeric")
	
	p = p[!is.na(p)]
	n = length(p)
	
	if (any(p > 1) | any(p < 0))  stop("These don't look like p-values")
	
	observed = qchisq(p,1,lower.tail=FALSE)
	expected = qchisq(1:n/n,1,lower.tail=FALSE)
	
	lambda = median(observed)/median(expected)
	return(lambda)
}


#' Get loci from GWAS file
#'
#' @description Use distance from lead significant SNP to estimate independet loci in GWAS summary stats. Uses -log10(p) derived from BETA/SE so does not need P as input.
#'
#' @return Returns a data frame of significant SNPs with locus annotation
#'
#' @author Luke Pilling
#'
#' @name get_loci
#'
#' @param gwas A data.frame. Contains the GWAS summary statistics.
#' @param snp_col A string. Default="SNP". The RSID/variantID column name.
#' @param chr_col A string. Default="CHR". The chromosome column name.
#' @param pos_col A string. Default="BP". The base pair/position column name.
#' @param maf_col A string. Default="MAF". The MAF/minor-allele-frequency column name.
#' @param beta_col A string. Default="BETA". The BETA column name.
#' @param se_col A string. Default="SE". The SE column name.
#' @param stat_col A string. Default="NA". The test statistic column name. (Only required if not providing beta+se, or neglog10p)
#' @param neglog10p_col A string. Default="NA". The -log10 p-value column name. (Only required if not providing beta+se, or stat)
#' @param n_bases An interger. Default=1e6. The distance between two significant SNPs, beyond which they are defined as in separate loci.
#' @param p_threshold A number. Default=5e-8. P-value threshold for statistical significance
#'
#' @examples
#' get_loci(gwas)
#'
#' @export
#'

get_loci = function(gwas,
                    snp_col       = "SNP",
                    chr_col       = "CHR",
                    pos_col       = "BP",
                    maf_col       = "MAF",
                    beta_col      = "BETA",
                    se_col        = "SE",
                    stat_col      = "NA",
                    neglog10p_col = "NA",
                    n_bases       = 1e6,
                    p_threshold   = 5e-8)  {

	## in case a tibble etc is passed...
	gwas = as.data.frame(gwas)

	## will use -log10 of the p-value in case any p-values were <5e-324 and are rounded to 0 in many software
	if (neglog10p_col != "NA")  gwas[,"P_neglog10"] = gwas[,neglog10p_col]
	if (stat_col != "NA")  gwas[,"stat"] = gwas[,stat_col]
	if (neglog10p_col == "NA" & stat_col == "NA")  {
		gwas[,"stat"] = gwas[,beta_col] / gwas[,se_col]
		gwas[,"P_neglog10"] = lukesRlib::get_p_neglog10( gwas[,"stat"] )
	}
	
	## determine "loci" 
	gwas_loci = NULL
	n_loci = 0
	
	## determine threshold in -log10 
	p_threshold_neglog10 = lukesRlib::get_p_neglog10(p_threshold, is_p=TRUE)
	
	# are any SNPs GWAS significant? i.e., -log10 of 5e-8
	if (any(gwas[,"P_neglog10"] > p_threshold_neglog10))
	{
	
		## creating GWAS hits file 
		gwas_loci = gwas[ gwas[,"P_neglog10"] > p_threshold_neglog10 , ]
		dim(gwas_loci)
		head(gwas_loci)
		
		## add empty "locus" column
		gwas_loci[,"locus"] = NA
		
		## order by BP and CHR
		gwas_loci = gwas_loci[ order(as.numeric(gwas_loci[,pos_col])) , ]
		gwas_loci = gwas_loci[ order(as.numeric(gwas_loci[,chr_col])) , ]
		
		## what chromosomes are in the GWAS results:
		chrs = unique(as.numeric(gwas_loci[,chr_col]))
		chrs = chrs[order(chrs)]
		chrs
		
		## start locus counting at 0
		locus = 0
		
		## for each chromosome, go down and create loci i.e. if a SNP is >X Mb from the previous SNP this is a new locus!
		for (i in chrs)
		{
			## new CHR, next locus:
			locus = locus + 1
		
			## start locus at 1... & get first position
			gwas_loci[gwas_loci[,chr_col] == i,"locus"][1] = locus
	
			## loop down GWAS results - if more than 1 SNP
			if (nrow(gwas_loci[gwas_loci[,chr_col] == i,]) > 1)  {
			for (j in 1:(nrow(gwas_loci[gwas_loci[,chr_col] == i,])-1) )
			{
				## by default put next SNP in this locus
				gwas_loci[gwas_loci[,chr_col] == i,"locus"][j+1] = locus
	
				## if distance to next SNP is >X Mb then next SNP is in new locus
				if (gwas_loci[gwas_loci[,chr_col] == i,pos_col][j+1] - gwas_loci[gwas_loci[,chr_col] == i,pos_col][j]  >  n_bases)
				{
					locus = locus + 1
					gwas_loci[gwas_loci[,chr_col] == i,"locus"][j+1] = locus
				}
			}
			}
		}
		
		#head(gwas_loci)
		#table(gwas_loci[,"locus"])
		
		n_loci = locus
		cat(paste0("N variants = ", nrow(gwas)), "\n")
		cat(paste0("N variants p<threshold = ", nrow(gwas_loci)), "\n")
		cat(paste0("N loci = ", n_loci, "\n"))
		
		######################################################
		## for each locus which is the "lead" SNP
		
		## create empty variable
		gwas_loci[,"lead"] = FALSE
		
		## for each locus:
		for (i in unique(gwas_loci[,"locus"]))
		{
			## which row(s) has the lowest p-value for this locus?
			r = which(gwas_loci[gwas_loci[,"locus"] == i,"P_neglog10"] == min(gwas_loci[gwas_loci[,"locus"] == i,"P_neglog10"]))
			#print(paste0("Locus ", i, " row(s) ", r))
			
			## some loci have 2 SNPs tied for significance and can have very high correlations! 
			## But which to choose as lead SNP?
			## use MAF (allele frequency), then BETA (effect size). 
			## If still tied, choose first.
			if (length(r)>1)
			{
				## get betas and mafs
				mafs = betas = rep(NA, length(r))
				for (j in 1:length(r))
				{
					mafs[j]  = gwas_loci[gwas_loci[,"locus"] == i,maf_col][r[j]]
					betas[j] = gwas_loci[gwas_loci[,"locus"] == i,beta_col][r[j]]
				}
				
				## only keep rows with highest maf or beta, in that order
				r = r[ which(mafs == max(mafs)) ]
				if (length(r)>1)  r = r[ which(betas == max(betas)) ]
				
				## still tied? just keep first one...
				if (length(r)>1)  r = r[ 1 ]
			}
			
			gwas_loci[gwas_loci[,"locus"] == i , "lead"][r] = TRUE
		}
		#table(gwas_loci[,"lead"])
	
	}
	
	gwas_loci
	
}



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
#' @param p_col Numeric. Default="P". The p-value column name.
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
                            snp_col       = "SNP",
                            chr_col       = "CHR",
                            pos_col       = "BP",
                            p_col         = "P",
                            build         = 37,
                            n_bases       = 1e5)  {
  
  # check input
  if (! build %in% c(37,38))  stop("Build has to be 37 or 38")
  gene_curated = snpsettest::gene.curated.GRCh37
  if (build == 38)  gene_curated = snpsettest::gene.curated.GRCh38
  
  # get SNPs actually in the genes
  gwas_map = gwas |> dplyr::mutate(id=!! rlang::sym(snp_col), chr=!! rlang::sym(chr_col), pos=!! rlang::sym(pos_col), pvalue = !! rlang::sym(p_col)) |> as.data.frame()
  gwas_mapped = snpsettest::map_snp_to_gene(gwas_map, gene_curated, extend_start=n_bases/1000, extend_end=n_bases/1000)
  gwas_mapped = gwas_mapped$map
  
  # merge gene symbols to IDs
  map = dplyr::left_join(gwas_mapped, gene_curated |> dplyr::select(gene.id, gene.name), by="gene.id")
  
  # get unique pos & distances
  map = map |> dplyr::mutate(dist_start = gene.start - pos, # positive integer = before start
                             dist_end   = gene.end - pos,   # positive integer = before end
                             dist       = NA)               # if variant is outside gene (both positive) update with distance
  map = map |> dplyr::mutate(dist = dplyr::case_when(
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
    
    mapped = data.frame(id, gene, dist)
    return(mapped)
    
  }
  
  # use `map` from {purrr} to make this quick and easy
  map2 = purrr::map(unique(map$id), \(id) get_nearest_gene1(id, map)) |> purrr::list_rbind()
  
  # merge original GWAS file with genes 
  gwas = gwas |> dplyr::rename(ID=!!snp_col)
  gwas_genes = dplyr::left_join(gwas, map2 |> dplyr::select(id, gene, dist), by=c("ID"="id"))
  gwas = gwas |> dplyr::rename(!! rlang::sym(snp_col) := ID)
  
  # return 
  return(gwas_genes)
}


