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
#' @param detect_headers Logical. Default=TRUE. Search input headers to see if BOLT-LMM, SAIGE, or REGENIE input (user therefore doesn't need to provide). If BOLT-LMM then automatically use p-value instead of SE.
#' @param snp_col A string. Default="SNP". The RSID/variantID column name.
#' @param chr_col A string. Default="CHR". The chromosome column name.
#' @param pos_col A string. Default="BP". The base pair/position column name.
#' @param maf_col A string. Default="MAF". The MAF/minor-allele-frequency column name.
#' @param beta_col A string. Default="BETA". The BETA column name.
#' @param se_col A string. Default="SE". The SE column name.
#' @param stat_col A string. Default="NA". The test statistic column name. (Only required if not providing beta+se, or neglog10p)
#' @param p_col A string. Default="NA". The p-value column name. (Only required if `use_pvalue==TRUE` i.e., using the provided p-value, e.g., P_BOLT_LMM)
#' @param neglog10p_col A string. Default="NA". The -log10 p-value column name. (Only required if not providing beta+se, or stat)
#' @param use_pvalue Logical. Default=FALSE. Use the provided p-value (in `p_col`) rather than computing from the test statistic? Useful for BOLT-LMM output
#' @param n_bases An interger. Default=5e5. The distance between two significant loci, beyond which they are defined as in separate loci.
#' @param p_threshold A number. Default=5e-8. P-value threshold for statistical significance
#' @param exclude_hla Logical. Default=FALSE. Treat HLA as one continuous locus
#' @param hla_pos A numeric vector of length 2. Default=c(25e6, 34e6). The HLA region on chromosome 6 to treat as one continuous locus if `exclude_hla==TRUE`
#' @param get_ld_indep Logical. Default=FALSE. Use Plink LD clumping to identify independent SNPs - see ieugwasr::ld_clump() docs
#' @param ld_pruning_r2 Numeric. Default=0.01. Pruning threshold for LD.
#' @param ld_clump_local Logical. Default=TRUE. If clumping using local installation (rather than IEU API) - see ieugwasr::ld_clump() docs
#' @param ld_plink_bin A string. Default="plink". Path to Plink v1.90 binary
#' @param ld_bfile A string. Default is to 5,000 random unrelated UK Biobank Europeans on my server :) needs a path to appropriate BIM/BED reference panel files on your server
#' @param verbose Logical. Default=FALSE. Be verbose
#'
#' @examples
#' gwas_loci = get_loci(gwas_example)
#'
#' head(gwas_loci)
#'
#' head(gwas_loci[ gwas_loci$lead==TRUE , ])
#'
#' @export
#'

get_loci = function(gwas,
                    detect_headers  = TRUE,
                    snp_col         = "SNP",
                    chr_col         = "CHR",
                    pos_col         = "BP",
                    maf_col         = "MAF",
                    beta_col        = "BETA",
                    se_col          = "SE",
                    stat_col        = "NA",
                    p_col           = "NA",
                    neglog10p_col   = "NA",
                    use_pvalue      = FALSE,
                    n_bases         = 5e5,
                    p_threshold     = 5e-8,
                    exclude_hla     = FALSE,
                    hla_pos         = c(25e6, 34e6),
                    get_ld_indep    = FALSE,
                    ld_pruning_r2   = 0.01,
                    ld_clump_local  = TRUE,
                    ld_plink_bin    = "plink",
                    ld_bfile        = "/indy/ukbiobank/data_14631/genetics/imputed_500k/5k_eur/ukb_imp_v3.qc_sub.5k_eur",
                    verbose         = FALSE
)  {

	cat(paste0("\nLocus size (bases) = ", n_bases, "\n"))
	cat(paste0("P-value threshold = ", p_threshold, "\n\n"))
	
	if (exclude_hla)  cat("\nHLA region will be treated as one continuous locus\n\n")
	
	## in case a tibble etc is passed...
	gwas = as.data.frame(gwas)
	if (verbose)  cat("GWAS has ", nrow(gwas), " variants\n")
	
	# check headers. Is this BOLT, SAIGE, or REGENIE output? If not, user needs to specify
	if (detect_headers)  {
		col_names = colnames(gwas)
		if ("P_BOLT_LMM" %in% col_names)  {
			cat("Detected BOLT-LMM input. Using default headers. Using p-value to get SEs. Disable with `detect_headers=FALSE`\n\n")
			snp_col  = "SNP"
			chr_col  = "CHR"
			pos_col  = "BP"
			maf_col  = "A1FREQ"
			beta_col = "BETA"
			se_col   = "SE"
			p_col    = "P_BOLT_LMM"
			use_pvalue = TRUE
		}
		if ("SNPID" %in% col_names & "AF_Allele2" %in% col_names)  {
			cat("Detected SAIGE input. Using default headers. Disable with `detect_headers=FALSE`\n\n")
			snp_col  = "SNPID"
			chr_col  = "CHR"
			pos_col  = "POS"
			maf_col  = "AF_Allele2"
			beta_col = "BETA"
			se_col   = "SE"
		}
		if ("ID" %in% col_names & "CHROM" %in% col_names & "GENPOS" %in% col_names)  {
			cat("Detected REGENIE input. Using default headers. Disable with `detect_headers=FALSE`\n\n")
			snp_col  = "ID"
			chr_col  = "CHROM"
			pos_col  = "GENPOS"
			maf_col  = "A1FREQ"
			beta_col = "BETA"
			se_col   = "SE"
		}
	}
	
	## will use -log10 of the p-value in case any p-values were <5e-324 and are rounded to 0 in many software
	if (verbose)  cat("Getting -log10 p-value\n")
	if (!use_pvalue & neglog10p_col != "NA")  gwas[,"P_neglog10"] = gwas[,neglog10p_col]
	if (!use_pvalue & stat_col != "NA")  gwas[,"stat_tmp"] = gwas[,stat_col]
	if (!use_pvalue & neglog10p_col == "NA" & stat_col == "NA")  {
		gwas[,"stat_tmp"] = gwas[,beta_col] / gwas[,se_col]
		gwas[,"P_neglog10"] = gwasRtools:::P_neglog10( gwas[,"stat_tmp"] )
	}
	if (use_pvalue) gwas[,"P_neglog10"] = gwasRtools:::P_neglog10( gwas[,p_col], is_p=TRUE )
	
	## determine "loci" 
	gwas_loci = NULL
	n_loci = 0
	
	## determine threshold in -log10 
	p_threshold_neglog10 = gwasRtools:::P_neglog10(p_threshold, is_p=TRUE)
	
	# are any SNPs GWAS significant? i.e., -log10 of 5e-8
	if (any(gwas[,"P_neglog10"] > p_threshold_neglog10))
	{
		
		## exclude if P_neglog10 is NA - implies problem with BETA or SE - suggest to user to provide the Z or P_neglog10 directly
		n_na = length(which( is.na(gwas[,"P_neglog10"]) ))
		if (n_na >= 1)  {
			cat(paste0("\n!!! Warning: excluding ", n_na, " variants where P_neglog10 is NA\n!!! this suggests a problem with the BETA or SE. Suggest providing the test statistic or P_neglog10 directly\n\n"))
			gwas = gwas[ !is.na(gwas[,"P_neglog10"]) , ]
		}
		
		## creating GWAS hits file 
		gwas_loci = gwas[ gwas[,"P_neglog10"] > p_threshold_neglog10 , ]
		
		if (verbose)  cat("Found ", nrow(gwas_loci), " with -log10 p-value below threshold\n")
		
		## add empty "locus" & lead column
		gwas_loci[,"locus"] = 0
		gwas_loci[,"lead"] = FALSE
		
		## order by BP and CHR
		gwas_loci[,chr_col] = as.numeric(gwas_loci[,chr_col])
		gwas_loci = gwas_loci[ order(as.numeric(gwas_loci[,pos_col])) , ]
		gwas_loci = gwas_loci[ order(as.numeric(gwas_loci[,chr_col])) , ]
		
		## what chromosomes are in the GWAS results:
		chrs = unique(as.numeric(gwas_loci[,chr_col]))
		chrs = chrs[order(chrs)]
		chrs
		
		## start locus counting at 0
		locus = 0
		
		if (verbose)  cat("Determining loci\n")
		
		# treat HLA region as one continuous locus?
		if (exclude_hla)  {
			gwas_loci_hla = gwas_loci[ gwas_loci[,chr_col]==6 & gwas_loci[,pos_col] > hla_pos[1] & gwas_loci[,pos_col] < hla_pos[2] ,]
			
			if (nrow( gwas_loci_hla )>0)  {
				gwas_loci = gwas_loci[ ! (gwas_loci[,chr_col]==6 & gwas_loci[,pos_col] > hla_pos[1] & gwas_loci[,pos_col] < hla_pos[2]) ,]
			}
		}
		
		## for each chromosome, go down and create loci 
		for (i in chrs)
		{
			
			## until all SNPs are assigned a locus, keep going
			while( any(gwas_loci[gwas_loci[,chr_col] == i,"locus"] == 0) ) 
			{
			
				## assign smallest p-value to new locus
				locus = locus + 1
				this_min_p = max(gwas_loci[gwas_loci[,chr_col] == i & gwas_loci[,"locus"] == 0,"P_neglog10"])
				gwas_loci[gwas_loci[,chr_col] == i & gwas_loci[,"locus"] == 0 & gwas_loci[,"P_neglog10"] == this_min_p, "locus"] = locus
				
				## if >1 pick based on MAF/BETA... otherwise just first one
				## which row(s) has the lowest p-value for this locus?
				r = which(gwas_loci[gwas_loci[,chr_col] == i & gwas_loci[,"locus"] == locus,"P_neglog10"] == this_min_p)
				if (length(r)>1)
				{
					## get betas and mafs
					mafs = betas = rep(NA, length(r))
					for (j in 1:length(r))
					{
						mafs[j]  = gwas_loci[gwas_loci[,chr_col] == i & gwas_loci[,"locus"] == locus, maf_col][r[j]]
						betas[j] = abs(gwas_loci[gwas_loci[,chr_col] == i & gwas_loci[,"locus"] == locus, beta_col][r[j]])
					}
					
					## only keep rows with highest maf or beta, in that order
					r = r[ which(mafs == max(mafs)) ]
					if (length(r)>1)  r = r[ which(betas == max(betas)) ]
					
					## still tied? just keep first one...
					if (length(r)>1)  r = r[ 1 ]
				}
				
				# define lead SNP
				gwas_loci[gwas_loci[,chr_col] == i & gwas_loci[,"locus"] == locus, "lead"][r] = TRUE
				
				## all variants +/- 500kb (n_bases) are in this new locus 
				lead_pos = gwas_loci[gwas_loci[,chr_col] == i & gwas_loci[,"locus"] == locus & gwas_loci[,"lead"] == TRUE, pos_col]
				gwas_loci[gwas_loci[,chr_col] == i & 
				          gwas_loci[,"locus"] == 0 &
				          gwas_loci[,pos_col] > lead_pos-n_bases &
				          gwas_loci[,pos_col] < lead_pos+n_bases 
				          ,"locus"] = locus
			
			} # end while
			
		} # end CHR loop
		
		# treat HLA region as one continuous locus?
		if (exclude_hla)  {
			
			if (nrow( gwas_loci_hla )>0)  {
				
				locus = locus + 1
				
				# make locus number all the same. get new "lead".
				gwas_loci_hla[,"locus"] = locus
				gwas_loci_hla[, "lead"] = FALSE
				gwas_loci_hla[ gwas_loci_hla[,"P_neglog10"] == max(gwas_loci_hla[,"P_neglog10"]), "lead"] = TRUE
				
				# add back to main object
				gwas_loci = rbind(gwas_loci, gwas_loci_hla)
				gwas_loci = gwas_loci[ order(as.numeric(gwas_loci[,pos_col])) , ]
				gwas_loci = gwas_loci[ order(as.numeric(gwas_loci[,chr_col])) , ]
				
			}
		}
		
		n_loci = locus
		if (verbose)  cat("Found ", n_loci, " loci\n")
		
		# reorder locus numbers 
		unique_loci = unique(gwas_loci[,"locus"])
		loci = gwas_loci[,"locus"]
		i = 1
		for (locus in unique_loci)  {
			gwas_loci[loci == locus, "locus"] = i
			i = i + 1
		}
		
		
		######################################################
		## for each CHR indentify independent SNPs using LD clumping
		## only if significant in CHRs 1:22
		
		if (get_ld_indep & any(gwas_loci[,chr_col] %in% 1:22) ) {
			
			# using the API?
			if (! ld_clump_local)  ld_bfile = ld_plink_bin = NULL
			
			# messages
			cat("** Performing LD clumping. Can take a few minutes\n")
			if (ld_clump_local)  cat("** Local Plink installation will be called -- output appears in the terminal screen\n")
			
			# get RSID, chr and p-value for ld_clump()
			for_clumping = data.frame(
				rsid = gwas_loci[,snp_col],
				chr  = gwas_loci[,chr_col],
				pval = 2*pnorm(-abs(gwas_loci[,beta_col]/gwas_loci[,se_col]))
			)
			
			# no X or Y for clumping step
			for_clumping = for_clumping[ for_clumping$chr %in% 1:22 , ]
			
			# treat HLA region as one continuous locus?
			if (exclude_hla)  {
				hla_rsids = gwas_loci[ gwas_loci[,chr_col]==6 & gwas_loci[,pos_col] > hla_pos[1] & gwas_loci[,pos_col] < hla_pos[2] , snp_col ]
				for_clumping = for_clumping[ ! for_clumping$rsid %in% hla_rsids , ]
			}
			
			# get independent SNPs based on LD
			ld_indep = NULL
			ld_indep = try({
				suppressMessages(
					ieugwasr::ld_clump(for_clumping,
					                   clump_kb  = n_bases/1000, 
					                   clump_r2  = ld_pruning_r2,
					                   plink_bin = ld_plink_bin, 
					                   bfile     = ld_bfile)$rsid
				)})
			
			# add additional indep SNPs to loci object
			gwas_loci[,"lead_dist"] = gwas_loci[,"lead"]
			gwas_loci[,"lead_ld"] = FALSE
			if (! is.null(ld_indep) )  gwas_loci[gwas_loci[,snp_col] %in% ld_indep,"lead_ld"] = TRUE
			
			# pick best "lead" SNP 
			#  if locus just has 1 (from lead_dist) use that 
			#  if LD clumping identified more that 1, use those 
			for (locus in unique(gwas_loci[,"locus"]))  {
				
				# if >1 variants identified from clumping, use distance-based variant 
				if (length(gwas_loci[ gwas_loci[,"locus"] == locus & gwas_loci[,"lead_ld"] == TRUE ,"lead_ld"]) > 1)  {
					
					gwas_loci[ gwas_loci[,"locus"] == locus ,"lead"] = gwas_loci[ gwas_loci[,"locus"] == locus ,"lead_ld"] 
					
				}
				
			}
			
		}
		
		cat(paste0("N variants = ", nrow(gwas), "\n"))
		cat(paste0("N variants p<threshold = ", nrow(gwas_loci), "\n"))
		cat(paste0("N loci = ", n_loci, "\n"))
		if (get_ld_indep) cat(paste0("N independent variants (LD R2 threshold ", ld_pruning_r2, ") = ", nrow(gwas_loci[gwas_loci$lead==TRUE,]), "\n\n"))

	}
	
	######################################################
	## return final object
	
	gwas_loci = gwas_loci[, ! colnames(gwas_loci) %in% c("stat_tmp","P_neglog10") ]
	gwas_loci
	
}
