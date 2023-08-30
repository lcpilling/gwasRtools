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
#' @param n_bases An interger. Default=5e5. The distance between two significant SNPs, beyond which they are defined as in separate loci.
#' @param p_threshold A number. Default=5e-8. P-value threshold for statistical significance
#' @param get_ld_indep Logical. Default=FALSE. Use Plink LD clumping to identify independent SNPs - see ieugwasr::ld_clump() docs
#' @param ld_pruning_r2 Numeric. Default=0.001. Pruning threshold for LD.
#' @param ld_clump_local Logical. Default=TRUE. If clumping using local installation (rather than IEU API) - see ieugwasr::ld_clump() docs
#' @param ld_plink_bin A string. Default="plink". Path to Plink v1.90 binary
#' @param ld_bfile A string. Default="/indy/data/ld/1kg_v3/EUR". Path to BIM/BED reference panel files
#'
#' @examples
#' get_loci(gwas)
#'
#' @export
#'

get_loci = function(gwas,
                    snp_col         = "SNP",
                    chr_col         = "CHR",
                    pos_col         = "BP",
                    maf_col         = "MAF",
                    beta_col        = "BETA",
                    se_col          = "SE",
                    stat_col        = "NA",
                    neglog10p_col   = "NA",
                    n_bases         = 5e5,
                    p_threshold     = 5e-8,
                    get_ld_indep    = FALSE,
                    ld_pruning_r2   = 0.001,
                    ld_clump_local  = TRUE,
                    ld_plink_bin    = "plink",
                    ld_bfile        = "/indy/data/ld/1kg_v3/EUR"
)  {

	## in case a tibble etc is passed...
	gwas = as.data.frame(gwas)

	## will use -log10 of the p-value in case any p-values were <5e-324 and are rounded to 0 in many software
	if (neglog10p_col != "NA")  gwas[,"P_neglog10"] = gwas[,neglog10p_col]
	if (stat_col != "NA")  gwas[,"stat"] = gwas[,stat_col]
	if (neglog10p_col == "NA" & stat_col == "NA")  {
		gwas[,"stat"] = gwas[,beta_col] / gwas[,se_col]
		gwas[,"P_neglog10"] = gwasRtools:::P_neglog10( gwas[,"stat"] )
	}
	
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
		dim(gwas_loci)
		head(gwas_loci)
		
		## sort by CHR and POS
		gwas_loci[,chr_col] = as.numeric(gwas_loci[,chr_col])
		gwas_loci = gwas_loci[order(gwas_loci[,chr_col], gwas_loci[,pos_col]), ]
		
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

		######################################################
		## for each locus which is the "lead" SNP
		
		## create empty variable
		gwas_loci[,"lead"] = FALSE
		
		## for each locus:
		for (i in unique(gwas_loci[,"locus"]))
		{
			## which row(s) has the lowest p-value for this locus?
			r = which(gwas_loci[gwas_loci[,"locus"] == i,"P_neglog10"] == max(gwas_loci[gwas_loci[,"locus"] == i,"P_neglog10"]))
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
	
		######################################################
		## for each CHR indentify independent SNPs using LD clumping
		
		if (get_ld_indep) {
			
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
			
			# get independent SNPs in each chromosome based on LD
			ld_indep = NULL
			#for(chr in unique(for_clumping[,"chr"])){
			#	print(chr)
				ld_indep = c(ld_indep, suppressMessages(
					ieugwasr::ld_clump(for_clumping,#[for_clumping[,"chr"]==chr,], 
					                   clump_kb  = n_bases/1000, 
					                   clump_r2  = ld_pruning_r2,
					                   plink_bin = ld_plink_bin, 
					                   bfile     = ld_bfile)$rsid))
			#}
			
			# add additional indep SNPs to loci object
			gwas_loci[,"ld_indep"] = FALSE
			gwas_loci[gwas_loci[,snp_col] %in% ld_indep,"ld_indep"] = TRUE
			
		}
		
		cat(paste0("\nN variants = ", nrow(gwas)), "\n")
		cat(paste0("N variants p<threshold = ", nrow(gwas_loci)), "\n")
		cat(paste0("N loci = ", n_loci, "\n"))
		if (get_ld_indep) cat(paste0("N LD indep variants = ", nrow(gwas_loci[gwas_loci$ld_indep==TRUE,]), "\n"))

	}
	
	######################################################
	## return final object
	
	gwas_loci = gwas_loci |> select(-stat, -P_neglog10)
	gwas_loci
	
}


