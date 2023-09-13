library(tidyverse)

# download LDL summary stats from GWAS catalog GCST90239658 :: Graham et al. 2021 Nature. Analysis of 1,320,016 Europeans
#  wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90239001-GCST90240000/GCST90239658/GCST90239658_buildGRCh37.tsv
gwas_example = vroom::vroom("GCST90239658_buildGRCh37.tsv")

gwas_example = gwas_example |>
	rename(
		SNP=variant_id, 
		CHR=chromosome,
		BP=base_pair_location,
		A1=other_allele,
		A2=effect_allele,
		MAF=effect_allele_frequency,
		BETA=beta, 
		SE=standard_error,
		P=p_value) |>
	select(-n, -N_studies)

# choose a few regions to include
gwas_example = gwas_example |>
	filter(
		(CHR==1  & BP>104275684 & BP<114275684) |  # CELSR2 locus
		(CHR==2  & BP>16029662 & BP<26029662) )  # APOE locus

# exclude missing & make sure data.frame
gwas_example = gwas_example |> na.omit() |> as.data.frame()

# save
usethis::use_data(gwas_example, overwrite = TRUE, compress = 'xz')

