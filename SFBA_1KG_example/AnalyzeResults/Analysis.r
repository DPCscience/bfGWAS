# rm(list=ls(all=TRUE))

####### Source Util Functions and set data directories
source("UtilFunc.r")

DataDir = "/net/fantasia/home/yjingj/GIT/SFBA_example/ExData/" # example data directory
OutDir = "/net/fantasia/home/yjingj/GIT/SFBA_example/AnalyzeResults/" # directory to save plots
ResultDir = "/net/fantasia/home/yjingj/GIT/SFBA_TEST" # result directory

######## Load and analyze GWAS results by SFBA

GWAS_Result = LoadEMdata( filename=paste(ResultDir, "/Eoutput/paramtemp5.txt", sep = "") )

#### Estiamted total number of causals
print(sum(GWAS_Result$pi))

#### Variants with association probabilities > 0.1
print(GWAS_Result[GWAS_Result$pi>0.1, ])

#### Variants with SVT pvalue < 5e-8
GWAS_Result[GWAS_Result$pval_LRT<5e-8, ]

#### True variants used to simulate this dataset
VCF_Data = fread(paste(DataDir, "vcfs/causalSNP.vcf", sep = ""), sep="\t", header=TRUE)
Causal_SNP_Result = GWAS_Result[VCF_Data$ID, ]
print(Causal_SNP_Result) # variants that were excluded from the analysis will show NA values

#### Check p-values for these true causal SNPs
print(Causal_SNP_Result[pval_LRT<5e-8, ])

# REMARK: you might want to check if the estimated variants with high association probability are in high LD with the true causal ones


######## Load results of hyper parameters ####### 
Hyp_Result <- LoadEMhyp( filename = paste(ResultDir, "/Eoutput/EM_result.txt", sep = "") )
anno_groups <- c("Coding", "UTR", "Promoter", "DHS", "Intronic", "Intergenic/Others")

CI_Table <- CItable(Hyp_Result[6, ], n_type=6, alpha=0.95, funcgroup = anno_groups)

Plot_GroupwisePP(hyp_table=CI_Table, pdfname=paste(OutDir, "PostProb_95CI.pdf", sep =""), 
						wid = 8, priorprob = 1e-6)

Plot_EffectVar(hyp_table=CI_Table, pdfname=paste(OutDir, "EffVar_95CI.pdf", sep =""), 
						wid = 8, priorprob = 1e-6)

######## Calculate number of variants per group ##########
n_type = 6
priorprob = 1e-6
n_vec = rep(NA, n_type)
for(k in 0:(n_type-1)){
	n_vec[k+1] = sum(GWAS_Result$func == k)
}

### Plot # of markers
Plot_MarkerNum(n_vec, anno_groups, pdfname = paste(OutDir, "MarkerCount.pdf", sep =""))


######## Compare hyper estimates to the genome-wide averages
pp_cols = (1:n_type) * 4
group_pp = as.vector( Hyp_Result[6, (pp_cols)])
no_asso_groups = (group_pp == priorprob)
group_pp[no_asso_groups] = 0
pp_se <- as.vector( Hyp_Result[6, (pp_cols + 1)])
pp_se[no_asso_groups] = 0

ncausal_vec = n_vec * group_pp
group_sigma2 <- as.vector( Hyp_Result[6, (pp_cols + 2)])
group_sigma2[no_asso_groups] = 0
sigma2_se <- as.vector( Hyp_Result[6, (pp_cols + 3)])
sigma2_se[no_asso_groups] = 0

### Construct comparison results
comp_group_sigma2 <-comp_group_pp <- data.frame(ratio_lcl = rep(NA, n_type), 
							ratio=rep(NA, n_type), 
							ratio_ucl = rep(NA, n_type), 
							pval = rep(NA, n_type))
for(i in 1:n_type){
	comp_group_pp[i, ] = comp_groupEst(est = group_pp, est_se = pp_se, n_vec = n_vec, i = i, conf = 0.95)
	comp_group_sigma2[i, ] = comp_groupEst(est = group_sigma2, est_se = sigma2_se, n_vec = ncausal_vec, i = i, conf = 0.95)
}

### Plot comparison results
Plot_Comp_PP(comp_group_pp, anno_groups, pdfname = paste(OutDir, "comp_group_pp.pdf", sep =""))
Plot_Comp_EffectVar(comp_group_sigma2, anno_groups, pdfname = paste(OutDir, "comp_group_sigma2.pdf", sep =""))

######## END ################










