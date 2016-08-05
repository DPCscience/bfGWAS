options(stringsAsFactors=F)
library(data.table)

#### Function to load SFBA results
LoadEMdata <- function(filename, header = FALSE){
    paramdata = fread(filename, sep = "\t", header = header)
    setnames(paramdata, c("ID", "chr", "bp", "maf", "func", "beta", "pi", "Zscore", "SE_beta", "LRT", "pval_LRT", "rank"))
    setkey(paramdata, ID)
    return(paramdata)
}

LoadEMhyp <- function(filename, header=FALSE){
    EMhyp_data <- read.table(filename, sep = "\t", header = header)
    n_type <- (dim(EMhyp_data)[2] - 3)/4
    temp_names <- c("Iter", "h", "loglike")
    for(i in 1:n_type){
        temp_names <- c(temp_names, paste(c("pi", "pi_se", "sigma2", "sigma2_se"), (i-1), sep="_"))
    }
    colnames(EMhyp_data)<- temp_names
    return(EMhyp_data)
}

getCI <- function(est, est_se, alpha){
    z_alpha <- -qnorm((1-alpha)/2)
    
    pi_low <- est - z_alpha * est_se
    if(pi_low < 0) pi_low = 0
    
    pi_up <- est + z_alpha * est_se
    return(c(pi_low, pi_up))
}


CItable <- function(v, n_type, alpha, funcgroup=c("nonsyn", "nonsyn","nonsyn", "others", "others", "others")){

    temp_table <- data.frame(matrix(NA, nrow = n_type, ncol = 6), 
        funcgroup= factor(funcgroup, levels = funcgroup) )

    colnames(temp_table)[1:6] <- c("pi", "pi_lcl", "pi_ucl", "v", "v_lcl", "v_ucl")
    
    for(i in 1:n_type){
        pi_hat <- v[paste("pi", (i-1), sep="_")]
        pi_se <- v[paste("pi_se", (i-1), sep="_")]

        sigma2_hat <- v[paste("sigma2", (i-1), sep="_")]
        sigma2_se <- v[paste("sigma2_se", (i-1), sep="_")]

        temp_table[i, 1:6] <- c(pi_hat, getCI(pi_hat, pi_se, alpha), sigma2_hat, getCI(sigma2_hat, sigma2_se, alpha))
    }
    return(temp_table)
}

####### Compare group-wise estimates #######
comp_groupEst <- function(est = group_pp, est_se = pp_se, n_vec, i = 1, conf = 0.95){

    est_1 = est[i]
    w = n_vec / sum(n_vec)
    est_all = sum(est * w )

    var_1 = est_se[i]^2
    var_all = sum(est_se^2 * w^2)

    if(est_all > 0){
    	if((est_1) >0 ){
    		log_ratio = unlist ( log(est_1 / est_all) )
    		log_ratio_se = sqrt(var_1 / est_1^2 + var_all / est_all^2)
    		pvalue = 2 * (1 - pnorm( abs(log_ratio), 0, log_ratio_se))
    		pvalue = unname( unlist(pvalue))

    		z = abs(qnorm((1 - conf)/2))
    		log_ratio_lcl = log_ratio - z * log_ratio_se
    		log_ratio_ucl = log_ratio + z * log_ratio_se
    		output = c(exp(unname( unlist(c(log_ratio_lcl, log_ratio, log_ratio_ucl))) ), pvalue)
    	}else {
    		print(paste("group", i, "has estimate", est_1))
    		output = rep(NA, 4)
    	}
	}else{
		print(paste("average estimate is", est_all))
		output = rep(NA, 4)
	}
    return(output)
}

####### Make plots

Plot_GroupwisePP <- function(hyp_table, pdfname="", size = 18, tit="Groupwise causal probabilities with 95% error bar", wid=6, priorprob = 1e-6){
    dodge <- position_dodge(width=0.9) 
    # set group estimates of these have causal probabilities equal to the priors as 0 
    hyp_table[hyp_table$pi == priorprob, 1:6] = 0
    p1 = ggplot(hyp_table, aes(y=pi, x = funcgroup, fill=funcgroup)) + geom_bar(position=dodge, stat="identity")+ guides (fill=FALSE) + 
         geom_errorbar(aes(ymax=pi_ucl, ymin=pi_lcl), width=0.3, position=dodge) + 
        labs(title = tit, x = "", y = "Causal Probability") + 
        theme_grey(base_size = size)
    pdf(pdfname, width = wid)
    print(p1)
    dev.off()

}

Plot_EffetVar <- function(hyp_table, pdfname="", size = 18, tit = "Effect-size variances with 95% error bar", wid=6 , priorprob = 1e-6){
    dodge <- position_dodge(width=0.9)
    hyp_table[hyp_table$pi == priorprob, 1:6] = 0
    p2 = ggplot(hyp_table, aes(y=v, x = funcgroup, fill=funcgroup)) + geom_bar(position=dodge, stat="identity") + guides(fill=FALSE) + 
         geom_errorbar(aes(ymax=v_ucl, ymin=v_lcl), width=0.3, position=dodge) + 
        labs(title = tit, x = "", y = "Effect-size Variance") + 
        theme_grey(base_size = size)
    pdf(pdfname, width = wid)
    print(p2)
    dev.off()
}

Plot_MarkerNum <- function(n_vec, anno_groups, pdfname = "", wid = 8, size = 18){
	
	p1 = ggplot(data.frame(count= as.vector(n_vec), funcgroup=factor(anno_groups, levels = anno_groups)), 
		aes(y=count, x = funcgroup, fill=funcgroup, label=count)) + 
		geom_bar(position=position_dodge(width=0.9), stat="identity") + 
		guides(fill=FALSE) + geom_text() +
	        labs(title = "# of Markers Per Group", x = "", y = "Count")+ 
	        theme_grey(base_size = size)

	pdf(pdfname, width = wid)       
	print(p1)
	dev.off()
}

Plot_Comp_PP <- function(comp_group_pp, anno_groups, pdfname = "", wid = 8, size = 18){
	dodge=position_dodge(width=0.9)
	p1 = ggplot(data.frame(comp_group_pp, funcgroup=factor(anno_groups, levels = anno_groups)), 
				aes(y=ratio, x = funcgroup, fill=funcgroup)) + 
		geom_bar(position=dodge, stat="identity") + 
		geom_errorbar(aes(ymax=ratio_ucl, ymin=ratio_lcl), width=0.3, 
							position=dodge) + 
		guides(fill=FALSE) + 
	        labs(title = "Compare groupwise PPs vs. average", x = "", y = expression(paste(pi[q], "/", pi[avg])) ) + 
	        theme_grey(base_size = size)

	pdf(pdfname, width = wid)
	print(p1)
	dev.off()
}

Plot_Comp_EffectVar <- function(comp_group_sigma2, anno_groups, pdfname = "", wid = 8, size = 18){
	dodge=position_dodge(width=0.9)

	p1 = ggplot(data.frame(comp_group_sigma2, funcgroup=factor(anno_groups, levels = anno_groups)), 
		aes(y=ratio, x = funcgroup, fill=funcgroup)) + 
	geom_bar(position=dodge, stat="identity") + 
	geom_errorbar(aes(ymax=ratio_ucl, ymin=ratio_lcl), width=0.3, 
				position=dodge) + 
	guides(fill=FALSE) + 
        labs(title = expression(paste("Compare group-wise ", sigma^2, "vs. average") ),
        		 x = "", y = expression(paste(sigma[q] ^2, "/", sigma[avg] ^2))) + 
        theme_grey(base_size = size)

    pdf(pdfname, width = wid)
	print(p1)
	dev.off()    
}




