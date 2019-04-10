NBumiCoexpression <- function(counts, fit, gene_list=NULL, method=c("both", "on", "off")) {
	# Set up
	if (is.null(gene_list)) {
		gene_list <- names(fit$vals$tjs)
	}
	pd_gene <- matrix(-1, nrow=length(gene_list), ncol=ncol(counts));
	name_gene <- rep("", length(gene_list))
	for (i in 1:length(gene_list)) {
		gid <- which(names(fit$vals$tjs) == gene_list[i])
		if (length(gid) == 0) {next;}
		mu_is <- fit$vals$tjs[gid]*fit$vals$tis/fit$vals$total
                p_is <- (1+mu_is/fit$sizes[gid])^(-fit$sizes[gid]);
                pd_gene[i,] <- p_is;
                name_gene[i] <- gene_list[i];

	}
	if (sum(name_gene == "") > 0) {
		warning(paste("Warning:", sum(name_gene == ""), "genes not found, check your gene list is correct."));
		exclude <- which(name_gene == "");
		pd_gene <- pd_gene[-exclude,]
		name_gene <- name_gene[-exclude]
	}
	rownames(pd_gene) <- name_gene
	lib.size <- fit$vals$tis;

	Z_mat <- matrix(-1, ncol=length(pd_gene), nrow=length(pd_gene));
	for(i in 1:nrow(pd_gene)) {
		for (j in (i):nrow(pd_gene)) {
			p_g1 <- pd_gene[i,];
			p_g2 <- pd_gene[j,];	
			expr_g1 <- counts[rownames(counts)==rownames(pd_gene)[i],]
			expr_g2 <- counts[rownames(counts)==rownames(pd_gene)[j],]
	
			if (method == "off" | method=="both") {
				# Both zero 
				expect_both_zero <- p_g1*p_g2
				expect_both_err <- expect_both_zero*(1-expect_both_zero)

				obs_both_zero <- sum(expr_g1==0 & expr_g2==0)
				Z <- (obs_both_zero - sum(expect_both_zero)) / sqrt(sum(expect_both_err))
				#p_val <- pnorm(-abs(Z))*2
			}
			if (method == "on" | method=="both") {
				# both nonzero
				obs_both_nonzero <- sum(expr_g1!=0 & expr_g2!=0)
				expect_both_nonzero <- (1-p_g1)*(1-p_g2)
				expect_non_err <- expect_both_nonzero*(1-expect_both_nonzero)
				Z <- (obs_both_nonzero - sum(expect_both_nonzero)) / sqrt(sum(expect_non_err))
				#p_val <- pnorm(-abs(Z))*2
			}

			if (method == "both") {
				# either
				obs_either <- obs_both_zero+obs_both_nonzero
				expect_either <- expect_both_zero+expect_both_nonzero
				expect_err <- expect_either*(1-expect_either)
				Z <- (obs_either - sum(expect_either)) / sqrt(sum(expect_err))
				#p_val <- pnorm(-abs(Z))*2
			}
			Z_mat[i,j] <- Z_mat[j,i] <- Z;
		}
	}
	rownames(Z_mat) <- names(pd_gene);
	colnames(Z_mat) <- names(pd_gene);
	return(Z_mat);
}

